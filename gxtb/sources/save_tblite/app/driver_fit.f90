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

!> Implementation of the driver entry points for parameter optimizations
module tblite_driver_fit
   use mctc_env, only : error_type, fatal_error, get_argument, wp
   use mctc_io, only : structure_type, read_structure, write_structure, &
      & filetype, get_filetype, to_symbol
   use mctc_io_resize, only : resize
   use tblite_cli, only : fit_config
   use tblite_context_type, only : context_type
   use tblite_output_format, only : format_string
   use tblite_param, only : param_record, count
   use tblite_fit_newuoa, only : newuoa
   use tblite_fit_settings, only : fit_settings
   use tblite_os, only : setenv, file_exists, delete_file
   implicit none
   private

   public :: main

   interface main
      module procedure :: fit_main
   end interface


contains


subroutine fit_main(config, error)
   type(fit_config), intent(in) :: config
   type(error_type), allocatable, intent(out) :: error

   type(fit_settings), target :: set
   type(context_type) :: ctx
   integer :: stat, npar
   real(wp), allocatable :: array(:)
   class(*), pointer :: handle
   real(wp), parameter :: conv = 1.0e-5_wp

   allocate(set%base)
   call set%base%load(config%param, error)
   if (allocated(error)) return

   call clear_metadata(set%base)

   call set%load(config%input, error)
   if (allocated(error)) return

   npar = count(set%mask)
   allocate(array(npar))
   if (set%relative) then
      array(:) = 1.0_wp
   else
      call set%base%dump(array, set%mask, error)
   end if
   if (allocated(error)) return

   call summary(ctx%unit, config, set)

   call setenv("TBLITE_EXE", config%prog, stat)
   call setenv("TBLITE_OUT", set%output, stat)
   call setenv("TBLITE_PAR", set%fitpar, stat)

   if (allocated(config%copy_input)) then
      call set%dump(config%copy_input, error)
      if (allocated(error)) return
      call ctx%message("[Info] Input settings dumped to '"//config%copy_input//"'")
   end if

   call write_param(set, array, error)
   if (allocated(error)) return

   if (config%dry_run) then
      call ctx%message("[Info] Initial parameters written to '"//set%fitpar//"'")
      return
   end if

   handle => set
   call newuoa(npar, 2*npar, array, set%trustr, conv, config%verbosity, set%max_iter*npar, &
      & eval, handle, error)
   if (allocated(error)) return

   call write_param(set, array, error)
   if (allocated(error)) return
   call ctx%message("[Info] Final parameters written to '"//set%fitpar//"'")
end subroutine fit_main


subroutine summary(unit, config, set)
   integer, intent(in) :: unit
   type(fit_config), intent(in) :: config
   type(fit_settings), intent(in) :: set

   write(unit, '(a, t30, a)') &
      & "Optimization method", set%method, &
      & "Max. function evaluations", format_string(count(set%mask)*set%max_iter, '(i0)'), &
      & "Number of parameters", format_string(count(set%mask), '(i0)'), &
      & "Initial trust radius", format_string(set%trustr, '(f8.6)'), &
      & "Initial parameter file", config%param, &
      & "Script command", "'"//set%script//"'"
end subroutine summary


function eval(n, x, h, error) result(f)
   integer, intent(in) :: n
   real(wp), intent(in) :: x(*)
   class(*), intent(in) :: h
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: f

   integer :: stat
   real(wp), allocatable :: actual(:), reference(:)

   f = huge(0.0_wp)

   select type(set => h)
   type is (fit_settings)
      call write_param(set, x(1:n), error)
      if (allocated(error)) return

      call delete_file(set%output)

      call execute_command_line(set%script, exitstat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Running '"//set%script//"' failed")
         return
      end if

      if (.not.file_exists(set%output)) then
         call fatal_error(error, "Output data '"//set%script//"' of script not found")
         return
      end if

      call read_data(set%output, actual, reference)

      f = sum((actual - reference)**2)
   end select

end function eval


subroutine read_data(file, actual, reference)
   character(len=*), intent(in) :: file
   real(wp), allocatable, intent(out) :: actual(:), reference(:)
   integer :: idata, stat, unit

   open(file=file, newunit=unit)
   call resize(actual)
   call resize(reference)
   idata = 0
   stat = 0
   do while(stat == 0)
      if (idata >= min(size(reference), size(actual))) then
         call resize(actual)
         call resize(reference)
      end if
      idata = idata + 1
      read(unit, *, iostat=stat) reference(idata), actual(idata)
      if (stat /= 0) idata = idata - 1
   end do
   close(unit)

   call resize(actual, idata)
   call resize(reference, idata)
end subroutine read_data


!> Make sure no stale meta data is present in parameter file
subroutine clear_metadata(param)
   type(param_record), intent(inout) :: param

   param%version = 0
   if (allocated(param%name)) deallocate(param%name)
   if (allocated(param%reference)) deallocate(param%reference)
end subroutine clear_metadata


subroutine write_param(set, array, error)
   type(fit_settings), intent(in) :: set
   real(wp), intent(in) :: array(:)
   type(error_type), allocatable, intent(out) :: error

   type(param_record) :: param
   real(wp), allocatable :: p(:)

   if (set%relative) then
      allocate(p(size(array)))
      call set%base%dump(p, set%mask, error)
      if (allocated(error)) return
      p(:) = array * p
      call param%load(p, set%base, set%mask, error)
      if (allocated(error)) return
   else
      call param%load(array, set%base, set%mask, error)
      if (allocated(error)) return
   end if

   call param%dump(set%fitpar, error)
   if (allocated(error)) return
end subroutine write_param


end module tblite_driver_fit
