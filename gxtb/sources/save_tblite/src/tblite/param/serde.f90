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

!> @file tblite/param/serde.f90
!> Provides a base class for serializable and deserializable data

!> Definition of a parameter record with serde properties.
!> Each record knows how to serialize and deserialize itself.
module tblite_param_serde
   use mctc_env, only : error_type, fatal_error
   use tblite_toml, only : toml_table, toml_error, toml_parse, toml_dump
   implicit none
   private


   !> Serializable and deserializable parameter record
   type, public, abstract :: serde_record
   contains
      !> Reading of parametrization data
      generic :: load => load_from_file, load_from_unit, load_from_toml
      !> Read parametrization data from file
      procedure, private :: load_from_file
      !> Read parametrization data from formatted unit
      procedure, private :: load_from_unit
      !> Read parametrization data from TOML data structure
      procedure(load_from_toml), deferred :: load_from_toml
      !> Writing of parametrization data
      generic :: dump => dump_to_file, dump_to_unit, dump_to_toml
      !> Write parametrization data to file
      procedure, private :: dump_to_file
      !> Write parametrization data to formatted unit
      procedure, private :: dump_to_unit
      !> Write parametrization data to TOML data structure
      procedure(dump_to_toml), deferred :: dump_to_toml
   end type serde_record


   abstract interface
      !> Read parametrization data from TOML data structure
      subroutine load_from_toml(self, table, error)
         import :: serde_record, toml_table, error_type
         !> Instance of the parametrization data
         class(serde_record), intent(inout) :: self
         !> Data structure
         type(toml_table), intent(inout) :: table
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine load_from_toml
      !> Write parametrization data to TOML datastructure
      subroutine dump_to_toml(self, table, error)
         import :: serde_record, toml_table, error_type
         !> Instance of the parametrization data
         class(serde_record), intent(in) :: self
         !> Data structure
         type(toml_table), intent(inout) :: table
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine dump_to_toml
   end interface

contains


!> Read parametrization data from file
subroutine load_from_file(self, file, error)
   !> Instance of the parametrization data
   class(serde_record), intent(inout) :: self
   !> File name
   character(len=*), intent(in) :: file
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: unit
   logical :: exist

   inquire(file=file, exist=exist)
   if (.not.exist) then
     call fatal_error(error, "Could not find parameter file '"//file//"'")
     return
   end if

   open(file=file, newunit=unit)
   call self%load(unit, error)
   close(unit)
end subroutine load_from_file


!> Read parametrization data from file
subroutine load_from_unit(self, unit, error)
   !> Instance of the parametrization data
   class(serde_record), intent(inout) :: self
   !> File name
   integer, intent(in) :: unit
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_error), allocatable :: parse_error
   type(toml_table), allocatable :: table

   call toml_parse(table, unit, parse_error)

   if (allocated(parse_error)) then
      allocate(error)
      call move_alloc(parse_error%message, error%message)
      return
   end if

   call self%load(table, error)
   if (allocated(error)) return

end subroutine load_from_unit


!> Write parametrization data to file
subroutine dump_to_file(self, file, error)
   !> Instance of the parametrization data
   class(serde_record), intent(in) :: self
   !> File name
   character(len=*), intent(in) :: file
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: unit

   open(file=file, newunit=unit)
   call self%dump(unit, error)
   close(unit)
   if (allocated(error)) return

end subroutine dump_to_file


!> Write parametrization data to file
subroutine dump_to_unit(self, unit, error)
   !> Instance of the parametrization data
   class(serde_record), intent(in) :: self
   !> Formatted unit
   integer, intent(in) :: unit
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(toml_error), allocatable :: ser_error

   table = toml_table()

   call self%dump(table, error)

   call toml_dump(table, unit, ser_error)
   if (allocated(ser_error)) then
      call fatal_error(error, ser_error%message)
   end if

end subroutine dump_to_unit


end module tblite_param_serde
