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

!> @file tblite/api/context.f90
!> Provides API exports for the #tblite_context handle.

!> API export for environment context setup
module tblite_api_context
   use, intrinsic :: iso_c_binding
   use mctc_env, only : error_type, fatal_error
   use tblite_context, only : context_type, context_logger, context_terminal
   use tblite_api_error, only : vp_error
   use tblite_api_version, only : namespace
   use tblite_api_utils, only : f_c_character
   implicit none
   private

   public :: vp_context
   public :: new_context_api, check_context_api, get_context_error_api, delete_context_api
   public :: set_context_logger_api, set_context_color_api, set_context_verbosity_api


   !> Void pointer to manage calculation context
   type :: vp_context
      !> Actual payload
      type(context_type) :: ptr
   end type vp_context


   abstract interface
      !> Interface for callbacks used in custom logger
      subroutine callback(verror, msg, len, udata) bind(C)
         import :: c_char, c_int, c_ptr, vp_error
         !> Error handle
         type(c_ptr), value :: verror
         !> Message payload to be displayed
         character(len=1, kind=c_char), intent(in) :: msg(*)
         !> Length of the message
         integer(c_int), value :: len
         !> Data pointer for callback
         type(c_ptr), value :: udata
      end subroutine callback
   end interface


   !> Custom logger for calculation context to display messages
   type, extends(context_logger) :: callback_logger
      !> Data pointer for callback
      type(c_ptr) :: udata = c_null_ptr
      !> Custom callback function to display messages
      procedure(callback), pointer, nopass :: callback => null()
      !> Error handle
      type(vp_error), pointer :: verror => null()
   contains
      !> Entry point for context instance to log message
      procedure :: message
   end type callback_logger


   logical, parameter :: debug = .false.


contains


!> Create new calculation context object
function new_context_api() &
      & result(vctx) &
      & bind(C, name=namespace//"new_context")
   type(vp_context), pointer :: ctx
   type(c_ptr) :: vctx

   if (debug) print '("[Info]", 1x, a)', "new_context"

   allocate(ctx)
   vctx = c_loc(ctx)

end function new_context_api


!> Create new calculation context object
function check_context_api(vctx) result(status) &
      & bind(C, name=namespace//"check_context")
   type(vp_context), pointer :: ctx
   type(c_ptr), value :: vctx
   integer(c_int) :: status

   if (debug) print '("[Info]", 1x, a)', "check_context"

   if (c_associated(vctx)) then
      call c_f_pointer(vctx, ctx)

      status = merge(1, 0, ctx%ptr%failed())
   else
      status = 2
   end if

end function check_context_api


!> Get error message from calculation environment
subroutine get_context_error_api(vctx, charptr, buffersize) &
      & bind(C, name=namespace//"get_context_error")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   character(kind=c_char), intent(inout) :: charptr(*)
   integer(c_int), intent(in), optional :: buffersize
   integer :: max_length
   type(error_type), allocatable :: error

   if (debug) print '("[Info]", 1x, a)', "get_context_error"

   if (c_associated(vctx)) then
      call c_f_pointer(vctx, ctx)

      if (present(buffersize)) then
         max_length = buffersize
      else
         max_length = huge(max_length) - 2
      end if

      call ctx%ptr%get_error(error)
      if (allocated(error)) then
         call f_c_character(error%message, charptr, max_length)
      end if
   end if

end subroutine get_context_error_api


subroutine set_context_color_api(vctx, color) &
      & bind(C, name=namespace//"set_context_color")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   integer(c_int), value :: color

   if (debug) print '("[Info]", 1x, a)', "set_context_color"

   if (c_associated(vctx)) then
      call c_f_pointer(vctx, ctx)

      ctx%ptr%terminal = context_terminal(color /= 0)
   end if
end subroutine set_context_color_api


subroutine set_context_verbosity_api(vctx, verbosity) &
      & bind(C, name=namespace//"set_context_verbosity")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   integer(c_int), value :: verbosity

   if (debug) print '("[Info]", 1x, a)', "set_context_verbosity"

   if (c_associated(vctx)) then
      call c_f_pointer(vctx, ctx)

      ctx%ptr%verbosity = verbosity
   end if
end subroutine set_context_verbosity_api


!> Delete context object
subroutine delete_context_api(vctx) &
      & bind(C, name=namespace//"delete_context")
   type(c_ptr), intent(inout) :: vctx
   type(vp_context), pointer :: ctx

   if (debug) print '("[Info]", 1x, a)', "delete_context"

   if (c_associated(vctx)) then
      call c_f_pointer(vctx, ctx)

      deallocate(ctx)
      vctx = c_null_ptr
   end if

end subroutine delete_context_api


!> Create a new custom logger for the calculation context
subroutine set_context_logger_api(vctx, vproc, vdata) &
      & bind(C, name=namespace//"set_context_logger")
   type(c_ptr), value :: vctx
   type(vp_context), pointer :: ctx
   type(c_funptr), value :: vproc
   procedure(callback), pointer :: fptr
   type(c_ptr), value :: vdata

   if (debug) print '("[Info]", 1x, a)', "set_context_logger"

   if (.not.c_associated(vctx)) return
   call c_f_pointer(vctx, ctx)

   if (c_associated(vproc)) then
      call c_f_procpointer(vproc, fptr)
      ctx%ptr%io = new_callback_logger(fptr, vdata)
   else
      if (allocated(ctx%ptr%io)) deallocate(ctx%ptr%io)
   end if

end subroutine set_context_logger_api


!> Create a new custom logger
function new_callback_logger(fptr, udata) result(self)
   !> Actual function used to display messages
   procedure(callback) :: fptr
   !> Data pointer used inside the callback function
   type(c_ptr), intent(in) :: udata
   !> New instance of the custom logger
   type(callback_logger) :: self

   self%callback => fptr
   self%udata = udata
end function new_callback_logger


!> Entry point for context type logger, transfers message from context to callback
subroutine message(self, msg, error)
   !> Instance of the custom logger with the actual logger callback function
   class(callback_logger), intent(inout) :: self
   !> Message payload from the calculation context
   character(len=*), intent(in) :: msg
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=1, kind=c_char) :: charptr(len(msg))
   integer(c_int) :: nchars

   if (.not.associated(self%verror)) then
      allocate(self%verror)
   end if

   charptr = transfer(msg, charptr)
   nchars = len(msg)
   call self%callback(c_loc(self%verror), charptr, nchars, self%udata)

   if (allocated(self%verror%ptr)) then
      call move_alloc(self%verror%ptr, error)
   end if
end subroutine message


end module tblite_api_context
