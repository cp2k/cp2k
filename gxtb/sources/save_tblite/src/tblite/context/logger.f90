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

!> @file tblite/context/logger.f90
!> Provides an abstract base class for implementing a logger

!> Logger to display strings
module tblite_context_logger
   use mctc_env, only : error_type
   implicit none
   private


   !> Base class defining the logger interface
   type, public, abstract :: context_logger
   contains
      !> Entry point for displaying a string in the logger
      procedure(message), deferred :: message
   end type context_logger


   abstract interface
      !> Entry point for displaying a string in the logger
      subroutine message(self, msg, error)
         import :: context_logger, error_type
         !> Instance of the logger
         class(context_logger), intent(inout) :: self
         !> String to display
         character(len=*), intent(in) :: msg
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine message
   end interface


end module tblite_context_logger
