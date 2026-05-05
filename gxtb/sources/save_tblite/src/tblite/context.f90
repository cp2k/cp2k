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

!> @dir tblite/context
!> Contains the individual components of the context manager.

!> @file tblite/context.f90
!> Provides reexports for the context manager.

!> Proxy module for the calculation context related data types.
!>
!> Provides access to the #tblite_context_type::context_type to manage the calculation.
!> The context stores error messages generated while running which can be
!> queried using the #tblite_context_type::context_type::failed procedure.
!> To access the actual errors the messages can be removed using the
!> #tblite_context_type::context_type::get_error procedure.
!>
!> Errors can be pushed to a context by using #tblite_context_type::context_type::set_error,
!> A context containing any error will automatically be marked as failed. Once all
!> error messages are removed from again with #tblite_context_type::context_type::get_error
!> the failed status is removed from the context as well.
!>
!> The context provides also access to output by providing the
!> #tblite_context_type::context_type::message procedure. Using this procedure
!> for printouts ensures they are correctly handled in case of API usage.
!>
!> An output verbosity is available in the context as the member verbosity, all procedures
!> with access to the context will default to the verbosity of the context unless
!> the verbosity level is overwritten by an argument.
!>
!> To cutomize the output the #tblite_context_logger::context_logger abstract base class
!> is available. It provides a #tblite_context_logger::context_logger::message procedure
!> which is used by the context to create output. This type can be used to create callbacks
!> for customizing or redirecting the output of the library, e.g. when used from C or Python.
!>
!> ANSI escape sequences are available via the #tblite_context_terminal::context_terminal
!> class and which is instantiated in the #tblite_context_type::context_type::terminal
!> member. The terminal instance provides access to various terminal formatting codes to
!> produce the respective ANSI escape sequences if the terminal is created with color
!> support, otherwise return empty strings. Detecting whether the terminal supports
!> ANSI escape sequences is left as excerise to the consumer of the
!> #tblite_context_terminal::context_terminal class.
!>
!> For convenience operators are provided to add escape codes as well as concatenate
!> them directly with strings.
module tblite_context
   use tblite_context_logger, only : context_logger
   use tblite_context_terminal, only : context_terminal, escape, operator(+), operator(//)
   use tblite_context_type, only : context_type
   implicit none
   public
end module tblite_context
