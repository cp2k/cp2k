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

!> @dir tblite/exchange
!> Contains the exchange related interactions.

!> @file tblite/exchange.f90
!> Reexports of the exchange interaction containers.

!> Proxy module for reexporting exchange-type interactions
module tblite_exchange
   use tblite_exchange_cache, only : exchange_cache
   use tblite_exchange_type, only : exchange_type
   use tblite_exchange_fock, only : exchange_fock, new_exchange_fock
   implicit none
   private

   public :: exchange_type, exchange_fock, new_exchange_fock, exchange_cache

end module tblite_exchange
