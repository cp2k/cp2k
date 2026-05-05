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

!> @dir tblite/container
!> Proxy module for general interaction containers.
!>
!> Provides access to the abstract base class [[container_type]] for declaring
!> a new interaction. Container extending the base class can provide (semi-)classical
!> energy contributions, i.e. independent of the density and/or density-dependent
!> contributions, which are included in the minimization of the electronic energy.
!>
!> The [[container_type]] itself is immutable once constructed, restart data must
!> be associated with a [[container_cache]] object. This cache object can be
!> populated in an [[container_type:update]] procedure. The update procedure is
!> allowed to invalidate existing caches and reallocate them as needed.
!>
!> Additionally, to ease compilation of several containers in an actual calculator
!> the [[container_list]] type is provided, which itself extends the [[container_type]]
!> base class to provide an uniform interface.
!>
!> New containers are added to the list by using the [[container_list:push_back]]
!> procedure. Since the list owns the containers, the container object should be
!> provided as allocatable class polymorphic [[container_type]], i.e. by transfer
!> using the intrinsic `move_alloc` procedure or by automatic LHS allocation of
!> the polymorphic object.
module tblite_container
   use tblite_container_cache, only : container_cache
   use tblite_container_list, only : container_list
   use tblite_container_type, only : container_type
   implicit none
   public
end module tblite_container
