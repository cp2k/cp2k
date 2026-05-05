! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module dftd4_model_reference_d4
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   use dftd4_data, only : get_hardness, get_effective_charge
   use dftd4_model_utils, only : zeta
   implicit none
   private

   public :: get_nref, set_refcn, set_refgw
   public :: set_refq_eeq, set_refalpha_eeq
   public :: set_refq_gfn2, set_refalpha_gfn2
   public :: set_refq_eeqbc, set_refalpha_eeqbc

   interface get_nref
      module procedure :: get_nref_sym
      module procedure :: get_nref_num
   end interface get_nref

   interface set_refcn
      module procedure :: set_refcn_sym
      module procedure :: set_refcn_num
   end interface set_refcn

   interface set_refgw
      module procedure :: set_refgw_sym
      module procedure :: set_refgw_num
   end interface set_refgw

   interface set_refq_eeq
      module procedure :: set_refq_eeq_sym
      module procedure :: set_refq_eeq_num
   end interface set_refq_eeq

   interface set_refalpha_eeq
      module procedure :: set_refalpha_eeq_sym
      module procedure :: set_refalpha_eeq_num
   end interface set_refalpha_eeq

   interface set_refq_gfn2
      module procedure :: set_refq_gfn2_sym
      module procedure :: set_refq_gfn2_num
   end interface set_refq_gfn2

   interface set_refalpha_gfn2
      module procedure :: set_refalpha_gfn2_sym
      module procedure :: set_refalpha_gfn2_num
   end interface set_refalpha_gfn2

   interface set_refq_eeqbc
      module procedure :: set_refq_eeqbc_sym
      module procedure :: set_refq_eeqbc_num
   end interface set_refq_eeqbc

   interface set_refalpha_eeqbc
      module procedure :: set_refalpha_eeqbc_sym
      module procedure :: set_refalpha_eeqbc_num
   end interface set_refalpha_eeqbc

   integer, parameter :: max_elem = 118

   integer, dimension(max_elem)      :: refn ! for D4
   real(wp),dimension(7,max_elem)    :: refq ! GFN2-xTB charges
   real(wp),dimension(7,max_elem)    :: refh ! GFN2-xTB charges
   real(wp),dimension(7,max_elem)    :: dftq,pbcq,gffq,clsq,eeqbcq !solq
   real(wp),dimension(7,max_elem)    :: dfth,pbch,gffh,clsh,eeqbch !solh
   real(wp),dimension(7,max_elem)    :: hcount
   real(wp),dimension(7,max_elem)    :: ascale
   real(wp),dimension(7,max_elem)    :: refcovcn
   real(wp),dimension(7,max_elem)    :: refcn
   integer, dimension(7,max_elem)    :: refsys
   real(wp),dimension(23,7,max_elem) :: alphaiw
   real(wp),dimension(17)       :: secq
   real(wp),dimension(17)       :: sscale
   real(wp),dimension(17)       :: seccn
   real(wp),dimension(17)       :: seccnd3
   real(wp),dimension(23,17)    :: secaiw

   include 'd4.inc'

contains


!> Get number of references for a given element symbol
elemental function get_nref_sym(sym) result(n)

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> Number of references
   integer :: n

   n = get_nref(to_number(sym))

end function get_nref_sym


!> Get number of references for a given atomic number
elemental function get_nref_num(num) result(n)

   !> Atomic number
   integer, intent(in) :: num

   !> Number of references
   integer :: n

   if (num > 0 .and. num <= size(refn)) then
      n = refn(num)
   else
      n = 0
   end if

end function get_nref_num


!> Set the reference coordination numbers for an element symbol
pure subroutine set_refcn_sym(cn, sym)

   !> Reference coordination number
   real(wp), intent(out) :: cn(:)

   !> Element symbol
   character(len=*), intent(in) :: sym

   call set_refcn(cn, to_number(sym))

end subroutine set_refcn_sym


!> Set the reference coordination numbers for an atomic number
pure subroutine set_refcn_num(cn, num)

   !> Reference coordination number
   real(wp), intent(out) :: cn(:)

   !> Atomic number
   integer, intent(in) :: num

   integer :: ref

   cn(:) = 0.0_wp
   if (num > 0 .and. num <= size(refn)) then
      ref = get_nref(num)
      cn(:ref) = refcovcn(:ref, num)
   end if

end subroutine set_refcn_num


!> Set the number of gaussian weights for an element symbol
pure subroutine set_refgw_sym(ngw, sym)

   !> Number of gaussian weights
   integer, intent(out) :: ngw(:)

   !> Element symbol
   character(len=*), intent(in) :: sym

   call set_refgw(ngw, to_number(sym))

end subroutine set_refgw_sym


!> Set the number of gaussian weights for an atomic number
pure subroutine set_refgw_num(ngw, num)

   !> Number of gaussian weights
   integer, intent(out) :: ngw(:)

   !> Atomic number
   integer, intent(in) :: num

   integer, parameter :: max_cn = 19
   integer :: icn, ir, ref
   integer :: cnc(0:max_cn)

   ngw(:) = 1
   if (num > 0 .and. num <= size(refn)) then
      ref = get_nref(num)
      cnc(:) = [1, spread(0, 1, max_cn)]
      do ir = 1, ref
         icn = min(nint(refcn(ir, num)), max_cn)
         cnc(icn) = cnc(icn) + 1
      end do
      do ir = 1, ref
         icn = cnc(min(nint(refcn(ir, num)), max_cn))
         ngw(ir) = icn*(icn+1)/2
      end do
   end if

end subroutine set_refgw_num


!> Set the reference partial charges for an element symbol
pure subroutine set_refq_eeq_sym(q, sym)

   !> Reference partial charge
   real(wp), intent(out) :: q(:)

   !> Element symbol
   character(len=*), intent(in) :: sym

   call set_refq_eeq(q, to_number(sym))

end subroutine set_refq_eeq_sym


!> Set the reference partial charges for an atomic number
pure subroutine set_refq_eeq_num(q, num)

   !> Reference partial charge
   real(wp), intent(out) :: q(:)

   !> Atomic number
   integer, intent(in) :: num

   integer :: ref

   q(:) = 0.0_wp
   if (num > 0 .and. num <= size(refn)) then
      ref = get_nref(num)
      q(:ref) = clsq(:ref, num)
   end if

end subroutine set_refq_eeq_num


!> Set the reference partial charges for an element symbol
pure subroutine set_refq_gfn2_sym(q, sym)

   !> Reference partial charge
   real(wp), intent(out) :: q(:)

   !> Element symbol
   character(len=*), intent(in) :: sym

   call set_refq_gfn2(q, to_number(sym))

end subroutine set_refq_gfn2_sym


!> Set the reference partial charges for an atomic number
pure subroutine set_refq_gfn2_num(q, num)

   !> Reference partial charge
   real(wp), intent(out) :: q(:)

   !> Atomic number
   integer, intent(in) :: num

   integer :: ref

   q(:) = 0.0_wp
   if (num > 0 .and. num <= size(refn)) then
      ref = get_nref(num)
      q(:ref) = refq(:ref, num)
   end if

end subroutine set_refq_gfn2_num


!> Set the reference partial charges for an element symbol
pure subroutine set_refq_eeqbc_sym(q, sym)

   !> Reference partial charge
   real(wp), intent(out) :: q(:)

   !> Element symbol
   character(len=*), intent(in) :: sym

   call set_refq_eeqbc(q, to_number(sym))

end subroutine set_refq_eeqbc_sym


!> Set the reference partial charges for an atomic number
pure subroutine set_refq_eeqbc_num(q, num)

   !> Reference partial charge
   real(wp), intent(out) :: q(:)

   !> Atomic number
   integer, intent(in) :: num

   integer :: ref

   q(:) = 0.0_wp
   if (num > 0 .and. num <= size(refn)) then
      ref = get_nref(num)
      q(:ref) = eeqbcq(:ref, num)
   end if

end subroutine set_refq_eeqbc_num


!> Set the reference polarizibility for an element symbol
pure subroutine set_refalpha_eeq_sym(alpha, ga, gc, sym)

   !> Reference polarizibility
   real(wp), intent(out) :: alpha(:, :)

   !> Maximum charge scaling height
   real(wp), intent(in) :: ga

   !> Charge scaling steepness
   real(wp), intent(in) :: gc

   !> Element symbol
   character(len=*), intent(in) :: sym

   call set_refalpha_eeq(alpha, ga, gc, to_number(sym))

end subroutine set_refalpha_eeq_sym


!> Set the reference polarizibility for an atomic number
pure subroutine set_refalpha_eeq_num(alpha, ga, gc, num)

   !> Reference polarizibility
   real(wp), intent(out) :: alpha(:, :)

   !> Maximum charge scaling height
   real(wp), intent(in) :: ga

   !> Charge scaling steepness
   real(wp), intent(in) :: gc

   !> Atomic number
   integer, intent(in) :: num

   integer :: ref
   integer :: ir, is
   real(wp) :: iz
   real(wp) :: aiw(23)

   alpha(:, :) = 0.0_wp
   if (num > 0 .and. num <= size(refn)) then
      ref = get_nref(num)
      do ir = 1, ref
         is = refsys(ir, num)
         if (abs(is) < 1e-12_wp) cycle

         iz = get_effective_charge(is)
         aiw = sscale(is)*secaiw(:, is) &
            &    * zeta(ga, get_hardness(is)*gc, iz, clsh(ir, num)+iz)
         alpha(:, ir) = max(ascale(ir, num)*(alphaiw(:, ir, num) &
            &            - hcount(ir, num)*aiw), 0.0_wp)
      end do
   end if

end subroutine set_refalpha_eeq_num


!> Set the reference polarizibility for an element symbol
pure subroutine set_refalpha_gfn2_sym(alpha, ga, gc, sym)

   !> Reference polarizibility
   real(wp), intent(out) :: alpha(:, :)

   !> Maximum charge scaling height
   real(wp), intent(in) :: ga

   !> Charge scaling steepness
   real(wp), intent(in) :: gc

   !> Element symbol
   character(len=*), intent(in) :: sym

   call set_refalpha_gfn2(alpha, ga, gc, to_number(sym))

end subroutine set_refalpha_gfn2_sym


!> Set the reference polarizibility for an atomic number
pure subroutine set_refalpha_gfn2_num(alpha, ga, gc, num)

   !> Reference polarizibility
   real(wp), intent(out) :: alpha(:, :)

   !> Maximum charge scaling height
   real(wp), intent(in) :: ga

   !> Charge scaling steepness
   real(wp), intent(in) :: gc

   !> Atomic number
   integer, intent(in) :: num

   integer :: ref
   integer :: ir, is
   real(wp) :: iz
   real(wp) :: aiw(23)

   alpha(:, :) = 0.0_wp
   if (num > 0 .and. num <= size(refn)) then
      ref = get_nref(num)
      do ir = 1, ref
         is = refsys(ir, num)
         if (abs(is) < 1e-12_wp) cycle

         iz = get_effective_charge(is)
         aiw = sscale(is)*secaiw(:, is) &
            &    * zeta(ga, get_hardness(is)*gc, iz, refh(ir, num)+iz)
         alpha(:, ir) = max(ascale(ir, num)*(alphaiw(:, ir, num) &
            &            - hcount(ir, num)*aiw), 0.0_wp)
      end do
   end if

end subroutine set_refalpha_gfn2_num


!> Set the reference polarizibility for an element symbol
pure subroutine set_refalpha_eeqbc_sym(alpha, ga, gc, sym)

   !> Reference polarizibility
   real(wp), intent(out) :: alpha(:, :)

   !> Maximum charge scaling height
   real(wp), intent(in) :: ga

   !> Charge scaling steepness
   real(wp), intent(in) :: gc

   !> Element symbol
   character(len=*), intent(in) :: sym

   call set_refalpha_eeqbc(alpha, ga, gc, to_number(sym))

end subroutine set_refalpha_eeqbc_sym


!> Set the reference polarizibility for an atomic number
pure subroutine set_refalpha_eeqbc_num(alpha, ga, gc, num)

   !> Reference polarizibility
   real(wp), intent(out) :: alpha(:, :)

   !> Maximum charge scaling height
   real(wp), intent(in) :: ga

   !> Charge scaling steepness
   real(wp), intent(in) :: gc

   !> Atomic number
   integer, intent(in) :: num

   integer :: ref
   integer :: ir, is
   real(wp) :: iz
   real(wp) :: aiw(23)

   alpha(:, :) = 0.0_wp
   if (num > 0 .and. num <= size(refn)) then
      ref = get_nref(num)
      do ir = 1, ref
         is = refsys(ir, num)
         if (abs(is) < 1e-12_wp) cycle

         iz = get_effective_charge(is)
         aiw = sscale(is)*secaiw(:, is) &
            &    * zeta(ga, get_hardness(is)*gc, iz, eeqbch(ir, num)+iz)
         alpha(:, ir) = max(ascale(ir, num)*(alphaiw(:, ir, num) &
            &            - hcount(ir, num)*aiw), 0.0_wp)
      end do
   end if

end subroutine set_refalpha_eeqbc_num

end module dftd4_model_reference_d4
