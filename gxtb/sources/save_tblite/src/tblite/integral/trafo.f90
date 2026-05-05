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

!> @file tblite/integral/trafo.f90
!> Provides transformation from cartesian to spherical harmonic basis functions

!> Implemenation of transformations from cartesian to spherical harmonic basis functions
!>
!> Spherical harmonics use standard ordering, *i.e.* [-l, ..., 0, ..., l].
module tblite_integral_trafo
   use mctc_env, only : wp
   implicit none
   private

   public :: transform0, transform1, transform2
   public :: adjoint_transform0, adjoint_transform1, adjoint_transform2


   real(wp), parameter :: s3 = sqrt(3.0_wp)
   real(wp), parameter :: s3_4 = s3 * 0.5_wp
   real(wp), parameter :: dtrafo(5, 6) = reshape([&
      ! -2      -1       0       1       2
      & 0.0_wp, 0.0_wp, -0.5_wp, 0.0_wp,   s3_4, & ! xx
      & 0.0_wp, 0.0_wp, -0.5_wp, 0.0_wp,  -s3_4, & ! yy
      & 0.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp, & ! zz
      &     s3, 0.0_wp,  0.0_wp, 0.0_wp, 0.0_wp, & ! xy
      & 0.0_wp, 0.0_wp,  0.0_wp,     s3, 0.0_wp, & ! xz
      & 0.0_wp,     s3,  0.0_wp, 0.0_wp, 0.0_wp],& ! yz
      & shape(dtrafo))

   real(wp), parameter :: d32 = 3.0_wp/2.0_wp
   real(wp), parameter :: s3_8 = sqrt(3.0_wp/8.0_wp)
   real(wp), parameter :: s5_8 = sqrt(5.0_wp/8.0_wp)
   real(wp), parameter :: s6 = sqrt(6.0_wp)
   real(wp), parameter :: s15 = sqrt(15.0_wp)
   real(wp), parameter :: s15_4 = sqrt(15.0_wp/4.0_wp)
   real(wp), parameter :: s45 = sqrt(45.0_wp)
   real(wp), parameter :: s45_8 = sqrt(45.0_wp/8.0_wp)
   real(wp), parameter :: ftrafo(7, 10) = reshape([&
      ! -3       -2       -1       0        1         2         3
      &  0.0_wp,  0.0_wp,  0.0_wp, 0.0_wp,   -s3_8,   0.0_wp,     s5_8, & ! xxx
      &   -s5_8,  0.0_wp,   -s3_8, 0.0_wp,  0.0_wp,   0.0_wp,   0.0_wp, & ! yyy
      &  0.0_wp,  0.0_wp,  0.0_wp, 1.0_wp,  0.0_wp,   0.0_wp,   0.0_wp, & ! zzz
      &   s45_8,  0.0_wp,   -s3_8, 0.0_wp,  0.0_wp,   0.0_wp,   0.0_wp, & ! xxy
      &  0.0_wp,  0.0_wp,  0.0_wp,   -d32,  0.0_wp,    s15_4,   0.0_wp, & ! xxz
      &  0.0_wp,  0.0_wp,  0.0_wp, 0.0_wp,   -s3_8,   0.0_wp,   -s45_8, & ! xyy
      &  0.0_wp,  0.0_wp,  0.0_wp,   -d32,  0.0_wp,   -s15_4,   0.0_wp, & ! yyz
      &  0.0_wp,  0.0_wp,  0.0_wp, 0.0_wp,      s6,   0.0_wp,   0.0_wp, & ! xzz
      &  0.0_wp,  0.0_wp,      s6, 0.0_wp,  0.0_wp,   0.0_wp,   0.0_wp, & ! yzz
      &  0.0_wp,     s15,  0.0_wp, 0.0_wp,  0.0_wp,   0.0_wp,   0.0_wp],& ! xyz
      & shape(ftrafo))

   real(wp), parameter :: d38 = 3.0_wp/8.0_wp
   real(wp), parameter :: d34 = 3.0_wp/4.0_wp
   real(wp), parameter :: s5_16 = sqrt(5.0_wp/16.0_wp)
   real(wp), parameter :: s10 = sqrt(10.0_wp)
   real(wp), parameter :: s10_8 = sqrt(10.0_wp/8.0_wp)
   real(wp), parameter :: s35_4 = sqrt(35.0_wp/4.0_wp)
   real(wp), parameter :: s35_8 = sqrt(35.0_wp/8.0_wp)
   real(wp), parameter :: s35_64 = sqrt(35.0_wp/64.0_wp)
   real(wp), parameter :: s45_4 = sqrt(45.0_wp/4.0_wp)
   real(wp), parameter :: s315_8 = sqrt(315.0_wp/8.0_wp)
   real(wp), parameter :: s315_16 = sqrt(315.0_wp/16.0_wp)
   real(wp), parameter :: gtrafo(9, 15) = reshape([&
      !  -4     -3     -2     -1       0    1      2       3        4
      &  0._wp, 0._wp, 0._wp, 0._wp,   d38, 0._wp,-s5_16,  0._wp,  s35_64, & ! xxxx
      &  0._wp, 0._wp, 0._wp, 0._wp,   d38, 0._wp, s5_16,  0._wp,  s35_64, & ! yyyy
      &  0._wp, 0._wp, 0._wp, 0._wp, 1._wp, 0._wp, 0._wp,  0._wp,   0._wp, & ! zzzz
      &  s35_4, 0._wp,-s10_8, 0._wp, 0._wp, 0._wp, 0._wp,  0._wp,   0._wp, & ! xxxy
      &  0._wp, 0._wp, 0._wp, 0._wp, 0._wp,-s45_8, 0._wp,  s35_8,   0._wp, & ! xxxz
      & -s35_4, 0._wp,-s10_8, 0._wp, 0._wp, 0._wp, 0._wp,  0._wp,   0._wp, & ! xyyy
      &  0._wp,-s35_8, 0._wp,-s45_8, 0._wp, 0._wp, 0._wp,  0._wp,   0._wp, & ! yyyz
      &  0._wp, 0._wp, 0._wp, 0._wp, 0._wp,   s10, 0._wp,  0._wp,   0._wp, & ! xzzz
      &  0._wp, 0._wp, 0._wp,   s10, 0._wp, 0._wp, 0._wp,  0._wp,   0._wp, & ! yzzz
      &  0._wp, 0._wp, 0._wp, 0._wp,   d34, 0._wp, 0._wp,  0._wp,-s315_16, & ! xxyy
      &  0._wp, 0._wp, 0._wp, 0._wp,-3._wp, 0._wp, s45_4,  0._wp,   0._wp, & ! xxzz
      &  0._wp, 0._wp, 0._wp, 0._wp,-3._wp, 0._wp,-s45_4,  0._wp,   0._wp, & ! yyzz
      &  0._wp,s315_8, 0._wp,-s45_8, 0._wp, 0._wp, 0._wp,  0._wp,   0._wp, & ! xxyz
      &  0._wp, 0._wp, 0._wp, 0._wp, 0._wp,-s45_8, 0._wp,-s315_8,   0._wp, & ! xyyz
      &  0._wp, 0._wp,   s45, 0._wp, 0._wp, 0._wp, 0._wp,  0._wp,   0._wp],& ! xyzz
      &  shape(gtrafo))

contains


pure subroutine transform0(lj, li, cart, sphr)
   integer, intent(in) :: li
   integer, intent(in) :: lj
   real(wp), intent(in) :: cart(:, :)
   real(wp), intent(out) :: sphr(:, :)

   select case(li)
   case(0, 1)
      select case(lj)
      case(0, 1)
         sphr = cart
      case(2)
         !sphr = matmul(dtrafo, cart)
         sphr(3, :) = cart(3, :) - 0.5_wp * (cart(1, :) + cart(2, :))
         sphr(4, :) = s3 * cart(5, :)
         sphr(2, :) = s3 * cart(6, :)
         sphr(5, :) = s3_4 * (cart(1, :) - cart(2, :))
         sphr(1, :) = s3 * cart(4, :)
      case(3)
         sphr = matmul(ftrafo, cart)
      case(4)
         sphr = matmul(gtrafo, cart)
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(2)
      select case(lj)
      case(0, 1)
         !sphr = matmul(cart, transpose(dtrafo))
         sphr(:, 3) = cart(:, 3) - 0.5_wp * (cart(:, 1) + cart(:, 2))
         sphr(:, 4) = s3 * cart(:, 5)
         sphr(:, 2) = s3 * cart(:, 6)
         sphr(:, 5) = s3_4 * (cart(:, 1) - cart(:, 2))
         sphr(:, 1) = s3 * cart(:, 4)
      case(2)
         !sphr = matmul(dtrafo, matmul(cart, transpose(dtrafo)))
         sphr(3, 3) = cart(3, 3) &
            & - 0.5_wp * (cart(3, 1) + cart(3, 2) + cart(1, 3) + cart(2, 3)) &
            & + 0.25_wp * (cart(1, 1) + cart(1, 2) + cart(2, 1) + cart(2, 2))
         sphr([4, 2, 1], 3) = s3 * cart([5, 6, 4], 3) &
            & - s3_4 * (cart([5, 6, 4], 1) + cart([5, 6, 4], 2))
         sphr(5, 3) = s3_4 * (cart(1, 3) - cart(2, 3)) &
            & - s3 * 0.25_wp * (cart(1, 1) - cart(2, 1) + cart(1, 2) - cart(2, 2))
         sphr(3, 4) = s3 * cart(3, 5) - s3_4 * (cart(1, 5) + cart(2, 5))
         sphr([4, 2, 1], 4) = 3 * cart([5, 6, 4], 5)
         sphr(5, 4) = 1.5_wp * (cart(1, 5) - cart(2, 5))
         sphr(3, 2) = s3 * cart(3, 6) - s3_4 * (cart(1, 6) + cart(2, 6))
         sphr([4, 2, 1], 2) = 3 * cart([5, 6, 4], 6)
         sphr(5, 2) = 1.5_wp * (cart(1, 6) - cart(2, 6))
         sphr(3, 5) = s3_4 * (cart(3, 1) - cart(3, 2)) &
            & - s3 * 0.25_wp * (cart(1, 1) - cart(1, 2) + cart(2, 1) - cart(2, 2))
         sphr([4, 2, 1], 5) = 1.5_wp * (cart([5, 6, 4], 1) - cart([5, 6, 4], 2))
         sphr(5, 5) = 0.75_wp * (cart(1, 1) - cart(2, 1) - cart(1, 2) + cart(2, 2))
         sphr(3, 1) = s3 * cart(3, 4) - s3_4 * (cart(1, 4) + cart(2, 4))
         sphr([4, 2, 1], 1) = 3 * cart([5, 6, 4], 4)
         sphr(5, 1) = 1.5_wp * (cart(1, 4) - cart(2, 4))
      case(3)
         sphr = matmul(ftrafo, matmul(cart, transpose(dtrafo)))
      case(4)
         sphr = matmul(gtrafo, matmul(cart, transpose(dtrafo)))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(3)
      select case(lj)
      case(0, 1)
         sphr = matmul(cart, transpose(ftrafo))
      case(2)
         sphr = matmul(dtrafo, matmul(cart, transpose(ftrafo)))
      case(3)
         sphr = matmul(ftrafo, matmul(cart, transpose(ftrafo)))
      case(4)
         sphr = matmul(gtrafo, matmul(cart, transpose(ftrafo)))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(4)
      select case(lj)
      case(0, 1)
         sphr = matmul(cart, transpose(gtrafo))
      case(2)
         sphr = matmul(dtrafo, matmul(cart, transpose(gtrafo)))
      case(3)
         sphr = matmul(ftrafo, matmul(cart, transpose(gtrafo)))
      case(4)
         sphr = matmul(gtrafo, matmul(cart, transpose(gtrafo)))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case default
      error stop "[Fatal] Moments higher than g are not supported"
   end select

end subroutine transform0

pure subroutine transform1(lj, li, cart, sphr)
   integer, intent(in) :: li
   integer, intent(in) :: lj
   real(wp), intent(in) :: cart(:, :, :)
   real(wp), intent(out) :: sphr(:, :, :)
   integer :: k

   do k = 1, size(cart, 1)
      call transform0(lj, li, cart(k, :, :), sphr(k, :, :))
   end do
end subroutine transform1

pure subroutine transform2(lj, li, cart, sphr)
   integer, intent(in) :: li
   integer, intent(in) :: lj
   real(wp), intent(in) :: cart(:, :, :, :)
   real(wp), intent(out) :: sphr(:, :, :, :)
   integer :: k, l

   do l = 1, size(cart, 2)
      do k = 1, size(cart, 1)
         call transform0(lj, li, cart(k, l, :, :), sphr(k, l, :, :))
      end do
   end do
end subroutine transform2


!> Adjoint transformation from the spherical harmonic to the cartesian basis
!> for a shell pair block. Applies quantities which behave contravariant w.r.t.
!> the basis functions (i.e. MO expansion coefficients or the density matrix)
pure subroutine adjoint_transform0(lj, li, sphr, cart, bra, ket)
   !> Angular momentum of ket shell i
   integer, intent(in) :: li
   !> Angular momentum of bra shell j
   integer, intent(in) :: lj
   !> Spherical harmonic representation of the integral [bra j, ket i]
   real(wp), intent(in) :: sphr(:, :)
   !> Cartesian representation of the integral [bra j, ket i]
   real(wp), intent(out) :: cart(:, :)
   !> Flag for transformation of the bra dimension
   logical, intent(in) :: bra
   !> Flag for transformation of the ket dimension
   logical, intent(in) :: ket

   if (.not. bra .and. .not. ket) then
      cart = sphr
      return
   end if

   ! Transform only ket dimension
   if (.not. bra .and. ket) then
      select case(li)
      case(0, 1)
         cart = sphr
      case(2)
         cart = matmul(sphr, dtrafo)
      case(3)
         cart = matmul(sphr, ftrafo)
      case(4)
         cart = matmul(sphr, gtrafo)
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select
      return
   end if

   ! Transform only bra dimension
   if (bra .and. .not. ket) then
      select case(lj)
      case(0, 1)
         cart = sphr
      case(2)
         cart = matmul(transpose(dtrafo), sphr)
      case(3)
         cart = matmul(transpose(ftrafo), sphr)
      case(4)
         cart = matmul(transpose(gtrafo), sphr)
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select
      return
   end if

   ! Transform both dimensions
   select case(lj)
   case(0, 1)
      select case(li)
      case(0, 1)
         cart = sphr
      case(2)
         cart = matmul(sphr, dtrafo)
      case(3)
         cart = matmul(sphr, ftrafo)
      case(4)
         cart = matmul(sphr, gtrafo)
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(2)
      select case(li)
      case(0, 1)
         cart = matmul(transpose(dtrafo), sphr)
      case(2)
         cart = matmul(transpose(dtrafo), matmul(sphr, dtrafo))
      case(3)
         cart = matmul(transpose(dtrafo), matmul(sphr, ftrafo))
      case(4)
         cart = matmul(transpose(dtrafo), matmul(sphr, gtrafo))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(3)
      select case(li)
      case(0, 1)
         cart = matmul(transpose(ftrafo), sphr)
      case(2)
         cart = matmul(transpose(ftrafo), matmul(sphr, dtrafo))
      case(3)
         cart = matmul(transpose(ftrafo), matmul(sphr, ftrafo))
      case(4)
         cart = matmul(transpose(ftrafo), matmul(sphr, gtrafo))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case(4)
      select case(li)
      case(0, 1)
         cart = matmul(transpose(gtrafo), sphr)
      case(2)
         cart = matmul(transpose(gtrafo), matmul(sphr, dtrafo))
      case(3)
         cart = matmul(transpose(gtrafo), matmul(sphr, ftrafo))
      case(4)
         cart = matmul(transpose(gtrafo), matmul(sphr, gtrafo))
      case default
         error stop "[Fatal] Moments higher than g are not supported"
      end select

   case default
      error stop "[Fatal] Moments higher than g are not supported"
   end select

end subroutine adjoint_transform0

!> Adjoint transformation from the spherical harmonic to the cartesian basis
!> for a vector of shell pair block. Applies quantities which behave contravariant
!> w.r.t. the basis functions (i.e. MO expansion coefficients or the density matrix)
pure subroutine adjoint_transform1(lj, li, sphr, cart, bra, ket)
   !> Angular momentum of ket shell i
   integer, intent(in) :: li
   !> Angular momentum of bra shell j
   integer, intent(in) :: lj
   !> Spherical harmonic representation of the integral [:, bra j, ket i]
   real(wp), intent(in) :: sphr(:, :, :)
   !> Cartesian representation of the integral [:, bra j, ket i]
   real(wp), intent(out) :: cart(:, :, :)
   !> Flag for transformation of the bra dimension
   logical, intent(in) :: bra
   !> Flag for transformation of the ket dimension
   logical, intent(in) :: ket

   integer :: k

   do k = 1, size(sphr, 1)
      call adjoint_transform0(lj, li, sphr(k, :, :), cart(k, :, :), bra, ket)
   end do
end subroutine adjoint_transform1

!> Adjoint transformation from the spherical harmonic to the cartesian basis
!> for a matrix of shell pair block. Applies quantities which behave contravariant
!> w.r.t. the basis functions (i.e. MO expansion coefficients or the density matrix)
pure subroutine adjoint_transform2(lj, li, sphr, cart, bra, ket)
   !> Angular momentum of ket shell i
   integer, intent(in) :: li
   !> Angular momentum of bra shell j
   integer, intent(in) :: lj
   !> Spherical harmonic representation of the integral [:, :, bra j, ket i]
   real(wp), intent(in) :: sphr(:, :, :, :)
   !> Cartesian representation of the integral [:, :, bra j, ket i]
   real(wp), intent(out) :: cart(:, :, :, :)
   !> Flag for transformation of the bra dimension
   logical, intent(in) :: bra
   !> Flag for transformation of the ket dimension
   logical, intent(in) :: ket

   integer :: k, l

   do l = 1, size(sphr, 2)
      do k = 1, size(sphr, 1)
         call adjoint_transform0(lj, li, sphr(k, l, :, :), cart(k, l, :, :), bra, ket)
      end do
   end do
end subroutine adjoint_transform2


end module tblite_integral_trafo
