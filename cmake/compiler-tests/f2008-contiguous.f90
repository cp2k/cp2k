!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!

program main
   implicit none

   ! test whether the compiler supports the CONTIGUOUS keyword and the
   ! IS_CONTIGUOUS intrinsic, both part of the Fortran 2008 standard
   integer, allocatable, target :: targ(:)
   integer, contiguous, pointer :: ptr(:)
   integer, target :: targ2d(10, 10)

   ! allocated data is always contiguous
   allocate (targ(10))
   ptr => targ

   ! IS_CONTIGUOUS has been available since GCC 9 and is required by CP2K
   if (.not. IS_CONTIGUOUS(targ)) error stop 1
   if (.not. IS_CONTIGUOUS(ptr)) error stop 1
   ! a strided array section is not contiguous
   if (IS_CONTIGUOUS(targ2d(1:10:2, 1))) error stop 2
end program
