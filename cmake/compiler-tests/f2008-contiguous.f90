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

   ! test whether the compiler supports the CONTIGUOUS keyword
   integer, allocatable, target :: targ(:)
   integer, contiguous, pointer :: ptr(:)

   ! allocated data is always contiguous
   allocate (targ(10))
   ptr => targ

   ! IS_CONTIGUOUS was implemented in gcc-9 and is therefore not tested for yet
end program
