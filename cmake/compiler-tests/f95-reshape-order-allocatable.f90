!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!
program test_reshape
   integer, dimension(4) :: x = [1, 2, 3, 4]
   integer, dimension(:), allocatable :: order

   allocate (order(2))
   order(:) = [2, 1]

   ! PGI <= 19.10 does not accept allocatables for the order parameter
   print *, reshape(x, shape=[2, 2], order=order)
end program
