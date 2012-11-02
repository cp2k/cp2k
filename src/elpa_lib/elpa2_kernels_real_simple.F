! --------------------------------------------------------------------------------------------------
!
! This file contains the compute intensive kernels for the Householder transformations.
!
! This is the small and simple version (no hand unrolling of loops etc.) but for some
! compilers this performs better than a sophisticated version with transformed and unrolled loops.
!
! It should be compiled with the highest possible optimization level.
! 
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".
!
! --------------------------------------------------------------------------------------------------

subroutine double_hh_trafo_2hv(q, hh, nb, nq, ldq, ldh)

   implicit none

   integer, intent(in) :: nb, nq, ldq, ldh
   real*8, intent(inout) :: q(ldq,*)
   real*8, intent(in) :: hh(ldh,*)

   real*8 s, h1, h2, tau1, tau2, x(nq), y(nq)
   integer i

   ! Calculate dot product of the two Householder vectors

   s = hh(2,2)*1
   do i=3,nb
      s = s+hh(i,2)*hh(i-1,1)
   enddo

   ! Do the Householder transformations

   x(1:nq) = q(1:nq,2)

   y(1:nq) = q(1:nq,1) + q(1:nq,2)*hh(2,2)

   do i=3,nb
      h1 = hh(i-1,1)
      h2 = hh(i,2)
      x(1:nq) = x(1:nq) + q(1:nq,i)*h1
      y(1:nq) = y(1:nq) + q(1:nq,i)*h2
   enddo

   x(1:nq) = x(1:nq) + q(1:nq,nb+1)*hh(nb,1)

   tau1 = hh(1,1)
   tau2 = hh(1,2)

   h1 = -tau1
   x(1:nq) = x(1:nq)*h1
   h1 = -tau2
   h2 = -tau2*s
   y(1:nq) = y(1:nq)*h1 + x(1:nq)*h2

   q(1:nq,1) = q(1:nq,1) + y(1:nq)
   q(1:nq,2) = q(1:nq,2) + x(1:nq) + y(1:nq)*hh(2,2)

   do i=3,nb
      h1 = hh(i-1,1)
      h2 = hh(i,2)
      q(1:nq,i) = q(1:nq,i) + x(1:nq)*h1 + y(1:nq)*h2
   enddo

   q(1:nq,nb+1) = q(1:nq,nb+1) + x(1:nq)*hh(nb,1)

end

! --------------------------------------------------------------------------------------------------
