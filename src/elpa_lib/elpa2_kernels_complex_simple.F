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

subroutine single_hh_trafo_complex(q, hh, nb, nq, ldq)

   implicit none

   integer, intent(in) :: nb, nq, ldq
   complex*16, intent(inout) :: q(ldq,*)
   complex*16, intent(in) :: hh(*)

   integer i
   complex*16 h1, tau1, x(nq)

   ! Just one Householder transformation

   x(1:nq) = q(1:nq,1)

   do i=2,nb
      x(1:nq) = x(1:nq) + q(1:nq,i)*conjg(hh(i))
   enddo

   tau1 = hh(1)
   x(1:nq) = x(1:nq)*(-tau1)

   q(1:nq,1) = q(1:nq,1) + x(1:nq)

   do i=2,nb
      q(1:nq,i) = q(1:nq,i) + x(1:nq)*hh(i)
   enddo

end

! --------------------------------------------------------------------------------------------------
