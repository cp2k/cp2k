! --------------------------------------------------------------------------------------------------
!
! This file contains the compute intensive kernels for the Householder transformations.
!
! *** Special IBM BlueGene/P version with BlueGene assembler instructions in Fortran ***
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

   real*8 s
   integer i

   ! Safety only:

   if(mod(ldq,4) /= 0) STOP 'double_hh_trafo: ldq not divisible by 4!'
   if(mod(loc(q),16) /= 0) STOP 'Q unaligned!'

   ! Calculate dot product of the two Householder vectors

   s = hh(2,2)*1
   do i=3,nb
      s = s+hh(i,2)*hh(i-1,1)
   enddo

   do i=1,nq-16,20
      call hh_trafo_kernel_10_bg(q(i   ,1), hh, nb, ldq, ldh, s)
      call hh_trafo_kernel_10_bg(q(i+10,1), hh, nb, ldq, ldh, s)
   enddo

   ! i > nq-16 now, i.e. at most 16 rows remain

   if(nq-i+1 > 12) then
      call hh_trafo_kernel_8_bg(q(i  ,1), hh, nb, ldq, ldh, s)
      call hh_trafo_kernel_8_bg(q(i+8,1), hh, nb, ldq, ldh, s)
   else if(nq-i+1 > 8) then
     call hh_trafo_kernel_8_bg(q(i  ,1), hh, nb, ldq, ldh, s)
     call hh_trafo_kernel_4_bg(q(i+8,1), hh, nb, ldq, ldh, s)
   else if(nq-i+1 > 4) then
     call hh_trafo_kernel_8_bg(q(i  ,1), hh, nb, ldq, ldh, s)
   else if(nq-i+1 > 0) then
     call hh_trafo_kernel_4_bg(q(i  ,1), hh, nb, ldq, ldh, s)
   endif

end

! --------------------------------------------------------------------------------------------------
! The following kernels perform the Householder transformation on Q for 10/8/4 rows.
! Please note that Q is declared complex*16 here.
! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_kernel_10_bg(q, hh, nb, ldq, ldh, s)


   implicit none

#if defined __parallel
   include 'mpif.h'
#endif   

   integer, intent(in) :: nb, ldq, ldh
   complex*16, intent(inout) :: q(ldq/2,*)
   real*8, intent(in) :: hh(ldh,*), s

   complex*16 x1, x2, x3, x4, x5, y1, y2, y3, y4, y5, q1, q2, q3, q4, q5, p1, p2, p3, p4, p5
   real*8 h1, h2
   integer i

!   complex*16 loadfp, fxcpmadd, fxpmul, fpadd, a, b
!   real*8 x
!   loadfp(a) = a
!   fxcpmadd(a,b,x) = a + b*x
!   fxpmul(a,x) = a*x
!   fpadd(a,b) = a+b
!
   call alignx(16,q)


   x1 = loadfp(q(1,2))
   x2 = loadfp(q(2,2))
   x3 = loadfp(q(3,2))
   x4 = loadfp(q(4,2))
   x5 = loadfp(q(5,2))

   h2 = hh(2,2)
   y1 = loadfp(q(1,1))
   y2 = loadfp(q(2,1))
   y3 = loadfp(q(3,1))
   y4 = loadfp(q(4,1))
   y5 = loadfp(q(5,1))
   y1 = fxcpmadd(y1,x1,h2)
   q1 = loadfp(q(1,3))
   y2 = fxcpmadd(y2,x2,h2)
   q2 = loadfp(q(2,3))
   y3 = fxcpmadd(y3,x3,h2)
   q3 = loadfp(q(3,3))
   y4 = fxcpmadd(y4,x4,h2)
   q4 = loadfp(q(4,3))
   y5 = fxcpmadd(y5,x5,h2)
   q5 = loadfp(q(5,3))

   h1 = hh(3-1,1)

   do i=3,nb,2

      h2 = hh(i,2)

      x1 = fxcpmadd(x1,q1,h1)
      x2 = fxcpmadd(x2,q2,h1)
      x3 = fxcpmadd(x3,q3,h1)
      x4 = fxcpmadd(x4,q4,h1)
      x5 = fxcpmadd(x5,q5,h1)

      h1 = hh(i  ,1)

      y1 = fxcpmadd(y1,q1,h2)
      q1 = loadfp(q(1,i+1))
      y2 = fxcpmadd(y2,q2,h2)
      q2 = loadfp(q(2,i+1))
      y3 = fxcpmadd(y3,q3,h2)
      q3 = loadfp(q(3,i+1))
      y4 = fxcpmadd(y4,q4,h2)
      q4 = loadfp(q(4,i+1))
      y5 = fxcpmadd(y5,q5,h2)
      q5 = loadfp(q(5,i+1))

      if(i==nb) exit

      h2 = hh(i+1,2)

      x1 = fxcpmadd(x1,q1,h1)
      x2 = fxcpmadd(x2,q2,h1)
      x3 = fxcpmadd(x3,q3,h1)
      x4 = fxcpmadd(x4,q4,h1)
      x5 = fxcpmadd(x5,q5,h1)

      h1 = hh(i+1,1)

      y1 = fxcpmadd(y1,q1,h2)
      q1 = loadfp(q(1,i+2))
      y2 = fxcpmadd(y2,q2,h2)
      q2 = loadfp(q(2,i+2))
      y3 = fxcpmadd(y3,q3,h2)
      q3 = loadfp(q(3,i+2))
      y4 = fxcpmadd(y4,q4,h2)
      q4 = loadfp(q(4,i+2))
      y5 = fxcpmadd(y5,q5,h2)
      q5 = loadfp(q(5,i+2))

   enddo

   x1 = fxcpmadd(x1,q1,h1)
   x2 = fxcpmadd(x2,q2,h1)
   x3 = fxcpmadd(x3,q3,h1)
   x4 = fxcpmadd(x4,q4,h1)
   x5 = fxcpmadd(x5,q5,h1)

   h1 = -hh(1,1) ! for below
   h2 = -hh(1,2)
   x1 = fxpmul(x1,h1)
   x2 = fxpmul(x2,h1)
   x3 = fxpmul(x3,h1)
   x4 = fxpmul(x4,h1)
   x5 = fxpmul(x5,h1)
   h1 = -hh(1,2)*s
   y1 = fxpmul(y1,h2)
   y2 = fxpmul(y2,h2)
   y3 = fxpmul(y3,h2)
   y4 = fxpmul(y4,h2)
   y5 = fxpmul(y5,h2)
   y1 = fxcpmadd(y1,x1,h1)
   q1 = loadfp(q(1,1))
   y2 = fxcpmadd(y2,x2,h1)
   q2 = loadfp(q(2,1))
   y3 = fxcpmadd(y3,x3,h1)
   q3 = loadfp(q(3,1))
   y4 = fxcpmadd(y4,x4,h1)
   q4 = loadfp(q(4,1))
   y5 = fxcpmadd(y5,x5,h1)
   q5 = loadfp(q(5,1))

   q1 = fpadd(q1,y1)
   p1 = loadfp(q(1,2))
   q2 = fpadd(q2,y2)
   p2 = loadfp(q(2,2))
   q3 = fpadd(q3,y3)
   p3 = loadfp(q(3,2))
   q4 = fpadd(q4,y4)
   p4 = loadfp(q(4,2))
   q5 = fpadd(q5,y5)
   p5 = loadfp(q(5,2))

   h2 = hh(2,2)

   call storefp(q(1,1),q1)
   p1 = fpadd(p1,x1)
   call storefp(q(2,1),q2)
   p2 = fpadd(p2,x2)
   call storefp(q(3,1),q3)
   p3 = fpadd(p3,x3)
   call storefp(q(4,1),q4)
   p4 = fpadd(p4,x4)
   call storefp(q(5,1),q5)
   p5 = fpadd(p5,x5)

   p1 = fxcpmadd(p1,y1,h2)
   q1 = loadfp(q(1,3))
   p2 = fxcpmadd(p2,y2,h2)
   q2 = loadfp(q(2,3))
   p3 = fxcpmadd(p3,y3,h2)
   q3 = loadfp(q(3,3))
   p4 = fxcpmadd(p4,y4,h2)
   q4 = loadfp(q(4,3))
   p5 = fxcpmadd(p5,y5,h2)
   q5 = loadfp(q(5,3))

   h1 = hh(3-1,1)

   do i=3,nb,2

      h2 = hh(i,2)

      call storefp(q(1,i-1),p1)
      q1 = fxcpmadd(q1,x1,h1)
      call storefp(q(2,i-1),p2)
      q2 = fxcpmadd(q2,x2,h1)
      call storefp(q(3,i-1),p3)
      q3 = fxcpmadd(q3,x3,h1)
      call storefp(q(4,i-1),p4)
      q4 = fxcpmadd(q4,x4,h1)
      call storefp(q(5,i-1),p5)
      q5 = fxcpmadd(q5,x5,h1)

      h1 = hh(i,1)

      q1 = fxcpmadd(q1,y1,h2)
      p1 = loadfp(q(1,i+1))
      q2 = fxcpmadd(q2,y2,h2)
      p2 = loadfp(q(2,i+1))
      q3 = fxcpmadd(q3,y3,h2)
      p3 = loadfp(q(3,i+1))
      q4 = fxcpmadd(q4,y4,h2)
      p4 = loadfp(q(4,i+1))
      q5 = fxcpmadd(q5,y5,h2)
      p5 = loadfp(q(5,i+1))

      if(i==nb) exit

      h2 = hh(i+1,2)

      call storefp(q(1,i),q1)
      p1 = fxcpmadd(p1,x1,h1)
      call storefp(q(2,i),q2)
      p2 = fxcpmadd(p2,x2,h1)
      call storefp(q(3,i),q3)
      p3 = fxcpmadd(p3,x3,h1)
      call storefp(q(4,i),q4)
      p4 = fxcpmadd(p4,x4,h1)
      call storefp(q(5,i),q5)
      p5 = fxcpmadd(p5,x5,h1)

      h1 = hh(i+1,1)

      p1 = fxcpmadd(p1,y1,h2)
      q1 = loadfp(q(1,i+2))
      p2 = fxcpmadd(p2,y2,h2)
      q2 = loadfp(q(2,i+2))
      p3 = fxcpmadd(p3,y3,h2)
      q3 = loadfp(q(3,i+2))
      p4 = fxcpmadd(p4,y4,h2)
      q4 = loadfp(q(4,i+2))
      p5 = fxcpmadd(p5,y5,h2)
      q5 = loadfp(q(5,i+2))

   enddo


   if(i==nb) then
      call storefp(q(1,nb),q1)
      p1 = fxcpmadd(p1,x1,h1)
      call storefp(q(2,nb),q2)
      p2 = fxcpmadd(p2,x2,h1)
      call storefp(q(3,nb),q3)
      p3 = fxcpmadd(p3,x3,h1)
      call storefp(q(4,nb),q4)
      p4 = fxcpmadd(p4,x4,h1)
      call storefp(q(5,nb),q5)
      p5 = fxcpmadd(p5,x5,h1)

      call storefp(q(1,nb+1),p1)
      call storefp(q(2,nb+1),p2)
      call storefp(q(3,nb+1),p3)
      call storefp(q(4,nb+1),p4)
      call storefp(q(5,nb+1),p5)
   else
      call storefp(q(1,nb),p1)
      q1 = fxcpmadd(q1,x1,h1)
      call storefp(q(2,nb),p2)
      q2 = fxcpmadd(q2,x2,h1)
      call storefp(q(3,nb),p3)
      q3 = fxcpmadd(q3,x3,h1)
      call storefp(q(4,nb),p4)
      q4 = fxcpmadd(q4,x4,h1)
      call storefp(q(5,nb),p5)
      q5 = fxcpmadd(q5,x5,h1)

      call storefp(q(1,nb+1),q1)
      call storefp(q(2,nb+1),q2)
      call storefp(q(3,nb+1),q3)
      call storefp(q(4,nb+1),q4)
      call storefp(q(5,nb+1),q5)
   endif


!contains
!
!   subroutine storefp(a,b)
!      complex*16 a, b
!
!      a = b
!   end subroutine
!   subroutine alignx(n, x)
!      integer n
!      complex*16 x(ldq/2,*)
!   end subroutine

end subroutine

! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_kernel_8_bg(q, hh, nb, ldq, ldh, s)


   implicit none

#if defined __parallel
   include 'mpif.h'
#endif   

   integer, intent(in) :: nb, ldq, ldh
   complex*16, intent(inout) :: q(ldq/2,*)
   real*8, intent(in) :: hh(ldh,*), s

   complex*16 x1, x2, x3, x4, y1, y2, y3, y4, q1, q2, q3, q4, p1, p2, p3, p4
   real*8 h1, h2
   integer i

!   complex*16 loadfp, fxcpmadd, fxpmul, fpadd, a, b
!   real*8 x
!   loadfp(a) = a
!   fxcpmadd(a,b,x) = a + b*x
!   fxpmul(a,x) = a*x
!   fpadd(a,b) = a+b

   call alignx(16,q)


   x1 = loadfp(q(1,2))
   x2 = loadfp(q(2,2))
   x3 = loadfp(q(3,2))
   x4 = loadfp(q(4,2))

   h2 = hh(2,2)
   y1 = loadfp(q(1,1))
   y2 = loadfp(q(2,1))
   y3 = loadfp(q(3,1))
   y4 = loadfp(q(4,1))
   y1 = fxcpmadd(y1,x1,h2)
   q1 = loadfp(q(1,3))
   y2 = fxcpmadd(y2,x2,h2)
   q2 = loadfp(q(2,3))
   y3 = fxcpmadd(y3,x3,h2)
   q3 = loadfp(q(3,3))
   y4 = fxcpmadd(y4,x4,h2)
   q4 = loadfp(q(4,3))

   h1 = hh(3-1,1)

   do i=3,nb,2

      h2 = hh(i,2)

      x1 = fxcpmadd(x1,q1,h1)
      x2 = fxcpmadd(x2,q2,h1)
      x3 = fxcpmadd(x3,q3,h1)
      x4 = fxcpmadd(x4,q4,h1)

      h1 = hh(i  ,1)

      y1 = fxcpmadd(y1,q1,h2)
      q1 = loadfp(q(1,i+1))
      y2 = fxcpmadd(y2,q2,h2)
      q2 = loadfp(q(2,i+1))
      y3 = fxcpmadd(y3,q3,h2)
      q3 = loadfp(q(3,i+1))
      y4 = fxcpmadd(y4,q4,h2)
      q4 = loadfp(q(4,i+1))

      if(i==nb) exit

      h2 = hh(i+1,2)

      x1 = fxcpmadd(x1,q1,h1)
      x2 = fxcpmadd(x2,q2,h1)
      x3 = fxcpmadd(x3,q3,h1)
      x4 = fxcpmadd(x4,q4,h1)

      h1 = hh(i+1,1)

      y1 = fxcpmadd(y1,q1,h2)
      q1 = loadfp(q(1,i+2))
      y2 = fxcpmadd(y2,q2,h2)
      q2 = loadfp(q(2,i+2))
      y3 = fxcpmadd(y3,q3,h2)
      q3 = loadfp(q(3,i+2))
      y4 = fxcpmadd(y4,q4,h2)
      q4 = loadfp(q(4,i+2))

   enddo

   x1 = fxcpmadd(x1,q1,h1)
   x2 = fxcpmadd(x2,q2,h1)
   x3 = fxcpmadd(x3,q3,h1)
   x4 = fxcpmadd(x4,q4,h1)

   h1 = -hh(1,1) ! for below
   h2 = -hh(1,2)
   x1 = fxpmul(x1,h1)
   x2 = fxpmul(x2,h1)
   x3 = fxpmul(x3,h1)
   x4 = fxpmul(x4,h1)
   h1 = -hh(1,2)*s
   y1 = fxpmul(y1,h2)
   y2 = fxpmul(y2,h2)
   y3 = fxpmul(y3,h2)
   y4 = fxpmul(y4,h2)
   y1 = fxcpmadd(y1,x1,h1)
   q1 = loadfp(q(1,1))
   y2 = fxcpmadd(y2,x2,h1)
   q2 = loadfp(q(2,1))
   y3 = fxcpmadd(y3,x3,h1)
   q3 = loadfp(q(3,1))
   y4 = fxcpmadd(y4,x4,h1)
   q4 = loadfp(q(4,1))

   q1 = fpadd(q1,y1)
   p1 = loadfp(q(1,2))
   q2 = fpadd(q2,y2)
   p2 = loadfp(q(2,2))
   q3 = fpadd(q3,y3)
   p3 = loadfp(q(3,2))
   q4 = fpadd(q4,y4)
   p4 = loadfp(q(4,2))

   h2 = hh(2,2)

   call storefp(q(1,1),q1)
   p1 = fpadd(p1,x1)
   call storefp(q(2,1),q2)
   p2 = fpadd(p2,x2)
   call storefp(q(3,1),q3)
   p3 = fpadd(p3,x3)
   call storefp(q(4,1),q4)
   p4 = fpadd(p4,x4)

   p1 = fxcpmadd(p1,y1,h2)
   q1 = loadfp(q(1,3))
   p2 = fxcpmadd(p2,y2,h2)
   q2 = loadfp(q(2,3))
   p3 = fxcpmadd(p3,y3,h2)
   q3 = loadfp(q(3,3))
   p4 = fxcpmadd(p4,y4,h2)
   q4 = loadfp(q(4,3))

   h1 = hh(3-1,1)

   do i=3,nb,2

      h2 = hh(i,2)

      call storefp(q(1,i-1),p1)
      q1 = fxcpmadd(q1,x1,h1)
      call storefp(q(2,i-1),p2)
      q2 = fxcpmadd(q2,x2,h1)
      call storefp(q(3,i-1),p3)
      q3 = fxcpmadd(q3,x3,h1)
      call storefp(q(4,i-1),p4)
      q4 = fxcpmadd(q4,x4,h1)

      h1 = hh(i,1)

      q1 = fxcpmadd(q1,y1,h2)
      p1 = loadfp(q(1,i+1))
      q2 = fxcpmadd(q2,y2,h2)
      p2 = loadfp(q(2,i+1))
      q3 = fxcpmadd(q3,y3,h2)
      p3 = loadfp(q(3,i+1))
      q4 = fxcpmadd(q4,y4,h2)
      p4 = loadfp(q(4,i+1))

      if(i==nb) exit

      h2 = hh(i+1,2)

      call storefp(q(1,i),q1)
      p1 = fxcpmadd(p1,x1,h1)
      call storefp(q(2,i),q2)
      p2 = fxcpmadd(p2,x2,h1)
      call storefp(q(3,i),q3)
      p3 = fxcpmadd(p3,x3,h1)
      call storefp(q(4,i),q4)
      p4 = fxcpmadd(p4,x4,h1)

      h1 = hh(i+1,1)

      p1 = fxcpmadd(p1,y1,h2)
      q1 = loadfp(q(1,i+2))
      p2 = fxcpmadd(p2,y2,h2)
      q2 = loadfp(q(2,i+2))
      p3 = fxcpmadd(p3,y3,h2)
      q3 = loadfp(q(3,i+2))
      p4 = fxcpmadd(p4,y4,h2)
      q4 = loadfp(q(4,i+2))

   enddo


   if(i==nb) then
      call storefp(q(1,nb),q1)
      p1 = fxcpmadd(p1,x1,h1)
      call storefp(q(2,nb),q2)
      p2 = fxcpmadd(p2,x2,h1)
      call storefp(q(3,nb),q3)
      p3 = fxcpmadd(p3,x3,h1)
      call storefp(q(4,nb),q4)
      p4 = fxcpmadd(p4,x4,h1)

      call storefp(q(1,nb+1),p1)
      call storefp(q(2,nb+1),p2)
      call storefp(q(3,nb+1),p3)
      call storefp(q(4,nb+1),p4)
   else
      call storefp(q(1,nb),p1)
      q1 = fxcpmadd(q1,x1,h1)
      call storefp(q(2,nb),p2)
      q2 = fxcpmadd(q2,x2,h1)
      call storefp(q(3,nb),p3)
      q3 = fxcpmadd(q3,x3,h1)
      call storefp(q(4,nb),p4)
      q4 = fxcpmadd(q4,x4,h1)

      call storefp(q(1,nb+1),q1)
      call storefp(q(2,nb+1),q2)
      call storefp(q(3,nb+1),q3)
      call storefp(q(4,nb+1),q4)
   endif


!contains
!
!   subroutine storefp(a,b)
!      complex*16 a, b
!
!      a = b
!   end subroutine
!   subroutine alignx(n, x)
!      integer n
!      complex*16 x(ldq/2,*)
!   end subroutine

end subroutine

! --------------------------------------------------------------------------------------------------

subroutine hh_trafo_kernel_4_bg(q, hh, nb, ldq, ldh, s)


   implicit none

#if defined __parallel
   include 'mpif.h'
#endif   

   integer, intent(in) :: nb, ldq, ldh
   complex*16, intent(inout) :: q(ldq/2,*)
   real*8, intent(in) :: hh(ldh,*), s

   complex*16 x1, x2, y1, y2, q1, q2, p1, p2
   real*8 h1, h2
   integer i

!   complex*16 loadfp, fxcpmadd, fxpmul, fpadd, a, b
!   real*8 x
!   loadfp(a) = a
!   fxcpmadd(a,b,x) = a + b*x
!   fxpmul(a,x) = a*x
!   fpadd(a,b) = a+b

   call alignx(16,q)


   x1 = loadfp(q(1,2))
   x2 = loadfp(q(2,2))

   h2 = hh(2,2)
   y1 = loadfp(q(1,1))
   y2 = loadfp(q(2,1))
   y1 = fxcpmadd(y1,x1,h2)
   q1 = loadfp(q(1,3))
   y2 = fxcpmadd(y2,x2,h2)
   q2 = loadfp(q(2,3))

   h1 = hh(3-1,1)

   do i=3,nb,2

      h2 = hh(i,2)

      x1 = fxcpmadd(x1,q1,h1)
      x2 = fxcpmadd(x2,q2,h1)

      h1 = hh(i  ,1)

      y1 = fxcpmadd(y1,q1,h2)
      q1 = loadfp(q(1,i+1))
      y2 = fxcpmadd(y2,q2,h2)
      q2 = loadfp(q(2,i+1))

      if(i==nb) exit

      h2 = hh(i+1,2)

      x1 = fxcpmadd(x1,q1,h1)
      x2 = fxcpmadd(x2,q2,h1)

      h1 = hh(i+1,1)

      y1 = fxcpmadd(y1,q1,h2)
      q1 = loadfp(q(1,i+2))
      y2 = fxcpmadd(y2,q2,h2)
      q2 = loadfp(q(2,i+2))

   enddo

   x1 = fxcpmadd(x1,q1,h1)
   x2 = fxcpmadd(x2,q2,h1)

   h1 = -hh(1,1) ! for below
   h2 = -hh(1,2)
   x1 = fxpmul(x1,h1)
   x2 = fxpmul(x2,h1)
   h1 = -hh(1,2)*s
   y1 = fxpmul(y1,h2)
   y2 = fxpmul(y2,h2)
   y1 = fxcpmadd(y1,x1,h1)
   q1 = loadfp(q(1,1))
   y2 = fxcpmadd(y2,x2,h1)
   q2 = loadfp(q(2,1))

   q1 = fpadd(q1,y1)
   p1 = loadfp(q(1,2))
   q2 = fpadd(q2,y2)
   p2 = loadfp(q(2,2))

   h2 = hh(2,2)

   call storefp(q(1,1),q1)
   p1 = fpadd(p1,x1)
   call storefp(q(2,1),q2)
   p2 = fpadd(p2,x2)

   p1 = fxcpmadd(p1,y1,h2)
   q1 = loadfp(q(1,3))
   p2 = fxcpmadd(p2,y2,h2)
   q2 = loadfp(q(2,3))

   h1 = hh(3-1,1)

   do i=3,nb,2

      h2 = hh(i,2)

      call storefp(q(1,i-1),p1)
      q1 = fxcpmadd(q1,x1,h1)
      call storefp(q(2,i-1),p2)
      q2 = fxcpmadd(q2,x2,h1)

      h1 = hh(i,1)

      q1 = fxcpmadd(q1,y1,h2)
      p1 = loadfp(q(1,i+1))
      q2 = fxcpmadd(q2,y2,h2)
      p2 = loadfp(q(2,i+1))

      if(i==nb) exit

      h2 = hh(i+1,2)

      call storefp(q(1,i),q1)
      p1 = fxcpmadd(p1,x1,h1)
      call storefp(q(2,i),q2)
      p2 = fxcpmadd(p2,x2,h1)

      h1 = hh(i+1,1)

      p1 = fxcpmadd(p1,y1,h2)
      q1 = loadfp(q(1,i+2))
      p2 = fxcpmadd(p2,y2,h2)
      q2 = loadfp(q(2,i+2))

   enddo


   if(i==nb) then
      call storefp(q(1,nb),q1)
      p1 = fxcpmadd(p1,x1,h1)
      call storefp(q(2,nb),q2)
      p2 = fxcpmadd(p2,x2,h1)

      call storefp(q(1,nb+1),p1)
      call storefp(q(2,nb+1),p2)
   else
      call storefp(q(1,nb),p1)
      q1 = fxcpmadd(q1,x1,h1)
      call storefp(q(2,nb),p2)
      q2 = fxcpmadd(q2,x2,h1)

      call storefp(q(1,nb+1),q1)
      call storefp(q(2,nb+1),q2)
   endif


!contains
!
!   subroutine storefp(a,b)
!      complex*16 a, b
!
!      a = b
!   end subroutine
!   subroutine alignx(n, x)
!      integer n
!      complex*16 x(ldq/2,*)
!   end subroutine

end subroutine

! --------------------------------------------------------------------------------------------------
