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

!> @file tblite/integral/diat_trafo.f90
!> Evaluation of the diatomic scaled overlap
module tblite_integral_diat_trafo
   use mctc_env, only : wp
   use tblite_blas, only: gemm

   implicit none
   private

   public :: diat_trafo, diat_trafo_grad

   !> Dimension of trafo matrix for highest angular momentum.
   integer, parameter :: sdim(0:6) = [1, 4, 9, 16, 25, 36, 49]

contains

   !> Transformation to the diatomic frame and back: 
   pure subroutine diat_trafo(block_overlap, vec, ksig, kpi, kdel, maxlj, maxli)
      !> Diatomic block of CGTOs to be transformed (+ scaled)
      real(wp),intent(inout)    :: block_overlap(:,:)
      !> Transformation vector for the diatomic frame (i.e. vector between the two centers)
      real(wp),intent(in)       :: vec(3)
      !> Scaling parameters for different bonding contributions
      real(wp),intent(in)       :: ksig, kpi, kdel
      !> Highest angular momentum of atom j (first index)
      integer,intent(in)        :: maxlj
      !> Highest angular momentum of atom i (second index)
      integer,intent(in)        :: maxli

      integer  :: dimj, dimi, maxl
      real(wp), allocatable :: trafomat(:,:), tmp(:,:), transformed_s(:,:)

      ! Select the dimensions of the transformation matrix
      dimj = sdim(maxlj)
      dimi = sdim(maxli)
      maxl = max(maxlj, maxli)

      allocate(trafomat(sdim(maxl), sdim(maxl)), transformed_s(dimj,dimi), &
      & tmp(dimj,dimi), source=0.0_wp)

      ! 1. Setup the transformation matrix
      call harmtr(maxl, vec, trafomat)

      ! 2. Transform the overlap submatrix to the diatomic frame: S' = O^T * S * O
      if (maxl > 0) then
         call gemm(amat=trafomat(1:dimj,1:dimj), bmat=block_overlap(1:dimj,1:dimi), &
            & cmat=tmp, transa='T', transb='N')
         call gemm(amat=tmp, bmat=trafomat(1:dimi,1:dimi), &
            & cmat=transformed_s, transa='N', transb='N')
      else
         transformed_s(1,1) = block_overlap(1,1)
      endif

      ! 3. Scale elements in the diatomic frame
      call scale_diatomic_frame(transformed_s, ksig, kpi, kdel, maxlj, maxli) 

      ! 4. Transform the overlap submatrix back to original frame: Ssc = O * Ssc' * O^T
      block_overlap = 0.0_wp
      if (maxl > 0) then
         call gemm(amat=trafomat(1:dimj,1:dimj), bmat=transformed_s, &
            & cmat=tmp, transa='N', transb='N')
         call gemm(amat=tmp, bmat=trafomat(1:dimi,1:dimi), &
            & cmat=block_overlap(1:dimj,1:dimi), transa='N', transb='T')
      else
         block_overlap(1,1) = transformed_s(1,1)
      endif

   end subroutine diat_trafo

   !> Gradient of the diatomic frame scaled overlap transformation: 
   pure subroutine diat_trafo_grad(block_overlap, block_doverlap, vec, ksig, kpi, kdel, maxlj, maxli)
      !> Diatomic block of CGTO overlap to be transformed (+ scaled)
      real(wp),intent(inout)    :: block_overlap(:,:)
      !> Derivative of diatomic block of CGTO overlap to be transformed (+ scaled)
      real(wp),intent(inout)    :: block_doverlap(:,:,:)
      !> Transformation vector for the diatomic frame (i.e. vector between the two centers)
      real(wp),intent(in)       :: vec(3)
      !> Scaling parameters for different bonding contributions
      real(wp),intent(in)       :: ksig, kpi, kdel
      !> Highest angular momentum of atom j (first index)
      integer,intent(in)        :: maxlj
      !> Highest angular momentum of atom i (second index)
      integer,intent(in)        :: maxli
      
      integer :: ic, dimj, dimi, maxl
      real(wp), allocatable :: trafomat(:,:,:), dtrafomat(:,:,:), tmp(:,:), &
      & interm_oso(:,:), interm_doso(:,:,:), interm_odso(:,:,:), interm_osdo(:,:,:)

      ! Select the dimensions of the transformation matrix
      dimj = sdim(maxlj)
      dimi = sdim(maxli)
      maxl = max(maxlj, maxli)

      allocate(trafomat(3,sdim(maxl), sdim(maxl)), dtrafomat(3,sdim(maxl), sdim(maxl)), &
      & interm_oso(dimj,dimi), interm_doso(3,dimj,dimi), interm_osdo(3,dimj, dimi), &
      & interm_odso(3,dimj,dimi), tmp(dimj,dimi), source=0.0_wp)

      ! 1. Setup the transformation matrix and its derivative for all directions.
      ! For the case vec || z-axis, the x- and y-derivatives are ill-defined 
      ! and have to be evaluated assuming either x or y orientation. 
      call d_harmtr(maxl, vec, trafomat, dtrafomat)

      ! 2. Transform the overlap submatrix to the diatomic frame: S' = O^T * S * O
      if (maxl > 0) then
         ! interm_oso = O^T * S * O
         call gemm(amat=trafomat(3,1:dimj,1:dimj), bmat=block_overlap(1:dimj,1:dimi), &
            & cmat=tmp,transa='T', transb='N')
         call gemm(amat=tmp, bmat=trafomat(3,1:dimi,1:dimi), &
            & cmat=interm_oso, transa='N', transb='N')
         do ic = 1, 3
            ! interm_doso = dO^T * S * O
            call gemm(amat=dtrafomat(ic,1:dimj,1:dimj), bmat=block_overlap(1:dimj,1:dimi), &
               & cmat=tmp, transa='T', transb='N')
            call gemm(amat=tmp, bmat=trafomat(ic,1:dimi,1:dimi), &
               & cmat=interm_doso(ic,:,:), transa='N', transb='N')
            ! interm_osdo = O^T * S * dO
            call gemm(amat=trafomat(ic,1:dimj,1:dimj), bmat=block_overlap(1:dimj,1:dimi), &
               & cmat=tmp, transa='T', transb='N')
            call gemm(amat=tmp, bmat=dtrafomat(ic,1:dimi,1:dimi), &
               & cmat=interm_osdo(ic,:,:), transa='N', transb='N')
            ! interm_odso = O^T * dS * O
            call gemm(amat=trafomat(ic,1:dimj,1:dimj), bmat=block_doverlap(ic,1:dimj,1:dimi), &
               & cmat=tmp, transa='T', transb='N')
            call gemm(amat=tmp, bmat=trafomat(ic,1:dimi,1:dimi), &
               & cmat=interm_odso(ic,:,:), transa='N', transb='N')
         end do
      else
         interm_oso(1,1) = block_overlap(1,1)
         interm_doso(:,1,1) = 0.0_wp
         interm_odso(:,1,1) = block_doverlap(:,1,1)
         interm_osdo(:,1,1) = 0.0_wp
      endif

      ! 3. Scale overlap and each dimension of the derivative in the diatomic frame
      call scale_diatomic_frame(interm_oso, ksig, kpi, kdel, maxlj, maxli) 
      do ic = 1, 3
         call scale_diatomic_frame(interm_doso(ic,:,:), ksig, kpi, kdel, maxlj, maxli)
         call scale_diatomic_frame(interm_osdo(ic,:,:), ksig, kpi, kdel, maxlj, maxli) 
         call scale_diatomic_frame(interm_odso(ic,:,:), ksig, kpi, kdel, maxlj, maxli) 
      end do

      ! 4. Transform diatomic frame quantities (S', (dOSO)', (OSdO)' and (OdSO)') back to original frame
      if (maxl > 0) then
         ! block_overlap = O * S' * O^T
         call gemm(amat=trafomat(3,1:dimj,1:dimj), bmat=interm_oso, &
            & cmat=tmp, transa='N', transb='N')
         call gemm(amat=tmp, bmat=trafomat(3,1:dimi,1:dimi), &
            & cmat=block_overlap(1:dimj,1:dimi), transa='N', transb='T')

         do ic = 1, 3
            ! block_doverlap = dO * S' * O^T + O * S' * dO^T
            call gemm(amat=dtrafomat(ic,1:dimj,1:dimj), bmat=interm_oso, &
               & cmat=tmp, transa='N', transb='N')
            call gemm(alpha=1.0_wp, amat=tmp, bmat=trafomat(ic,1:dimi,1:dimi), &
               & cmat=block_doverlap(ic,1:dimj,1:dimi), transa='N', transb='T')

            call gemm(amat=trafomat(ic,1:dimj,1:dimj), bmat=interm_oso, &
               & cmat=tmp, transa='N', transb='N')
            call gemm(alpha=1.0_wp, amat=tmp, bmat=dtrafomat(ic,1:dimi,1:dimi), beta=1.0_wp, &
               & cmat=block_doverlap(ic,1:dimj,1:dimi), transa='N', transb='T')
            
            ! block_doverlap += O * (dOSO)' * O^T + O * (OSdO)' * O^T 
            call gemm(amat=trafomat(ic,1:dimj,1:dimj), bmat=interm_doso(ic,:,:), &
               & cmat=tmp, transa='N', transb='N')
            call gemm(alpha=1.0_wp, amat=tmp, bmat=trafomat(ic,1:dimi,1:dimi), beta=1.0_wp, &
               & cmat=block_doverlap(ic,1:dimj,1:dimi), transa='N', transb='T')

            call gemm(amat=trafomat(ic,1:dimj,1:dimj), bmat=interm_odso(ic,:,:), &
               & cmat=tmp, transa='N', transb='N')
            call gemm(amat=tmp, bmat=trafomat(ic,1:dimi,1:dimi), beta=1.0_wp, &
               & cmat=block_doverlap(ic,1:dimj,1:dimi), transa='N', transb='T')
            
            ! block_doverlap += O * (OdSO)' * O^T 
            call gemm(amat=trafomat(ic,1:dimj,1:dimj), bmat=interm_osdo(ic,:,:), &
               & cmat=tmp, transa='N', transb='N')
            call gemm(alpha=1.0_wp, amat=tmp, bmat=trafomat(ic,1:dimi,1:dimi), beta=1.0_wp, &
               & cmat=block_doverlap(ic,1:dimj,1:dimi), transa='N', transb='T')
         end do
      else
         block_overlap(1,1) = interm_oso(1,1)
         block_doverlap(:,1,1) = interm_odso(:,1,1)
      endif

   end subroutine diat_trafo_grad


   pure subroutine harmtr(maxl,vec,trafomat)
      !> Maximum angular momentum
      integer, intent(in)  :: maxl
      !> Normalized vector from atom k to atom l
      real(wp), intent(in) :: vec(3)
      !> Transformation matrix
      real(wp), intent(out) :: trafomat(sdim(maxl),sdim(maxl))
      real(wp) :: cos2p, cos2t, cosp, cost, sin2p, sin2t, sinp, sint, sqrt3, len
      real(wp) :: norm_vec(3)

      trafomat = 0.0_wp

      ! -----------------------------
      !  s functions (trafomat(1x1))
      ! -----------------------------
      
      trafomat(1,1) = 1.0_wp

      if ( maxl == 0 ) return

      ! Normalize the vector
      len = sqrt(sum(vec**2))
      norm_vec = vec / len

      ! Prepare spherical coordinats
      cost = norm_vec(3)
      if ( abs(cost) .eq. 1.0_wp ) then
         ! Here, phi is arbitrary as the vector is parallel to the z-axis.
         ! We choose the x-axis as the arbitrary direction.
         sint = 0.0_wp
         cosp = 1.0_wp
         sinp = 0.0_wp
      else if ( abs(cost) .eq. 0.0_wp ) then
         sint = 1.0_wp
         cosp = norm_vec(1)
         sinp = norm_vec(2)
      else
         sint = SQRT(norm_vec(1)**2 + norm_vec(2)**2)
         cosp = norm_vec(1)/SINT
         sinp = norm_vec(2)/SINT
      endif

      ! -----------------------------
      !  p functions (trafomat(4x4))
      ! -----------------------------

      ! Adapted to tblite ordering from MSINDO
      ! 1st index:
      ! (2,:)_MSINDO -> (px,:) -> (4,:)_tblite
      ! (3,:)_MSINDO -> (py,:) -> (2,:)_tblite
      ! (4,:)_MSINDO -> (pz,:) -> (3,:)_tblite
      trafomat(4,3) = SINT*COSP
      trafomat(2,3) = SINT*SINP
      trafomat(3,3) = COST
      trafomat(4,4) = COST*COSP
      trafomat(2,4) = COST*SINP
      trafomat(3,4) = -SINT
      trafomat(4,2) = -SINP
      trafomat(2,2) = COSP
      trafomat(3,2) = 0.0_wp

      if ( maxl <= 1 ) return

      ! -----------------------------
      !  d functions (trafomat(9x9)) 
      ! -----------------------------

      COS2T = COST**2 - SINT**2
      SIN2T = 2.0_wp * SINT*COST
      COS2P = COSP**2 - SINP**2
      SIN2P = 2.0_wp * SINP*COSP
      SQRT3 = SQRT(3.0_wp)

      ! Changed from MSINDO ordering (0,-1,1,-2,2) to tblite ordering (-2,-1,0,1,2) of d-functions 
      ! (5,:)_MSINDO -> (dz2,:) -> (7,:)_tblite
      ! (6,:)_MSINDO -> (dxz,:) -> (8,:)_tblite
      ! (7,:)_MSINDO -> (dyz,:) -> (6,:)_tblite
      ! (8,:)_MSINDO -> (dx2-y2,:) -> (5,:)_tblite
      ! (9,:)_MSINDO -> (dxy,:) -> (9,:)_tblite
      trafomat(7,7) = (3.0_wp * COST**2 - 1.0_wp) * 0.5_wp
      trafomat(8,7) = SQRT3*SIN2T*COSP*0.5_wp
      trafomat(6,7) = SQRT3*SIN2T*SINP*0.5_wp
      trafomat(9,7) = SQRT3*SINT**2*COS2P*0.5_wp
      trafomat(5,7) = SQRT3*SINT**2*SIN2P*0.5_wp
      trafomat(7,8) = -SQRT3*SIN2T*0.5_wp
      trafomat(8,8) = COS2T*COSP
      trafomat(6,8) = COS2T*SINP
      trafomat(9,8) = SIN2T*COS2P*0.5_wp
      trafomat(5,8) = SIN2T*SIN2P*0.5_wp
      trafomat(7,6) = 0.0_wp
      trafomat(8,6) = -COST*SINP
      trafomat(6,6) = COST*COSP
      trafomat(9,6) = -SINT*SIN2P
      trafomat(5,6) = SINT*COS2P
      trafomat(7,9) = SQRT3*SINT**2 * 0.5_wp
      trafomat(8,9) = -SIN2T*COSP*0.5_wp
      trafomat(6,9) = -SIN2T*SINP*0.5_wp
      trafomat(9,9) = (1.0_wp + COST**2) * COS2P * 0.5_wp
      trafomat(5,9) = (1.0_wp + COST**2) * SIN2P * 0.5_wp
      trafomat(7,5) = 0.0_wp
      trafomat(8,5) = SINT*SINP
      trafomat(6,5) = -SINT*COSP
      trafomat(9,5) = -COST*SIN2P
      trafomat(5,5) = COST*COS2P

      if ( maxl <= 2 ) return

      ! -------------------------------
      !  f functions (trafomat(16x16)) 
      ! -------------------------------

      ! f-functions are not transformed to the diatomic frame!
      trafomat(10,10) = 1.0_wp
      trafomat(11,11) = 1.0_wp
      trafomat(12,12) = 1.0_wp
      trafomat(13,13) = 1.0_wp
      trafomat(14,14) = 1.0_wp
      trafomat(15,15) = 1.0_wp
      trafomat(16,16) = 1.0_wp

      if ( maxl <= 3 ) return

      ! -------------------------------
      !  g functions (trafomat(25x25)) 
      ! -------------------------------

      ! g-functions are not transformed to the diatomic frame!
      trafomat(17,17) = 1.0_wp
      trafomat(18,18) = 1.0_wp
      trafomat(19,19) = 1.0_wp
      trafomat(20,20) = 1.0_wp
      trafomat(21,21) = 1.0_wp
      trafomat(22,22) = 1.0_wp
      trafomat(23,23) = 1.0_wp
      trafomat(24,24) = 1.0_wp
      trafomat(25,25) = 1.0_wp

      if ( maxl <= 4 ) return

   end subroutine harmtr

   pure subroutine d_harmtr(maxl,vec,trafomat, dtrafomat)
      !> Maximum angular momentum
      integer, intent(in)  :: maxl
      !> Normalized vector from atom k to atom l
      real(wp), intent(in) :: vec(3)
      !> Transformation matrix
      real(wp), intent(out) :: trafomat(3,sdim(maxl),sdim(maxl))
      !> Derivative of transformation matrix
      real(wp), intent(out) :: dtrafomat(3,sdim(maxl),sdim(maxl))
      
      real(wp), parameter              :: eps = 1.0e-08_wp

      !> Derivative of transformation matrix w.r.t. theta (x- or z-direction)
      real(wp) :: trafomat_dt(sdim(maxl),sdim(maxl))
      !> Derivative of transformation matrix w.r.t. theta (y-direction)
      real(wp) :: trafomat_dty(sdim(maxl),sdim(maxl))
      !> Derivative of transformation matrix w.r.t. phi (x- or z-direction)
      real(wp) :: trafomat_dp(sdim(maxl),sdim(maxl))
      !> Derivative of transformation matrix w.r.t. phi (y-direction)
      real(wp) :: trafomat_dpy(sdim(maxl),sdim(maxl))

      ! Intermediate variables for the trigonometric functions
      ! Separte version for y for the case: vec || z-axis
      real(wp) :: cos2p, cos2t, cosp, cost, sin2p, sin2t, sinp
      real(wp) :: cos2py, cospy, sin2py, sinpy
      real(wp) :: dcos2t, dsin2t, dcos2p, dsin2p
      real(wp) :: dcos2py, dsin2py
      real(wp) :: dpdx, dpdy, dpdz, dtdx, dtdy, dtdz
      real(wp) :: norm_vec(3), sint, sqrt3, len

      trafomat = 0.0_wp
      trafomat_dt = 0.0_wp
      trafomat_dty = 0.0_wp
      trafomat_dp = 0.0_wp
      trafomat_dpy = 0.0_wp

      ! -----------------------------
      !  s functions (trafomat(1x1)) 
      ! -----------------------------

      trafomat(:,1,1) = 1.0_wp
      dtrafomat(:,1,1) = 0.0_wp

      if ( maxl == 0 ) return

      ! Normalize the vector
      len = sqrt(sum(vec**2))
      norm_vec = vec / len

      ! Prepare spherical coordinats
      cost = norm_vec(3)
      if ( abs(cost) .eq. 1.0_wp ) then
         sint = 0.0_wp
         ! Here, phi is arbitrary as the vector is parallel to the z-axis.
         ! In turn, the derivative is ill defined and has to be evaluated
         ! assuming either x or y orientation for the infinitesimal change. 
         ! This is only require for the p- and d-functions, which are transformed.
         cosp = 1.0_wp
         sinp = 0.0_wp
         if ( maxl == 1 .or. maxl == 2 ) then 
            cospy = 0.0_wp
            sinpy = -1.0_wp
         else
            cospy = cosp
            sinpy = sinp
         end if
      else if ( abs(cost) .eq. 0.0_wp ) then
         sint = 1.0_wp
         cosp = norm_vec(1)
         sinp = norm_vec(2)
         cospy = cosp
         sinpy = sinp
      else
         sint = SQRT(norm_vec(1)**2 + norm_vec(2)**2)
         cosp = norm_vec(1)/SINT
         sinp = norm_vec(2)/SINT
         cospy = cosp
         sinpy = sinp
      endif

      ! Prepare sperical coordinate derivative
      ! In the case of (exactly) vec || z-axis, the phi derivative vanishes.
      if( norm_vec(1)**2 + norm_vec(2)**2 .eq. 0.0_wp ) then
         dpdx = 0.0_wp
         dpdy = 0.0_wp
      else
         dpdx = -sinp / sqrt(vec(1)**2 + vec(2)**2) 
         dpdy = cospy / sqrt(vec(1)**2 + vec(2)**2)
      end if 
      dpdz = 0.0_wp
      dtdx = cost * cosp / len
      dtdy = cost * sinpy / len
      dtdz = -sint / len

      ! -----------------------------
      !  p functions (trafomat(4x4)) 
      ! -----------------------------

      ! Adapted to tblite ordering from MSINDO
      ! 1st index:
      ! (2,:)_MSINDO -> (px,:) -> (4,:)_tblite
      ! (3,:)_MSINDO -> (py,:) -> (2,:)_tblite
      ! (4,:)_MSINDO -> (pz,:) -> (3,:)_tblite

      ! x-direction
      trafomat(1,4,3) = SINT*COSP
      trafomat(1,2,3) = SINT*SINP
      trafomat(1,3,3) = COST
      trafomat(1,4,4) = COST*COSP
      trafomat(1,2,4) = COST*SINP
      trafomat(1,3,4) = -SINT
      trafomat(1,4,2) = -SINP
      trafomat(1,2,2) = COSP
      trafomat(1,3,2) = 0.0_wp
      ! y-direction
      trafomat(2,4,3) = SINT*COSPY
      trafomat(2,2,3) = SINT*SINPY
      trafomat(2,3,3) = COST
      trafomat(2,4,4) = COST*COSPY
      trafomat(2,2,4) = COST*SINPY
      trafomat(2,3,4) = -SINT
      trafomat(2,4,2) = -SINPY
      trafomat(2,2,2) = COSPY
      trafomat(2,3,2) = 0.0_wp
      ! z-direction equal to x-direction
      trafomat(3,:,:) = trafomat(1,:,:)

      ! derivative w.r.t. theta (x- or z-direction)
      trafomat_dt(4,3) = COST*COSP
      trafomat_dt(2,3) = COST*SINP
      trafomat_dt(3,3) = -SINT
      trafomat_dt(4,4) = -SINT*COSP
      trafomat_dt(2,4) = -SINT*SINP
      trafomat_dt(3,4) = -COST
      trafomat_dt(4,2) = 0.0_wp
      trafomat_dt(2,2) = 0.0_wp
      trafomat_dt(3,2) = 0.0_wp
      ! y-direction
      trafomat_dty(4,3) = COST*COSPY
      trafomat_dty(2,3) = COST*SINPY
      trafomat_dty(3,3) = -SINT
      trafomat_dty(4,4) = -SINT*COSPY
      trafomat_dty(2,4) = -SINT*SINPY
      trafomat_dty(3,4) = -COST
      trafomat_dty(4,2) = 0.0_wp
      trafomat_dty(2,2) = 0.0_wp
      trafomat_dty(3,2) = 0.0_wp

      ! derivative w.r.t. phi (x- or z-direction)
      trafomat_dp(4,3) = -SINT*SINP
      trafomat_dp(2,3) = SINT*COSP
      trafomat_dp(3,3) = 0.0_wp
      trafomat_dp(4,4) = -COST*SINP
      trafomat_dp(2,4) = COST*COSP
      trafomat_dp(3,4) = 0.0_wp
      trafomat_dp(4,2) = -COSP
      trafomat_dp(2,2) = -SINP
      trafomat_dp(3,2) = 0.0_wp
      ! y-direction
      trafomat_dpy(4,3) = -SINT*SINPY
      trafomat_dpy(2,3) = SINT*COSPY
      trafomat_dpy(3,3) = 0.0_wp
      trafomat_dpy(4,4) = -COST*SINPY
      trafomat_dpy(2,4) = COST*COSPY
      trafomat_dpy(3,4) = 0.0_wp
      trafomat_dpy(4,2) = -COSPY
      trafomat_dpy(2,2) = -SINPY
      trafomat_dpy(3,2) = 0.0_wp
           

      if ( maxl <= 1 ) then 
         dtrafomat(1, 1:4, 1:4) = dpdx * trafomat_dp(1:4, 1:4) + dtdx * trafomat_dt(1:4, 1:4) 
         dtrafomat(2, 1:4, 1:4) = dpdy * trafomat_dpy(1:4, 1:4) + dtdy * trafomat_dty(1:4, 1:4) 
         dtrafomat(3, 1:4, 1:4) = dtdz * trafomat_dt(1:4, 1:4)
         return
      end if

      ! -----------------------------
      !  d functions (trafomat(9x9)) 
      ! -----------------------------

      SQRT3 = SQRT(3.0_wp)
      COS2T = COST**2 - SINT**2
      SIN2T = 2.0_wp * SINT*COST
      COS2P = COSP**2 - SINP**2
      SIN2P = 2.0_wp * SINP*COSP
      COS2PY = COSPY**2 - SINPY**2
      SIN2PY = 2.0_wp * SINPY*COSPY
      
      DCOS2T = -2.0_wp * SIN2T
      DSIN2T =  2.0_wp * COS2T
      DCOS2P = -2.0_wp * SIN2P
      DSIN2P =  2.0_wp * COS2P
      DCOS2PY = -2.0_wp * SIN2PY
      DSIN2PY =  2.0_wp * COS2PY


      ! Changed from MSINDO ordering (0,-1,1,-2,2) to tblite ordering (-2,-1,0,1,2) of d-functions 
      ! (5,:)_MSINDO -> (dz2,:) -> (7,:)_tblite
      ! (6,:)_MSINDO -> (dxz,:) -> (8,:)_tblite
      ! (7,:)_MSINDO -> (dyz,:) -> (6,:)_tblite
      ! (8,:)_MSINDO -> (dx2-y2,:) -> (5,:)_tblite
      ! (9,:)_MSINDO -> (dxy,:) -> (9,:)_tblite
      ! x-direction
      trafomat(1,7,7) = (3.0_wp * COST**2 - 1.0_wp) * 0.5_wp
      trafomat(1,8,7) = SQRT3*SIN2T*COSP*0.5_wp
      trafomat(1,6,7) = SQRT3*SIN2T*SINP*0.5_wp
      trafomat(1,9,7) = SQRT3*SINT**2*COS2P*0.5_wp
      trafomat(1,5,7) = SQRT3*SINT**2*SIN2P*0.5_wp
      trafomat(1,7,8) = -SQRT3*SIN2T*0.5_wp
      trafomat(1,8,8) = COS2T*COSP
      trafomat(1,6,8) = COS2T*SINP
      trafomat(1,9,8) = SIN2T*COS2P*0.5_wp
      trafomat(1,5,8) = SIN2T*SIN2P*0.5_wp
      trafomat(1,7,6) = 0.0_wp
      trafomat(1,8,6) = -COST*SINP
      trafomat(1,6,6) = COST*COSP
      trafomat(1,9,6) = -SINT*SIN2P
      trafomat(1,5,6) = SINT*COS2P
      trafomat(1,7,9) = SQRT3*SINT**2 * 0.5_wp
      trafomat(1,8,9) = -SIN2T*COSP*0.5_wp
      trafomat(1,6,9) = -SIN2T*SINP*0.5_wp
      trafomat(1,9,9) = (1.0_wp + COST**2) * COS2P * 0.5_wp
      trafomat(1,5,9) = (1.0_wp + COST**2) * SIN2P * 0.5_wp
      trafomat(1,7,5) = 0.0_wp
      trafomat(1,8,5) = SINT*SINP
      trafomat(1,6,5) = -SINT*COSP
      trafomat(1,9,5) = -COST*SIN2P
      trafomat(1,5,5) = COST*COS2P
      ! y-direction
      trafomat(2,7,7) = (3.0_wp * COST**2 - 1.0_wp) * 0.5_wp
      trafomat(2,8,7) = SQRT3*SIN2T*COSPY*0.5_wp
      trafomat(2,6,7) = SQRT3*SIN2T*SINPY*0.5_wp
      trafomat(2,9,7) = SQRT3*SINT**2*COS2PY*0.5_wp
      trafomat(2,5,7) = SQRT3*SINT**2*SIN2PY*0.5_wp
      trafomat(2,7,8) = -SQRT3*SIN2T*0.5_wp
      trafomat(2,8,8) = COS2T*COSPY
      trafomat(2,6,8) = COS2T*SINPY
      trafomat(2,9,8) = SIN2T*COS2PY*0.5_wp
      trafomat(2,5,8) = SIN2T*SIN2PY*0.5_wp
      trafomat(2,7,6) = 0.0_wp
      trafomat(2,8,6) = -COST*SINPY
      trafomat(2,6,6) = COST*COSPY
      trafomat(2,9,6) = -SINT*SIN2PY
      trafomat(2,5,6) = SINT*COS2PY
      trafomat(2,7,9) = SQRT3*SINT**2 * 0.5_wp
      trafomat(2,8,9) = -SIN2T*COSPY*0.5_wp
      trafomat(2,6,9) = -SIN2T*SINPY*0.5_wp
      trafomat(2,9,9) = (1.0_wp + COST**2) * COS2PY * 0.5_wp
      trafomat(2,5,9) = (1.0_wp + COST**2) * SIN2PY * 0.5_wp
      trafomat(2,7,5) = 0.0_wp
      trafomat(2,8,5) = SINT*SINPY
      trafomat(2,6,5) = -SINT*COSPY
      trafomat(2,9,5) = -COST*SIN2PY
      trafomat(2,5,5) = COST*COS2PY
      ! z-direction
      trafomat(3,:,:) = trafomat(1,:,:)
      
      ! derivative w.r.t. theta (x- or z-direction)
      trafomat_dt(7,7) = -3.0_wp*SIN2T*0.5_wp
      trafomat_dt(8,7) = SQRT3*DSIN2T*COSP*0.5_wp
      trafomat_dt(6,7) = SQRT3*DSIN2T*SINP*0.5_wp
      trafomat_dt(9,7) = SQRT3*SIN2T*COS2P*0.5_wp
      trafomat_dt(5,7) = SQRT3*SIN2T*SIN2P*0.5_wp
      trafomat_dt(7,8) = -SQRT3*DSIN2T*0.5_wp
      trafomat_dt(8,8) = DCOS2T*COSP
      trafomat_dt(6,8) = DCOS2T*SINP
      trafomat_dt(9,8) = DSIN2T*COS2P*0.5_wp
      trafomat_dt(5,8) = DSIN2T*SIN2P*0.5_wp
      trafomat_dt(7,6) = 0.0_wp
      trafomat_dt(8,6) = SINT*SINP
      trafomat_dt(6,6) = -SINT*COSP
      trafomat_dt(9,6) = -COST*SIN2P
      trafomat_dt(5,6) = COST*COS2P
      trafomat_dt(7,9) = SQRT3*SIN2T*0.5_wp
      trafomat_dt(8,9) = -DSIN2T*COSP*0.5_wp
      trafomat_dt(6,9) = -DSIN2T*SINP*0.5_wp
      trafomat_dt(9,9) = -SIN2T*COS2P*0.5_wp
      trafomat_dt(5,9) = -SIN2T*SIN2P*0.5_wp
      trafomat_dt(7,5) = 0.0_wp
      trafomat_dt(8,5) = COST*SINP
      trafomat_dt(6,5) = -COST*COSP
      trafomat_dt(9,5) = SINT*SIN2P
      trafomat_dt(5,5) = -SINT*COS2P
      ! y-direction
      trafomat_dty(7,7) = -3.0_wp*SIN2T*0.5_wp
      trafomat_dty(8,7) = SQRT3*DSIN2T*COSPY*0.5_wp
      trafomat_dty(6,7) = SQRT3*DSIN2T*SINPY*0.5_wp
      trafomat_dty(9,7) = SQRT3*SIN2T*COS2PY*0.5_wp
      trafomat_dty(5,7) = SQRT3*SIN2T*SIN2PY*0.5_wp
      trafomat_dty(7,8) = -SQRT3*DSIN2T*0.5_wp
      trafomat_dty(8,8) = DCOS2T*COSPY
      trafomat_dty(6,8) = DCOS2T*SINPY
      trafomat_dty(9,8) = DSIN2T*COS2PY*0.5_wp
      trafomat_dty(5,8) = DSIN2T*SIN2PY*0.5_wp
      trafomat_dty(7,6) = 0.0_wp
      trafomat_dty(8,6) = SINT*SINPY
      trafomat_dty(6,6) = -SINT*COSPY
      trafomat_dty(9,6) = -COST*SIN2PY
      trafomat_dty(5,6) = COST*COS2PY
      trafomat_dty(7,9) = SQRT3*SIN2T*0.5_wp
      trafomat_dty(8,9) = -DSIN2T*COSPY*0.5_wp
      trafomat_dty(6,9) = -DSIN2T*SINPY*0.5_wp
      trafomat_dty(9,9) = -SIN2T*COS2PY*0.5_wp
      trafomat_dty(5,9) = -SIN2T*SIN2PY*0.5_wp
      trafomat_dty(7,5) = 0.0_wp
      trafomat_dty(8,5) = COST*SINPY
      trafomat_dty(6,5) = -COST*COSPY
      trafomat_dty(9,5) = SINT*SIN2PY
      trafomat_dty(5,5) = -SINT*COS2PY

      ! derivative w.r.t. phi (x- or z-direction)
      trafomat_dp(7,7) = 0.0_wp
      trafomat_dp(8,7) = -SQRT3*SIN2T*SINP*0.5_wp
      trafomat_dp(6,7) = SQRT3*SIN2T*COSP*0.5_wp
      trafomat_dp(9,7) = SQRT3*SINT**2*DCOS2P*0.5_wp
      trafomat_dp(5,7) = SQRT3*SINT**2*DSIN2P*0.5_wp
      trafomat_dp(7,8) = 0.0_wp
      trafomat_dp(8,8) = -COS2T*SINP
      trafomat_dp(6,8) = COS2T*COSP
      trafomat_dp(9,8) = SIN2T*DCOS2P*0.5_wp
      trafomat_dp(5,8) = SIN2T*DSIN2P*0.5_wp
      trafomat_dp(7,6) = 0.0_wp
      trafomat_dp(8,6) = -COST*COSP
      trafomat_dp(6,6) = -COST*SINP
      trafomat_dp(9,6) = -SINT*DSIN2P
      trafomat_dp(5,6) = SINT*DCOS2P
      trafomat_dp(7,9) = 0.0_wp
      trafomat_dp(8,9) = SIN2T*SINP*0.5_wp
      trafomat_dp(6,9) = -SIN2T*COSP*0.5_wp
      trafomat_dp(9,9) = (1.0_wp + COST**2)*DCOS2P*0.5_wp
      trafomat_dp(5,9) = (1.0_wp + COST**2)*DSIN2P*0.5_wp
      trafomat_dp(7,5) = 0.0_wp
      trafomat_dp(8,5) = SINT*COSP
      trafomat_dp(6,5) = SINT*SINP
      trafomat_dp(9,5) = -COST*DSIN2P
      trafomat_dp(5,5) = COST*DCOS2P
      ! y-direction
      trafomat_dpy(7,7) = 0.0_wp
      trafomat_dpy(8,7) = -SQRT3*SIN2T*SINP*0.5_wp
      trafomat_dpy(6,7) = SQRT3*SIN2T*COSP*0.5_wp
      trafomat_dpy(9,7) = SQRT3*SINT**2*DCOS2P*0.5_wp
      trafomat_dpy(5,7) = SQRT3*SINT**2*DSIN2P*0.5_wp
      trafomat_dpy(7,8) = 0.0_wp
      trafomat_dpy(8,8) = -COS2T*SINP
      trafomat_dpy(6,8) = COS2T*COSP
      trafomat_dpy(9,8) = SIN2T*DCOS2P*0.5_wp
      trafomat_dpy(5,8) = SIN2T*DSIN2P*0.5_wp
      trafomat_dpy(7,6) = 0.0_wp
      trafomat_dpy(8,6) = -COST*COSP
      trafomat_dpy(6,6) = -COST*SINP
      trafomat_dpy(9,6) = -SINT*DSIN2P
      trafomat_dpy(5,6) = SINT*DCOS2P
      trafomat_dpy(7,9) = 0.0_wp
      trafomat_dpy(8,9) = SIN2T*SINP*0.5_wp
      trafomat_dpy(6,9) = -SIN2T*COSP*0.5_wp
      trafomat_dpy(9,9) = (1.0_wp + COST**2)*DCOS2P*0.5_wp
      trafomat_dpy(5,9) = (1.0_wp + COST**2)*DSIN2P*0.5_wp
      trafomat_dpy(7,5) = 0.0_wp
      trafomat_dpy(8,5) = SINT*COSP
      trafomat_dpy(6,5) = SINT*SINP
      trafomat_dpy(9,5) = -COST*DSIN2P
      trafomat_dpy(5,5) = COST*DCOS2P

      if ( maxl <= 2 ) then 
         ! Transform to cartesian coordinates
         dtrafomat(1, 1:9, 1:9) = dpdx * trafomat_dp(1:9, 1:9) + dtdx * trafomat_dt(1:9, 1:9)
         dtrafomat(2, 1:9, 1:9) = dpdy * trafomat_dpy(1:9, 1:9) + dtdy * trafomat_dty(1:9, 1:9)
         dtrafomat(3, 1:9, 1:9) = dtdz * trafomat_dt(1:9, 1:9)
         return
      end if

      ! -------------------------------
      !  f functions (trafomat(16x16)) 
      ! -------------------------------

      ! f-functions are not transformed to the diatomic frame!
      ! Hence, the derivative remains 0!
      trafomat(:,10,10) = 1.0_wp
      trafomat(:,11,11) = 1.0_wp
      trafomat(:,12,12) = 1.0_wp
      trafomat(:,13,13) = 1.0_wp
      trafomat(:,14,14) = 1.0_wp
      trafomat(:,15,15) = 1.0_wp
      trafomat(:,16,16) = 1.0_wp

      if ( maxl <= 3 ) then 
         ! Transform to cartesian coordinates
         dtrafomat(1, 1:16, 1:16) = dpdx * trafomat_dp(1:16, 1:16) + dtdx * trafomat_dt(1:16, 1:16)
         dtrafomat(2, 1:16, 1:16) = dpdy * trafomat_dpy(1:16, 1:16) + dtdy * trafomat_dty(1:16, 1:16)
         dtrafomat(3, 1:16, 1:16) = dtdz * trafomat_dt(1:16, 1:16)
         return
      end if

      ! -------------------------------
      !  g functions (trafomat(25x25)) 
      ! -------------------------------

      ! g-functions are not transformed to the diatomic frame!
      ! Hence, the derivative remains 0!
      trafomat(:,17,17) = 1.0_wp
      trafomat(:,18,18) = 1.0_wp
      trafomat(:,19,19) = 1.0_wp
      trafomat(:,20,20) = 1.0_wp
      trafomat(:,21,21) = 1.0_wp
      trafomat(:,22,22) = 1.0_wp
      trafomat(:,23,23) = 1.0_wp
      trafomat(:,24,24) = 1.0_wp
      trafomat(:,25,25) = 1.0_wp

      if ( maxl <= 4 ) then 
         ! Transform to cartesian coordinates
         dtrafomat(1, 1:25, 1:25) = dpdx * trafomat_dp(1:25, 1:25) + dtdx * trafomat_dt(1:25, 1:25)
         dtrafomat(2, 1:25, 1:25) = dpdy * trafomat_dpy(1:25, 1:25) + dtdy * trafomat_dty(1:25, 1:25)
         dtrafomat(3, 1:25, 1:25) = dtdz * trafomat_dt(1:25, 1:25)
         return
      end if

   end subroutine d_harmtr

   pure subroutine scale_diatomic_frame(diat_mat, ksig, kpi, kdel, maxlj, maxli)
      !> Block matrix in the diatomic frame to be scaled
      real(wp),intent(inout)    :: diat_mat(:,:)
      !> Scaling parameters for different bonding contributions
      real(wp),intent(in)       :: ksig, kpi, kdel
      !> Highest angular momentum of atom j (first index)
      integer,intent(in)        :: maxlj
      !> Highest angular momentum of atom i (second index)
      integer,intent(in)        :: maxli

      integer :: maxl

      maxl = max(maxlj, maxli)

      diat_mat(1,1) = diat_mat(1,1)*ksig ! Sigma bond s   <-> s
      if(maxlj > 0) then
         diat_mat(3,1) = diat_mat(3,1)*ksig ! Sigma bond pz  <-> s 
      end if
      if(maxli > 0)  then
         diat_mat(1,3) = diat_mat(1,3)*ksig ! Sigma bond s   <-> pz
      end if
      if(maxlj > 0 .and. maxli > 0) then 
         diat_mat(3,3) = diat_mat(3,3)*ksig ! Sigma bond pz  <-> pz
         diat_mat(4,4) = diat_mat(4,4)*kpi  ! Pi    bond px  <-> px
         diat_mat(2,2) = diat_mat(2,2)*kpi  ! Pi    bond py  <-> py
         if(maxlj > 1) then
            diat_mat(7,3) = diat_mat(7,3)*ksig ! Sigma bond dz2 <-> pz
            diat_mat(8,4) = diat_mat(8,4)*kpi  ! Pi    bond dxz <-> px
            diat_mat(6,2) = diat_mat(6,2)*kpi  ! Pi    bond dyz <-> py
         end if
         if(maxli > 1) then         
            diat_mat(3,7) = diat_mat(3,7)*ksig ! Sigma bond pz  <-> dz2
            diat_mat(4,8) = diat_mat(4,8)*kpi  ! Pi    bond px  <-> dxz
            diat_mat(2,6) = diat_mat(2,6)*kpi  ! Pi    bond py  <-> dyz
         end if
      end if
      if (maxlj > 1) then
         diat_mat(7,1) = diat_mat(7,1)*ksig ! Sigma bond dz2 <-> s
      end if
      if (maxli > 1) then
         diat_mat(1,7) = diat_mat(1,7)*ksig ! Sigma bond s   <-> dz2   
      end if
      if (maxlj > 1 .and. maxli > 1) then
         diat_mat(7,7) = diat_mat(7,7)*ksig ! Sigma bond dz2 <-> dz2
         diat_mat(8,8) = diat_mat(8,8)*kpi  ! Pi    bond dxz <-> dxz
         diat_mat(6,6) = diat_mat(6,6)*kpi  ! Pi    bond dyz <-> dyz
         diat_mat(9,9) = diat_mat(9,9)*kdel ! Delta bond dx2-y2 <-> dx2-y2
         diat_mat(5,5) = diat_mat(5,5)*kdel ! Delta bond dxy <-> dxy
      endif
      ! f- and g-functions remain unscaled

   end subroutine scale_diatomic_frame

end module tblite_integral_diat_trafo
