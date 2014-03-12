!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \note
!> This module contains routines necessary to operate on plane waves on NVIDIA
!> GPUs using CUDA. It depends at execution time on NVIDIA's CUFFT library.
!> \par History
!>      BGL (06-Mar-2008)  : Created
!>      AG  (18-May-2012)  : Refacturing:
!>                           - added explicit interfaces to C routines
!>                           - enable double precision complex transformations
!>      AG  (11-Sept-2012) : Modifications:
!>                          - use pointers if precision mapping is not required
!>                          - use OMP for mapping
!> \author Benjamin G. Levine
! *****************************************************************************
MODULE pw_methods_cuda
#if !defined (__HAS_NO_ISO_C_BINDING)
  USE ISO_C_BINDING
#endif
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE pw_types,                        ONLY: RECIPROCALSPACE,&
                                             REALSPACE,&
                                             pw_type
  USE pw_grid_types,                   ONLY: FULLSPACE
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE fft_tools,                       ONLY: get_fft_scratch,&
                                             release_fft_scratch,&
                                             fft_scratch_type,&
                                             fft_scratch_sizes,&
                                             cube_transpose_1,&
                                             cube_transpose_2,&
                                             yz_to_x, x_to_yz,&
                                             xz_to_yz, yz_to_xz
  USE message_passing,                 ONLY: mp_comm_compare,&
                                             mp_environ,&
                                             mp_rank_compare
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

#include "cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pw_cuda_r3dc1d_3d
  PUBLIC :: pw_cuda_c1dr3d_3d
  PUBLIC :: pw_cuda_r3dc1d_3d_ps
  PUBLIC :: pw_cuda_c1dr3d_3d_ps

 ! Explicit interfaces to double precision complex transformation
 ! routines (C+CUDA). For more details, see: ../cuda_tools!
  INTERFACE pw_cuda_cfffg_cu
    SUBROUTINE pw_cuda_cfffg_z(din, zout, ghatmap, npts, ngpts, scale)&
      BIND(C, name="pw_cuda_cfffg_z_")
      USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: din
    TYPE(C_PTR), VALUE                       :: zout
    TYPE(C_PTR), INTENT(IN), VALUE           :: ghatmap
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: ngpts
    REAL(KIND=C_DOUBLE), INTENT(IN), VALUE   :: scale

    END SUBROUTINE pw_cuda_cfffg_z
  END INTERFACE

  INTERFACE pw_cuda_sfffc_cu
    SUBROUTINE pw_cuda_sfffc_z(zin, dout, ghatmap, npts, ngpts, nmaps, scale)&
      BIND(C, name="pw_cuda_sfffc_z_")
      USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: zin
    TYPE(C_PTR), VALUE                       :: dout
    TYPE(C_PTR), INTENT(IN), VALUE           :: ghatmap
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: ngpts, nmaps
    REAL(KIND=C_DOUBLE), INTENT(IN), VALUE   :: scale

    END SUBROUTINE pw_cuda_sfffc_z
  END INTERFACE

  INTERFACE pw_cuda_cff_cu
    SUBROUTINE pw_cuda_cff_z(din, zout, npts)&
      BIND(C, name="pw_cuda_cff_z_")
      USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: din
    TYPE(C_PTR), VALUE                       :: zout
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts

    END SUBROUTINE pw_cuda_cff_z
  END INTERFACE

  INTERFACE pw_cuda_ffc_cu
    SUBROUTINE pw_cuda_ffc_z(zin, dout, npts)&
      BIND(C, name="pw_cuda_ffc_z_")
      USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: zin
    TYPE(C_PTR), VALUE                       :: dout
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts

    END SUBROUTINE pw_cuda_ffc_z
  END INTERFACE

  INTERFACE pw_cuda_cf_cu
    SUBROUTINE pw_cuda_cf_z(din, zout, npts)&
      BIND(C, name="pw_cuda_cf_z_")
      USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: din
    TYPE(C_PTR), VALUE                       :: zout
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts

    END SUBROUTINE pw_cuda_cf_z
  END INTERFACE

  INTERFACE pw_cuda_fc_cu
    SUBROUTINE pw_cuda_fc_z(zin, dout, npts)&
      BIND(C, name="pw_cuda_fc_z_")
      USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: zin
    TYPE(C_PTR), VALUE                       :: dout
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts

    END SUBROUTINE pw_cuda_fc_z
  END INTERFACE

  INTERFACE pw_cuda_f_cu
    SUBROUTINE pw_cuda_f_z(zin, zout, dir, n, m)&
      BIND(C, name="pw_cuda_f_z_")
      USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: zin
    TYPE(C_PTR), VALUE                       :: zout
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: dir, n, m

    END SUBROUTINE pw_cuda_f_z
  END INTERFACE

  INTERFACE pw_cuda_fg_cu
    SUBROUTINE pw_cuda_fg_z(zin, zout, ghatmap, npts, mmax, ngpts, scale)&
      BIND(C, name="pw_cuda_fg_z_")
      USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: zin
    TYPE(C_PTR), VALUE                       :: zout
    TYPE(C_PTR), INTENT(IN), VALUE           :: ghatmap
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: mmax, ngpts
    REAL(KIND=C_DOUBLE), INTENT(IN), VALUE   :: scale

    END SUBROUTINE pw_cuda_fg_z
  END INTERFACE

  INTERFACE pw_cuda_sf_cu
    SUBROUTINE pw_cuda_sf_z(zin, zout, ghatmap, npts, mmax, ngpts, nmaps, scale)&
      BIND(C, name="pw_cuda_sf_z_")
      USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: zin
    TYPE(C_PTR), VALUE                       :: zout
    TYPE(C_PTR), INTENT(IN), VALUE           :: ghatmap
    INTEGER(KIND=C_INT), DIMENSION(*), &
      INTENT(IN)                             :: npts
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: mmax, ngpts, nmaps
    REAL(KIND=C_DOUBLE), INTENT(IN), VALUE   :: scale

    END SUBROUTINE pw_cuda_sf_z
  END INTERFACE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_methods_cuda'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

CONTAINS
! *****************************************************************************
!> \brief perform an fft followed by a gather on the gpu
!> \author Benjamin G Levine
! *****************************************************************************
  SUBROUTINE pw_cuda_r3dc1d_3d(pw1, pw2, scale, error)
    TYPE(pw_type), TARGET, INTENT(IN)        :: pw1
    TYPE(pw_type), TARGET, INTENT(INOUT)     :: pw2
    REAL(KIND=dp)                            :: scale
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_r3dc1d_3d', &
      routineP = moduleN//':'//routineN

#if defined (__PW_CUDA)
    INTEGER                                  :: handle, ngpts
    INTEGER                                  :: l1, l2, l3
    INTEGER, DIMENSION(:), POINTER           :: npts

    REAL(KIND=dp), POINTER                   :: ptr_pwin
    COMPLEX(KIND=dp), POINTER                :: ptr_pwout
    INTEGER, POINTER                         :: ptr_ghatmap
    CALL timeset(routineN,handle)

    ngpts = SIZE(pw2%pw_grid%gsq)
    l1 = LBOUND(pw1%cr3d,1)
    l2 = LBOUND(pw1%cr3d,2)
    l3 = LBOUND(pw1%cr3d,3)
    npts => pw1%pw_grid%npts

   ! pointers to data arrays
    ptr_pwin => pw1%cr3d(l1,l2,l3)
    ptr_pwout => pw2%cc(1)

   ! pointer to map array
    ptr_ghatmap => pw2%pw_grid%g_hatmap(1,1)

   ! invoke the combined transformation
    CALL pw_cuda_cfffg_cu(c_loc(ptr_pwin), c_loc(ptr_pwout), c_loc(ptr_ghatmap), npts, ngpts, scale)

    pw2 % in_space = RECIPROCALSPACE

    CALL timestop(handle)
#endif
  END SUBROUTINE pw_cuda_r3dc1d_3d

! *****************************************************************************
!> \brief perform an scatter followed by a fft on the gpu
!> \author Benjamin G Levine
! *****************************************************************************
  SUBROUTINE pw_cuda_c1dr3d_3d(pw1, pw2, scale, error)
    TYPE(pw_type), TARGET, INTENT(IN)        :: pw1
    TYPE(pw_type), TARGET, INTENT(INOUT)     :: pw2
    REAL(KIND=dp)                            :: scale
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_c1dr3d_3d', &
      routineP = moduleN//':'//routineN

#if defined (__PW_CUDA)
    INTEGER                                  :: handle, ngpts, nmaps
    INTEGER                                  :: l1, l2, l3
    INTEGER, DIMENSION(:), POINTER           :: npts

    COMPLEX(KIND=dp), POINTER                :: ptr_pwin
    REAL(KIND=dp), POINTER                   :: ptr_pwout
    INTEGER, POINTER                         :: ptr_ghatmap

    CALL timeset(routineN,handle)

    ngpts = SIZE(pw1%pw_grid%gsq)
    l1 = LBOUND(pw2%cr3d,1)
    l2 = LBOUND(pw2%cr3d,2)
    l3 = LBOUND(pw2%cr3d,3)
    npts => pw1%pw_grid%npts

   ! pointers to data arrays
    ptr_pwin => pw1%cc(1)
    ptr_pwout => pw2%cr3d(l1,l2,l3)

   ! pointer to map array
    nmaps = SIZE(pw1%pw_grid%g_hatmap,2)
    ptr_ghatmap => pw1%pw_grid%g_hatmap(1,1)

   ! invoke the combined transformation
    CALL pw_cuda_sfffc_cu(c_loc(ptr_pwin), c_loc(ptr_pwout), c_loc(ptr_ghatmap), npts, ngpts, nmaps, scale)

    pw2 % in_space = REALSPACE

    CALL timestop(handle)
#endif
  END SUBROUTINE pw_cuda_c1dr3d_3d

! *****************************************************************************
!> \brief perform an parallel fft followed by a gather on the gpu
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_r3dc1d_3d_ps(pw1, pw2, scale, error)
    TYPE(pw_type), TARGET, INTENT(IN)        :: pw1
    TYPE(pw_type), TARGET, INTENT(INOUT)     :: pw2
    REAL(KIND=dp)                            :: scale
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_r3dc1d_3d_ps', &
      routineP = moduleN//':'//routineN

#if defined (__PW_CUDA)
    INTEGER                                  :: handle, ierr
    INTEGER                                  :: i, ngpts, stat, iout
    INTEGER                                  :: lg, mg, mmax, lmax, rp, mx1, my1, mx2, mz2, n1, n2, nmax
    INTEGER                                  :: g_pos, numtask, numtask_r, numtask_g
    INTEGER, DIMENSION(2)                    :: r_pos, r_dim
    INTEGER                                  :: gs_group, rs_group
    INTEGER, DIMENSION(:), POINTER           :: n, nloc
    INTEGER, DIMENSION(:), POINTER           :: nyzray !nyzray(0:)
    INTEGER, DIMENSION(:,:,:), POINTER       :: yzp    !yzp(:,:,0:)
    INTEGER, DIMENSION(:,:,:,:), POINTER     :: bo     !bo(:,:,0:,:)
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      POINTER                                :: grays
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: tbuf
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: pbuf, rbuf, sbuf, qbuf
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: p2p
    TYPE(fft_scratch_sizes)                  :: fft_scratch_size
    TYPE(fft_scratch_type), POINTER          :: fft_scratch
    
    CALL timeset(routineN,handle)

   ! dimensions
    n => pw1%pw_grid%npts
    nloc => pw1%pw_grid%npts_local
    grays => pw1%pw_grid%grays
    ngpts = nloc(1) * nloc(2) * nloc(3)

    !..transform
    IF ( pw1%pw_grid%para%ray_distribution ) THEN
       gs_group =  pw1%pw_grid%para%group
       rs_group =  pw1%pw_grid%para%rs_group
       yzp      => pw1%pw_grid%para%yzp
       nyzray   => pw1%pw_grid%para%nyzray
       bo       => pw1%pw_grid%para%bo

       CALL mp_environ(numtask_g, g_pos, gs_group)
       CALL mp_environ(numtask_r, r_dim, r_pos, rs_group)
       IF ( numtask_g /= numtask_r ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                         "Real space and G space groups are different.")
       END IF
       numtask = numtask_r
       CALL mp_comm_compare(rs_group, gs_group, iout)
       IF ( iout >3 ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                         "Real space and G space groups are different.")
       END IF

       lg   = SIZE(grays, 1)
       mg   = SIZE(grays, 2)
       mmax = MAX(mg, 1)
       lmax = MAX(lg, (ngpts / mmax + 1))

       ALLOCATE(p2p(0:numtask - 1), STAT = ierr)
       IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "p2p",int_size*numtask)

       CALL mp_rank_compare(gs_group, rs_group, p2p)

       rp   = p2p(g_pos)
       mx2  = bo(2,1,rp,2) - bo(1,1,rp,2) + 1
       mz2  = bo(2,3,rp,2) - bo(1,3,rp,2) + 1
       n1   = MAXVAL(bo(2,1,:,1) - bo(1,1,:,1) + 1)
       n2   = MAXVAL(bo(2,2,:,1) - bo(1,2,:,1) + 1)
       nmax = MAX((2*n2)/numtask, 2) * mx2*mz2
       nmax = MAX(nmax, n1*MAXVAL(nyzray))

       fft_scratch_size%nx        = nloc(1)
       fft_scratch_size%ny        = nloc(2)
       fft_scratch_size%nz        = nloc(3)
       fft_scratch_size%lmax      = lmax
       fft_scratch_size%mmax      = mmax
       fft_scratch_size%mx1       = bo(2,1,rp,1) - bo(1,1,rp,1) + 1
       fft_scratch_size%mx2       = mx2
       fft_scratch_size%my1       = bo(2,2,rp,1) - bo(1,2,rp,1) + 1
       fft_scratch_size%mz2       = mz2
       fft_scratch_size%lg        = lg
       fft_scratch_size%mg        = mg
       fft_scratch_size%nbx       = MAXVAL(bo(2,1,:,2))
       fft_scratch_size%nbz       = MAXVAL(bo(2,3,:,2))
       fft_scratch_size%mcz1      = MAXVAL(bo(2,3,:,1) - bo(1,3,:,1) + 1)
       fft_scratch_size%mcx2      = MAXVAL(bo(2,1,:,2) - bo(1,1,:,2) + 1)
       fft_scratch_size%mcz2      = MAXVAL(bo(2,3,:,2) - bo(1,3,:,2) + 1)
       fft_scratch_size%nmax      = nmax
       fft_scratch_size%nmray     = MAXVAL(nyzray)
       fft_scratch_size%nyzray    = nyzray(g_pos)
       fft_scratch_size%gs_group  = gs_group
       fft_scratch_size%rs_group  = rs_group
       fft_scratch_size%g_pos     = g_pos
       fft_scratch_size%r_pos     = r_pos
       fft_scratch_size%r_dim     = r_dim
       fft_scratch_size%numtask   = numtask

       IF (r_dim(2) > 1) THEN
          !
          ! real space is distributed over x and y coordinate
          ! we have two stages of communication
          !
          IF (r_dim(1) == 1) &
             CALL stop_program(routineN,moduleN,__LINE__,&
                  "This processor distribution is not supported.")

          CALL get_fft_scratch(fft_scratch, tf_type = 300, n = n, fft_sizes = fft_scratch_size, error = error)

          ! assign buffers
          qbuf => fft_scratch%p2buf
          rbuf => fft_scratch%p3buf
          pbuf => fft_scratch%p4buf
          sbuf => fft_scratch%p5buf

          ! FFT along z
          CALL pw_cuda_cf(pw1, qbuf, error)

          ! Exchange data ( transpose of matrix )
          CALL cube_transpose_2(qbuf, rs_group, bo(:,:,:,1), bo(:,:,:,2), rbuf, fft_scratch, error)

          ! FFT along y
          ! use the inbuild fft-lib
          ! CALL fft_1dm(fft_scratch%fft_plan(2), rbuf, pbuf, 1.0_dp, stat)
          ! or cufft (works faster, but is only faster if plans are stored)
          CALL pw_cuda_f(rbuf, pbuf, +1, n(2), mx2*mz2, error)

          ! Exchange data ( transpose of matrix ) and sort
          CALL xz_to_yz(pbuf, rs_group, r_dim, g_pos, p2p, yzp, nyzray, &
               bo(:,:,:,2), sbuf, fft_scratch, error)

          ! FFT along x
          CALL pw_cuda_fg(sbuf, pw2, scale, error)

          CALL release_fft_scratch(fft_scratch, error)

       ELSE
          !
          ! real space is only distributed over x coordinate
          ! we have one stage of communication, after the transform of
          ! direction x
          !

          CALL get_fft_scratch(fft_scratch, tf_type = 200, n = n, fft_sizes = fft_scratch_size, error = error)

          ! assign buffers
          tbuf => fft_scratch%tbuf
          sbuf => fft_scratch%r1buf

          ! FFT along y and z
          CALL pw_cuda_cff(pw1, tbuf, error)

          ! Exchange data ( transpose of matrix ) and sort
          CALL yz_to_x(tbuf, gs_group, g_pos, p2p, yzp, nyzray, &
               bo(:,:,:,2), sbuf, fft_scratch, error)

          ! FFT along x
          CALL pw_cuda_fg(sbuf, pw2, scale, error)

          CALL release_fft_scratch(fft_scratch,error)

       ENDIF

       DEALLOCATE ( p2p, STAT = ierr )
       IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"p2p")

!--------------------------------------------------------------------------
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
          "Not implemented (no ray_distr.) in: pw_cuda_c1dr3d_3d_ps.")
       !CALL fft3d ( dir, n, pwin, grays, pw1%pw_grid%para%rs_group, &
       !     pw1%pw_grid%para%bo, scale = scale, debug=test )
    END IF

    pw2 % in_space = RECIPROCALSPACE

    CALL timestop(handle)
#endif
  END SUBROUTINE pw_cuda_r3dc1d_3d_ps

! *****************************************************************************
!> \brief perform an parallel scatter followed by a fft on the gpu
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_c1dr3d_3d_ps(pw1, pw2, scale, error)
    TYPE(pw_type), TARGET, INTENT(IN)        :: pw1
    TYPE(pw_type), TARGET, INTENT(INOUT)     :: pw2
    REAL(KIND=dp)                            :: scale
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_c1dr3d_3d_ps', &
      routineP = moduleN//':'//routineN

#if defined (__PW_CUDA)
    INTEGER                                  :: handle, ierr

    INTEGER                                  :: i, ngpts, stat, iout
    INTEGER                                  :: lg, mg, mmax, lmax, rp, mx1, my1, mx2, mz2, n1, n2, nmax, mcz1, mcx2, mcz2
    INTEGER                                  :: g_pos, numtask, numtask_r, numtask_g
    INTEGER, DIMENSION(2)                    :: r_pos, r_dim
    INTEGER                                  :: gs_group, rs_group
    INTEGER, DIMENSION(:), POINTER           :: n, nloc
    INTEGER, DIMENSION(:), POINTER           :: nyzray !nyzray(0:)
    INTEGER, DIMENSION(:,:,:), POINTER       :: yzp    !yzp(:,:,0:)
    INTEGER, DIMENSION(:,:,:,:), POINTER     :: bo     !bo(:,:,0:,:)
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      POINTER                                :: grays
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: tbuf
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: pbuf, rbuf, sbuf, qbuf
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: p2p
    TYPE(fft_scratch_sizes)                  :: fft_scratch_size
    TYPE(fft_scratch_type), POINTER          :: fft_scratch
    
    CALL timeset(routineN,handle)

    ! dimensions
    n => pw1%pw_grid%npts
    nloc => pw1%pw_grid%npts_local
    grays => pw1%pw_grid%grays
    ngpts = nloc(1) * nloc(2) * nloc(3)

    !..transform
    IF ( pw1%pw_grid%para%ray_distribution ) THEN
       gs_group =  pw1%pw_grid%para%group
       rs_group =  pw1%pw_grid%para%rs_group
       yzp      => pw1%pw_grid%para%yzp
       nyzray   => pw1%pw_grid%para%nyzray
       bo       => pw1%pw_grid%para%bo

       CALL mp_environ(numtask_g, g_pos, gs_group)
       CALL mp_environ(numtask_r, r_dim, r_pos, rs_group)
       IF ( numtask_g /= numtask_r ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Real space and G space groups are different.")
       END IF
       numtask = numtask_r
       CALL mp_comm_compare(rs_group, gs_group, iout)
       IF ( iout >3 ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                         "Real space and G space groups are different.")
       END IF

       lg   = SIZE(grays, 1)
       mg   = SIZE(grays, 2)
       mmax = MAX(mg, 1)
       lmax = MAX(lg, (ngpts / mmax + 1))

       ALLOCATE(p2p(0:numtask - 1), STAT = ierr)
       IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "p2p",int_size*numtask)

       CALL mp_rank_compare(gs_group, rs_group, p2p)

       rp   = p2p(g_pos)
       mx2  = bo(2,1,rp,2) - bo(1,1,rp,2) + 1
       mz2  = bo(2,3,rp,2) - bo(1,3,rp,2) + 1
       n1   = MAXVAL(bo(2,1,:,1) - bo(1,1,:,1) + 1)
       n2   = MAXVAL(bo(2,2,:,1) - bo(1,2,:,1) + 1)
       nmax = MAX((2*n2)/numtask, 2) * mx2*mz2
       nmax = MAX(nmax, n1*MAXVAL(nyzray))

       fft_scratch_size%nx       = nloc(1)
       fft_scratch_size%ny       = nloc(2)
       fft_scratch_size%nz       = nloc(3)
       fft_scratch_size%lmax     = lmax
       fft_scratch_size%mmax     = mmax
       fft_scratch_size%mx1      = bo(2,1,rp,1) - bo(1,1,rp,1) + 1
       fft_scratch_size%mx2      = mx2
       fft_scratch_size%my1      = bo(2,2,rp,1) - bo(1,2,rp,1) + 1
       fft_scratch_size%mz2      = mz2
       fft_scratch_size%lg       = lg
       fft_scratch_size%mg       = mg
       fft_scratch_size%nbx      = MAXVAL(bo(2,1,:,2))
       fft_scratch_size%nbz      = MAXVAL(bo(2,3,:,2))
       fft_scratch_size%mcz1     = MAXVAL(bo(2,3,:,1) - bo(1,3,:,1) + 1)
       fft_scratch_size%mcx2     = MAXVAL(bo(2,1,:,2) - bo(1,1,:,2) + 1)
       fft_scratch_size%mcz2     = MAXVAL(bo(2,3,:,2) - bo(1,3,:,2) + 1)
       fft_scratch_size%nmax     = nmax
       fft_scratch_size%nmray    = MAXVAL(nyzray)
       fft_scratch_size%nyzray   = nyzray(g_pos)
       fft_scratch_size%gs_group = gs_group
       fft_scratch_size%rs_group = rs_group
       fft_scratch_size%g_pos    = g_pos
       fft_scratch_size%r_pos    = r_pos
       fft_scratch_size%r_dim    = r_dim
       fft_scratch_size%numtask  = numtask

       IF (r_dim(2) > 1) THEN
          !
          ! real space is distributed over x and y coordinate
          ! we have two stages of communication
          !
          IF (r_dim(1) == 1) &
             CALL stop_program(routineN,moduleN,__LINE__,&
                  "This processor distribution is not supported.")

          CALL get_fft_scratch(fft_scratch, tf_type = 300, n = n, fft_sizes = fft_scratch_size, error = error)

          ! assign buffers
          pbuf => fft_scratch%p7buf
          qbuf => fft_scratch%p4buf
          rbuf => fft_scratch%p3buf
          sbuf => fft_scratch%p2buf

          ! FFT along x
          CALL pw_cuda_sf(pw1, pbuf, scale, error)

          ! Exchange data ( transpose of matrix ) and sort
          IF (pw1%pw_grid%grid_span /= FULLSPACE) CALL zero_c(SIZE(qbuf),qbuf)
          CALL yz_to_xz(pbuf, rs_group, r_dim, g_pos, p2p, yzp, nyzray, &
               bo(:,:,:,2), qbuf, fft_scratch, error)

          ! FFT along y
          ! use the inbuild fft-lib
          ! CALL fft_1dm(fft_scratch%fft_plan(5), qbuf, rbuf, 1.0_dp, stat)
          ! or cufft (works faster, but is only faster if plans are stored)
          CALL pw_cuda_f(qbuf, rbuf, -1, n(2), mx2*mz2, error)

          ! Exchange data ( transpose of matrix )
          IF (pw1%pw_grid%grid_span /= FULLSPACE) CALL zero_c(SIZE(sbuf),sbuf)
          CALL cube_transpose_1(rbuf, rs_group, bo(:,:,:,2), bo(:,:,:,1), sbuf, fft_scratch, error)

          ! FFT along z
          CALL pw_cuda_fc(sbuf, pw2, error)

          CALL release_fft_scratch(fft_scratch, error)

       ELSE
          !
          ! real space is only distributed over x coordinate
          ! we have one stage of communication, after the transform of
          ! direction x
          !

          CALL get_fft_scratch(fft_scratch, tf_type = 200, n = n, fft_sizes = fft_scratch_size, error = error)

          ! assign buffers
          sbuf => fft_scratch%r1buf
          tbuf => fft_scratch%tbuf

          ! FFT along x
          CALL pw_cuda_sf(pw1, sbuf, scale, error)

          ! Exchange data ( transpose of matrix ) and sort
          IF (pw1%pw_grid%grid_span /= FULLSPACE) CALL zero_c(SIZE(tbuf),tbuf)
          CALL x_to_yz (sbuf, gs_group, g_pos, p2p, yzp, nyzray, &
               bo(:,:,:,2), tbuf, fft_scratch, error)

          ! FFT along y and z
          CALL pw_cuda_ffc(tbuf, pw2, error)

          CALL release_fft_scratch(fft_scratch,error)

       ENDIF

       DEALLOCATE ( p2p, STAT = ierr )
       IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"p2p")

!--------------------------------------------------------------------------
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
          "Not implemented (no ray_distr.) in: pw_cuda_c1dr3d_3d_ps.")
       !CALL fft3d ( dir, n, pwin, grays, pw1%pw_grid%para%rs_group, &
       !     pw1%pw_grid%para%bo, scale = scale, debug=test )
    END IF

    pw2 % in_space = REALSPACE

    CALL timestop(handle)
#endif
  END SUBROUTINE pw_cuda_c1dr3d_3d_ps

! *****************************************************************************
!> \brief perform a parallel real_to_complex copy followed by a 2D-FFT on the gpu
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_cff (pw1, pwbuf, error)
    TYPE(pw_type), TARGET, INTENT(IN)        :: pw1
    COMPLEX(KIND=dp), DIMENSION(:,:,:), &
      POINTER, INTENT(INOUT)                 :: pwbuf
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_cff', &
      routineP = moduleN//':'//routineN

#if defined (__PW_CUDA)
    INTEGER                                  :: l1, l2, l3, handle
    INTEGER, DIMENSION(:), POINTER           :: npts
    REAL(KIND=dp), POINTER                   :: ptr_pwin
    COMPLEX(KIND=dp), POINTER                :: ptr_pwout

    CALL timeset(routineN,handle)

   ! dimensions
    npts => pw1%pw_grid%npts_local
    l1 = LBOUND(pw1%cr3d,1)
    l2 = LBOUND(pw1%cr3d,2)
    l3 = LBOUND(pw1%cr3d,3)

   ! pointers to data arrays
    ptr_pwin => pw1%cr3d(l1,l2,l3)
    ptr_pwout => pwbuf(1,1,1)

   ! invoke the combined transformation
    CALL pw_cuda_cff_cu(c_loc(ptr_pwin), c_loc(ptr_pwout), npts)

    CALL timestop(handle)
#endif
  END SUBROUTINE pw_cuda_cff

! *****************************************************************************
!> \brief perform a parallel 2D-FFT followed by a complex_to_real copy on the gpu
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_ffc (pwbuf, pw2, error)
    COMPLEX(KIND=dp), DIMENSION(:,:,:), &
      POINTER, INTENT(IN)                    :: pwbuf
    TYPE(pw_type), TARGET, INTENT(INOUT)     :: pw2
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_ffc', &
      routineP = moduleN//':'//routineN

#if defined (__PW_CUDA)
    INTEGER                                  :: l1, l2, l3, handle
    INTEGER, DIMENSION(:), POINTER           :: npts
    COMPLEX(KIND=dp), POINTER                :: ptr_pwin
    REAL(KIND=dp), POINTER                   :: ptr_pwout

    CALL timeset(routineN,handle)

   ! dimensions
    npts => pw2%pw_grid%npts_local
    l1 = LBOUND(pw2%cr3d,1)
    l2 = LBOUND(pw2%cr3d,2)
    l3 = LBOUND(pw2%cr3d,3)

   ! pointers to data arrays
    ptr_pwin => pwbuf(1,1,1)
    ptr_pwout => pw2%cr3d(l1,l2,l3)

   ! invoke the combined transformation
    CALL pw_cuda_ffc_cu(c_loc(ptr_pwin), c_loc(ptr_pwout), npts)

    CALL timestop(handle)
#endif
  END SUBROUTINE pw_cuda_ffc

! *****************************************************************************
!> \brief perform a parallel real_to_complex copy followed by a 1D-FFT on the gpu
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_cf (pw1, pwbuf, error)
    TYPE(pw_type), TARGET, INTENT(IN)        :: pw1
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      POINTER, INTENT(INOUT)                 :: pwbuf
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_cf', &
      routineP = moduleN//':'//routineN

#if defined (__PW_CUDA)
    INTEGER                                  :: l1, l2, l3, handle
    INTEGER, DIMENSION(:), POINTER           :: npts
    REAL(KIND=dp), POINTER                   :: ptr_pwin
    COMPLEX(KIND=dp), POINTER                :: ptr_pwout

    CALL timeset(routineN,handle)

   ! dimensions
    npts => pw1%pw_grid%npts_local
    l1 = LBOUND(pw1%cr3d,1)
    l2 = LBOUND(pw1%cr3d,2)
    l3 = LBOUND(pw1%cr3d,3)

   ! pointers to data arrays
    ptr_pwin => pw1%cr3d(l1,l2,l3)
    ptr_pwout => pwbuf(1,1)

   ! invoke the combined transformation
    CALL pw_cuda_cf_cu(c_loc(ptr_pwin), c_loc(ptr_pwout), npts)

    CALL timestop(handle)
#endif
  END SUBROUTINE pw_cuda_cf

! *****************************************************************************
!> \brief perform a parallel 1D-FFT followed by a complex_to_real copy on the gpu
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_fc (pwbuf, pw2, error)
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      POINTER, INTENT(IN)                    :: pwbuf
    TYPE(pw_type), TARGET, INTENT(INOUT)     :: pw2
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_fc', &
      routineP = moduleN//':'//routineN

#if defined (__PW_CUDA)
    INTEGER                                  :: l1, l2, l3, handle
    INTEGER, DIMENSION(:), POINTER           :: npts
    COMPLEX(KIND=dp), POINTER                :: ptr_pwin
    REAL(KIND=dp), POINTER                   :: ptr_pwout

    CALL timeset(routineN,handle)

    npts => pw2%pw_grid%npts_local
    l1 = LBOUND(pw2%cr3d,1)
    l2 = LBOUND(pw2%cr3d,2)
    l3 = LBOUND(pw2%cr3d,3)

   ! pointers to data arrays
    ptr_pwin => pwbuf(1,1)
    ptr_pwout => pw2%cr3d(l1,l2,l3)

   ! invoke the combined transformation
    CALL pw_cuda_fc_cu(c_loc(ptr_pwin), c_loc(ptr_pwout), npts)

    CALL timestop(handle)
#endif
  END SUBROUTINE pw_cuda_fc

! *****************************************************************************
!> \brief perform a parallel 1D-FFT on the gpu
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_f(pwbuf1, pwbuf2, dir, n, m, error)
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      POINTER, INTENT(IN)                    :: pwbuf1
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      POINTER, INTENT(INOUT)                 :: pwbuf2
    INTEGER, INTENT(IN)                      :: dir
    INTEGER, INTENT(IN)                      :: n
    INTEGER, INTENT(IN)                      :: m
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_f', &
      routineP = moduleN//':'//routineN

#if defined (__PW_CUDA)
    INTEGER                                  :: handle
    COMPLEX(KIND=dp), POINTER                :: ptr_pwin
    COMPLEX(KIND=dp), POINTER                :: ptr_pwout

    CALL timeset(routineN,handle)

   ! pointers to data arrays
    ptr_pwin => pwbuf1(1,1)
    ptr_pwout => pwbuf2(1,1)

   ! invoke the combined transformation
    CALL pw_cuda_f_cu(c_loc(ptr_pwin), c_loc(ptr_pwout), dir, n, m)

    CALL timestop(handle)
#endif
  END SUBROUTINE pw_cuda_f
! *****************************************************************************
!> \brief perform a parallel 1D-FFT followed by a gather on the gpu
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_fg (pwbuf, pw2, scale, error)
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      POINTER, INTENT(IN)                    :: pwbuf
    TYPE(pw_type), TARGET, INTENT(INOUT)     :: pw2
    REAL(KIND=dp), INTENT(IN)                :: scale
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_fg', &
      routineP = moduleN//':'//routineN

#if defined (__PW_CUDA)
    INTEGER                                  :: mmax, mg, handle, ngpts
    INTEGER, DIMENSION(:), POINTER           :: npts
    
    COMPLEX(KIND=dp), POINTER                :: ptr_pwin
    COMPLEX(KIND=dp), POINTER                :: ptr_pwout
    INTEGER, POINTER                         :: ptr_ghatmap

    CALL timeset(routineN,handle)

    ngpts = SIZE(pw2%pw_grid%gsq)
    npts => pw2%pw_grid%npts

    mg = SIZE(pw2%pw_grid%grays, 2)
    mmax = MAX(mg, 1)

   ! pointers to data arrays
    ptr_pwin => pwbuf(1,1)
    ptr_pwout => pw2%cc(1)

   ! pointer to map array
    ptr_ghatmap => pw2%pw_grid%g_hatmap(1,1)

   ! invoke the combined transformation
    CALL pw_cuda_fg_cu(c_loc(ptr_pwin), c_loc(ptr_pwout), c_loc(ptr_ghatmap), npts, mmax, ngpts, scale)
  
    CALL timestop(handle)
#endif
  END SUBROUTINE pw_cuda_fg

! *****************************************************************************
!> \brief perform a parallel scatter followed by a 1D-FFT on the gpu
!> \author Andreas Gloess
! *****************************************************************************
  SUBROUTINE pw_cuda_sf (pw1, pwbuf, scale, error)
    TYPE(pw_type), TARGET, INTENT(IN)        :: pw1
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      POINTER, INTENT(INOUT)                 :: pwbuf
    REAL(KIND=dp), INTENT(IN)                :: scale
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_cuda_sf', &
      routineP = moduleN//':'//routineN

#if defined (__PW_CUDA)
    INTEGER                                  :: mmax, mg, ngpts, handle, nmaps
    INTEGER, DIMENSION(:), POINTER           :: npts
    
    COMPLEX(KIND=dp), DIMENSION(:), &
      ALLOCATABLE, TARGET                    :: pwin
    COMPLEX(KIND=dp), DIMENSION(:,:), &
      ALLOCATABLE, TARGET                    :: pwout
    COMPLEX(KIND=dp), POINTER                :: ptr_pwin
    COMPLEX(KIND=dp), POINTER                :: ptr_pwout
    INTEGER, POINTER                         :: ptr_ghatmap

    CALL timeset(routineN,handle)

    ngpts = SIZE(pw1%pw_grid%gsq)
    npts => pw1%pw_grid%npts

    mg = SIZE(pw1%pw_grid%grays, 2)
    mmax = MAX(mg, 1)

   ! pointers to data arrays
    ptr_pwin => pw1%cc(1)
    ptr_pwout => pwbuf(1,1)

   ! pointer to map array
    nmaps = SIZE(pw1%pw_grid%g_hatmap,2)
    ptr_ghatmap => pw1%pw_grid%g_hatmap(1,1)

   ! invoke the combined transformation
    CALL pw_cuda_sf_cu(c_loc(ptr_pwin), c_loc(ptr_pwout), c_loc(ptr_ghatmap), npts, mmax, ngpts, nmaps, scale)
  
    CALL timestop(handle)
#endif
  END SUBROUTINE pw_cuda_sf
END MODULE pw_methods_cuda

