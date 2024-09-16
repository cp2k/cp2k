!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \note
!> This module contains routines necessary to operate on plane waves on GPUs
!     >  independently of the GPU platform.
!> \par History
!>      BGL (06-Mar-2008)  : Created
!>      AG  (18-May-2012)  : Refacturing:
!>                           - added explicit interfaces to C routines
!>                           - enable double precision complex transformations
!>      AG  (11-Sept-2012) : Modifications:
!>                          - use pointers if precision mapping is not required
!>                          - use OMP for mapping
!>      MT  (Jan 2022)     : Modifications
!>                          - use a generic interface for fft calls to GPUs
!>                          - Support both Nvidia and AMD GPUs. Other GPUs manufacturers
!>                            can be added easily.
!> \author Benjamin G. Levine
! **************************************************************************************************
MODULE pw_gpu
   USE ISO_C_BINDING,                   ONLY: C_DOUBLE,&
                                              C_INT,&
                                              C_LOC,&
                                              C_PTR
   USE kinds,                           ONLY: dp
   USE pw_grid_types,                   ONLY: pw_grid_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pw_gpu_sf, pw_gpu_f, pw_gpu_fc, pw_gpu_ffc, pw_gpu_cf, pw_gpu_cff, pw_gpu_fg
   PUBLIC :: pw_gpu_cfffg, pw_gpu_sfffc
   PUBLIC :: pw_gpu_init, pw_gpu_finalize

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_gpu'
   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief Allocates resources on the gpu device for gpu fft acceleration
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE pw_gpu_init()
      INTEGER                                            :: dummy
      INTERFACE
         SUBROUTINE pw_gpu_init_c() BIND(C, name="pw_gpu_init")
         END SUBROUTINE pw_gpu_init_c
      END INTERFACE

      MARK_USED(dummy) ! TODO: fix fpretty
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_init_c()
#else
      ! Nothing to do.
#endif
   END SUBROUTINE pw_gpu_init

! **************************************************************************************************
!> \brief Releases resources on the gpu device for gpu fft acceleration
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE pw_gpu_finalize()
      INTEGER                                            :: dummy
      INTERFACE
         SUBROUTINE pw_gpu_finalize_c() BIND(C, name="pw_gpu_finalize")
         END SUBROUTINE pw_gpu_finalize_c
      END INTERFACE

      MARK_USED(dummy) ! TODO: fix fpretty
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_finalize_c()
#else
      ! Nothing to do.
#endif
   END SUBROUTINE pw_gpu_finalize

! **************************************************************************************************
!> \brief perform an fft followed by a gather on the gpu
!> \param zin ...
!> \param zout ...
!> \param pw_grid_out ...
!> \author Benjamin G Levine
! **************************************************************************************************
   SUBROUTINE pw_gpu_cfffg(zin, zout, pw_grid_out)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         TARGET                                          :: zin
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(INOUT), &
         TARGET                                          :: zout
      TYPE(pw_grid_type), INTENT(IN), TARGET             :: pw_grid_out

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_cfffg'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwout
      INTEGER                                            :: handle, l1, l2, l3, ngpts
      INTEGER, DIMENSION(:), POINTER                     :: npts
      INTEGER, POINTER                                   :: ptr_ghatmap
      REAL(KIND=dp)                                      :: scale
      REAL(KIND=dp), POINTER                             :: ptr_pwin
      INTERFACE
         SUBROUTINE pw_gpu_cfffg_c(din, zout, ghatmap, npts, ngpts, scale) BIND(C, name="pw_gpu_cfffg")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: din
            TYPE(C_PTR), VALUE                           :: zout
            TYPE(C_PTR), INTENT(IN), VALUE               :: ghatmap
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
            INTEGER(KIND=C_INT), INTENT(IN), VALUE       :: ngpts
            REAL(KIND=C_DOUBLE), INTENT(IN), VALUE       :: scale

         END SUBROUTINE pw_gpu_cfffg_c
      END INTERFACE

      CALL timeset(routineN, handle)

      scale = 1.0_dp/REAL(pw_grid_out%ngpts, KIND=dp)

      ngpts = SIZE(pw_grid_out%gsq)
      l1 = LBOUND(zin, 1)
      l2 = LBOUND(zin, 2)
      l3 = LBOUND(zin, 3)
      npts => pw_grid_out%npts

      ! pointers to data arrays
      ptr_pwin => zin(l1, l2, l3)
      ptr_pwout => zout(1)

      ! pointer to map array
      ptr_ghatmap => pw_grid_out%g_hatmap(1, 1)

      ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_cfffg_c(c_loc(ptr_pwin), c_loc(ptr_pwout), c_loc(ptr_ghatmap), npts, ngpts, scale)
#else
      CPABORT("Compiled without pw offloading.")
#endif

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_cfffg

! **************************************************************************************************
!> \brief perform an scatter followed by a fft on the gpu
!> \param zin ...
!> \param zout ...
!> \param pw_grid_in ...
!> \author Benjamin G Levine
! **************************************************************************************************
   SUBROUTINE pw_gpu_sfffc(zin, zout, pw_grid_in)
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN), TARGET :: zin
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         TARGET                                          :: zout
      TYPE(pw_grid_type), INTENT(IN), TARGET             :: pw_grid_in

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_sfffc'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwin
      INTEGER                                            :: handle, l1, l2, l3, ngpts, nmaps
      INTEGER, DIMENSION(:), POINTER                     :: npts
      INTEGER, POINTER                                   :: ptr_ghatmap
      REAL(KIND=dp)                                      :: scale
      REAL(KIND=dp), POINTER                             :: ptr_pwout
      INTERFACE
         SUBROUTINE pw_gpu_sfffc_c(zin, dout, ghatmap, npts, ngpts, nmaps, scale) BIND(C, name="pw_gpu_sfffc")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: zin
            TYPE(C_PTR), VALUE                           :: dout
            TYPE(C_PTR), INTENT(IN), VALUE               :: ghatmap
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
            INTEGER(KIND=C_INT), INTENT(IN), VALUE       :: ngpts, nmaps
            REAL(KIND=C_DOUBLE), INTENT(IN), VALUE       :: scale
         END SUBROUTINE pw_gpu_sfffc_c
      END INTERFACE

      CALL timeset(routineN, handle)

      scale = 1.0_dp

      ngpts = SIZE(pw_grid_in%gsq)
      l1 = LBOUND(zout, 1)
      l2 = LBOUND(zout, 2)
      l3 = LBOUND(zout, 3)
      npts => pw_grid_in%npts

      ! pointers to data arrays
      ptr_pwin => zin(1)
      ptr_pwout => zout(l1, l2, l3)

      ! pointer to map array
      nmaps = SIZE(pw_grid_in%g_hatmap, 2)
      ptr_ghatmap => pw_grid_in%g_hatmap(1, 1)

      ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_sfffc_c(c_loc(ptr_pwin), c_loc(ptr_pwout), c_loc(ptr_ghatmap), npts, ngpts, nmaps, scale)
#else
      CPABORT("Compiled without pw offloading")
#endif

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_sfffc

! **************************************************************************************************
!> \brief perform a parallel real_to_complex copy followed by a 2D-FFT on the gpu
!> \param zin ...
!> \param pwbuf ...
!> \param pw_grid ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_cff(zin, pwbuf, pw_grid)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         TARGET                                          :: zin
      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(INOUT), TARGET                           :: pwbuf
      TYPE(pw_grid_type), INTENT(IN), TARGET             :: pw_grid

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_cff'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwout
      INTEGER                                            :: handle, l1, l2, l3
      INTEGER, DIMENSION(:), POINTER                     :: npts
      REAL(KIND=dp), POINTER                             :: ptr_pwin
      INTERFACE
         SUBROUTINE pw_gpu_cff_c(din, zout, npts) BIND(C, name="pw_gpu_cff")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: din
            TYPE(C_PTR), VALUE                           :: zout
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
         END SUBROUTINE pw_gpu_cff_c
      END INTERFACE

      CALL timeset(routineN, handle)

      ! dimensions
      npts => pw_grid%npts_local
      l1 = LBOUND(zin, 1)
      l2 = LBOUND(zin, 2)
      l3 = LBOUND(zin, 3)

      ! pointers to data arrays
      ptr_pwin => zin(l1, l2, l3)
      ptr_pwout => pwbuf(1, 1, 1)

      ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_cff_c(c_loc(ptr_pwin), c_loc(ptr_pwout), npts)
#else
      CPABORT("Compiled without pw offloading")
#endif

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_cff

! **************************************************************************************************
!> \brief perform a parallel 2D-FFT followed by a complex_to_real copy on the gpu
!> \param pwbuf ...
!> \param zout ...
!> \param pw_grid ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_ffc(pwbuf, zout, pw_grid)
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         TARGET                                          :: pwbuf
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         TARGET                                          :: zout
      TYPE(pw_grid_type), INTENT(IN), TARGET             :: pw_grid

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_ffc'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwin
      INTEGER                                            :: handle, l1, l2, l3
      INTEGER, DIMENSION(:), POINTER                     :: npts
      REAL(KIND=dp), POINTER                             :: ptr_pwout
      INTERFACE
         SUBROUTINE pw_gpu_ffc_c(zin, dout, npts) BIND(C, name="pw_gpu_ffc")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: zin
            TYPE(C_PTR), VALUE                           :: dout
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
         END SUBROUTINE pw_gpu_ffc_c
      END INTERFACE

      CALL timeset(routineN, handle)

      ! dimensions
      npts => pw_grid%npts_local
      l1 = LBOUND(zout, 1)
      l2 = LBOUND(zout, 2)
      l3 = LBOUND(zout, 3)

      ! pointers to data arrays
      ptr_pwin => pwbuf(1, 1, 1)
      ptr_pwout => zout(l1, l2, l3)

      ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_ffc_c(c_loc(ptr_pwin), c_loc(ptr_pwout), npts)
#else
      CPABORT("Compiled without pw offloading")
#endif

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_ffc

! **************************************************************************************************
!> \brief perform a parallel real_to_complex copy followed by a 1D-FFT on the gpu
!> \param zin ...
!> \param pwbuf ...
!> \param pw_grid ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_cf(zin, pwbuf, pw_grid)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN), &
         TARGET                                          :: zin
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         TARGET                                          :: pwbuf
      TYPE(pw_grid_type), INTENT(IN), TARGET             :: pw_grid

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_cf'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwout
      INTEGER                                            :: handle, l1, l2, l3
      INTEGER, DIMENSION(:), POINTER                     :: npts
      REAL(KIND=dp), POINTER                             :: ptr_pwin
      INTERFACE
         SUBROUTINE pw_gpu_cf_c(din, zout, npts) BIND(C, name="pw_gpu_cf")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: din
            TYPE(C_PTR), VALUE                           :: zout
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
         END SUBROUTINE pw_gpu_cf_c
      END INTERFACE

      CALL timeset(routineN, handle)

      ! dimensions
      npts => pw_grid%npts_local
      l1 = LBOUND(zin, 1)
      l2 = LBOUND(zin, 2)
      l3 = LBOUND(zin, 3)

      ! pointers to data arrays
      ptr_pwin => zin(l1, l2, l3)
      ptr_pwout => pwbuf(1, 1)

      ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_cf_c(c_loc(ptr_pwin), c_loc(ptr_pwout), npts)
#else
      CPABORT("Compiled without pw offloading")
#endif
      CALL timestop(handle)
   END SUBROUTINE pw_gpu_cf

! **************************************************************************************************
!> \brief perform a parallel 1D-FFT followed by a complex_to_real copy on the gpu
!> \param pwbuf ...
!> \param zout ...
!> \param pw_grid ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_fc(pwbuf, zout, pw_grid)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         TARGET                                          :: pwbuf
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         TARGET                                          :: zout
      TYPE(pw_grid_type), INTENT(IN), TARGET             :: pw_grid

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_fc'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwin
      INTEGER                                            :: handle, l1, l2, l3
      INTEGER, DIMENSION(:), POINTER                     :: npts
      REAL(KIND=dp), POINTER                             :: ptr_pwout
      INTERFACE
         SUBROUTINE pw_gpu_fc_c(zin, dout, npts) BIND(C, name="pw_gpu_fc")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: zin
            TYPE(C_PTR), VALUE                           :: dout
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
         END SUBROUTINE pw_gpu_fc_c
      END INTERFACE

      CALL timeset(routineN, handle)

      npts => pw_grid%npts_local
      l1 = LBOUND(zout, 1)
      l2 = LBOUND(zout, 2)
      l3 = LBOUND(zout, 3)

      ! pointers to data arrays
      ptr_pwin => pwbuf(1, 1)
      ptr_pwout => zout(l1, l2, l3)

      ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
      CALL pw_gpu_fc_c(c_loc(ptr_pwin), c_loc(ptr_pwout), npts)
#else
      CPABORT("Compiled without pw offloading")
#endif

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_fc

! **************************************************************************************************
!> \brief perform a parallel 1D-FFT on the gpu
!> \param pwbuf1 ...
!> \param pwbuf2 ...
!> \param dir ...
!> \param n ...
!> \param m ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_f(pwbuf1, pwbuf2, dir, n, m)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         TARGET                                          :: pwbuf1
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         TARGET                                          :: pwbuf2
      INTEGER, INTENT(IN)                                :: dir, n, m

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_f'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwin, ptr_pwout
      INTEGER                                            :: handle
      INTERFACE
         SUBROUTINE pw_gpu_f_c(zin, zout, dir, n, m) BIND(C, name="pw_gpu_f")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: zin
            TYPE(C_PTR), VALUE                           :: zout
            INTEGER(KIND=C_INT), INTENT(IN), VALUE       :: dir, n, m
         END SUBROUTINE pw_gpu_f_c
      END INTERFACE

      CALL timeset(routineN, handle)

      IF (n*m /= 0) THEN
         ! pointers to data arrays
         ptr_pwin => pwbuf1(1, 1)
         ptr_pwout => pwbuf2(1, 1)

         ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
         CALL pw_gpu_f_c(c_loc(ptr_pwin), c_loc(ptr_pwout), dir, n, m)
#else
         MARK_USED(dir)
         CPABORT("Compiled without pw offloading")
#endif
      END IF

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_f
! **************************************************************************************************
!> \brief perform a parallel 1D-FFT followed by a gather on the gpu
!> \param pwbuf ...
!> \param zout ...
!> \param scale ...
!> \param pw_grid ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_fg(pwbuf, zout, scale, pw_grid)
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         TARGET                                          :: pwbuf
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN), TARGET :: zout
      REAL(KIND=dp), INTENT(IN)                          :: scale
      TYPE(pw_grid_type), INTENT(IN), TARGET             :: pw_grid

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_fg'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwin, ptr_pwout
      INTEGER                                            :: handle, mg, mmax, ngpts
      INTEGER, DIMENSION(:), POINTER                     :: npts
      INTEGER, POINTER                                   :: ptr_ghatmap
      INTERFACE
         SUBROUTINE pw_gpu_fg_c(zin, zout, ghatmap, npts, mmax, ngpts, scale) BIND(C, name="pw_gpu_fg")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: zin
            TYPE(C_PTR), VALUE                           :: zout
            TYPE(C_PTR), INTENT(IN), VALUE               :: ghatmap
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
            INTEGER(KIND=C_INT), INTENT(IN), VALUE       :: mmax, ngpts
            REAL(KIND=C_DOUBLE), INTENT(IN), VALUE       :: scale

         END SUBROUTINE pw_gpu_fg_c
      END INTERFACE

      CALL timeset(routineN, handle)

      ngpts = SIZE(pw_grid%gsq)
      npts => pw_grid%npts

      IF ((npts(1) /= 0) .AND. (ngpts /= 0)) THEN
         mg = SIZE(pw_grid%grays, 2)
         mmax = MAX(mg, 1)

         ! pointers to data arrays
         ptr_pwin => pwbuf(1, 1)
         ptr_pwout => zout(1)

         ! pointer to map array
         ptr_ghatmap => pw_grid%g_hatmap(1, 1)

         ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
         CALL pw_gpu_fg_c(c_loc(ptr_pwin), c_loc(ptr_pwout), c_loc(ptr_ghatmap), npts, mmax, ngpts, scale)
#else
         MARK_USED(scale)
         CPABORT("Compiled without pw offloading")
#endif
      END IF

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_fg

! **************************************************************************************************
!> \brief perform a parallel scatter followed by a 1D-FFT on the gpu
!> \param zin ...
!> \param pwbuf ...
!> \param scale ...
!> \param pw_grid ...
!> \author Andreas Gloess
! **************************************************************************************************
   SUBROUTINE pw_gpu_sf(zin, pwbuf, scale, pw_grid)
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN), TARGET :: zin
      COMPLEX(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         TARGET                                          :: pwbuf
      REAL(KIND=dp), INTENT(IN)                          :: scale
      TYPE(pw_grid_type), INTENT(IN), TARGET             :: pw_grid

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_gpu_sf'

      COMPLEX(KIND=dp), POINTER                          :: ptr_pwin, ptr_pwout
      INTEGER                                            :: handle, mg, mmax, ngpts, nmaps
      INTEGER, DIMENSION(:), POINTER                     :: npts
      INTEGER, POINTER                                   :: ptr_ghatmap
      INTERFACE
         SUBROUTINE pw_gpu_sf_c(zin, zout, ghatmap, npts, mmax, ngpts, nmaps, scale) BIND(C, name="pw_gpu_sf")
            IMPORT
            TYPE(C_PTR), INTENT(IN), VALUE               :: zin
            TYPE(C_PTR), VALUE                           :: zout
            TYPE(C_PTR), INTENT(IN), VALUE               :: ghatmap
            INTEGER(KIND=C_INT), DIMENSION(*), INTENT(IN):: npts
            INTEGER(KIND=C_INT), INTENT(IN), VALUE       :: mmax, ngpts, nmaps
            REAL(KIND=C_DOUBLE), INTENT(IN), VALUE       :: scale

         END SUBROUTINE pw_gpu_sf_c
      END INTERFACE

      CALL timeset(routineN, handle)

      ngpts = SIZE(pw_grid%gsq)
      npts => pw_grid%npts

      IF ((npts(1) /= 0) .AND. (ngpts /= 0)) THEN
         mg = SIZE(pw_grid%grays, 2)
         mmax = MAX(mg, 1)

         ! pointers to data arrays
         ptr_pwin => zin(1)
         ptr_pwout => pwbuf(1, 1)

         ! pointer to map array
         nmaps = SIZE(pw_grid%g_hatmap, 2)
         ptr_ghatmap => pw_grid%g_hatmap(1, 1)

         ! invoke the combined transformation
#if defined(__OFFLOAD) && !defined(__NO_OFFLOAD_PW)
         CALL pw_gpu_sf_c(c_loc(ptr_pwin), c_loc(ptr_pwout), c_loc(ptr_ghatmap), npts, mmax, ngpts, nmaps, scale)
#else
         MARK_USED(scale)
         CPABORT("Compiled without pw offloading")
#endif
      END IF

      CALL timestop(handle)
   END SUBROUTINE pw_gpu_sf

END MODULE pw_gpu

