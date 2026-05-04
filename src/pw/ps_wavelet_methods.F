!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Definition and initialisation of the ps_wavelet data type.
!> \history 01.2014 Renamed from ps_wavelet_types to disentangle dependencies (Ole Schuett)
!> \author Florian Schiffmann (09.2007,fschiff)
! **************************************************************************************************
MODULE ps_wavelet_methods

   USE bibliography,                    ONLY: Genovese2006,&
                                              Genovese2007,&
                                              cite_reference
   USE kinds,                           ONLY: dp
   USE ps_wavelet_kernel,               ONLY: createKernel
   USE ps_wavelet_types,                ONLY: WAVELET0D,&
                                              WAVELET2D,&
                                              ps_wavelet_release,&
                                              ps_wavelet_type
   USE ps_wavelet_util,                 ONLY: F_FFT_dimensions,&
                                              PSolver,&
                                              P_FFT_dimensions,&
                                              S_FFT_dimensions
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_poisson_types,                ONLY: pw_poisson_parameter_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE util,                            ONLY: get_limit
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ps_wavelet_methods'

! *** Public data types ***

   PUBLIC :: ps_wavelet_create, &
             cp2k_distribution_to_z_slices, &
             z_slices_to_cp2k_distribution, &
             ps_wavelet_solve

CONTAINS

! **************************************************************************************************
!> \brief creates the ps_wavelet_type which is needed for the link to
!>      the Poisson Solver of Luigi Genovese
!> \param poisson_params ...
!> \param wavelet wavelet to create
!> \param pw_grid the grid that is used to create the wavelet kernel
!> \author Flroian Schiffmann
! **************************************************************************************************
   SUBROUTINE ps_wavelet_create(poisson_params, wavelet, pw_grid)
      TYPE(pw_poisson_parameter_type), INTENT(IN)        :: poisson_params
      TYPE(ps_wavelet_type), POINTER                     :: wavelet
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ps_wavelet_create'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: hx, hy, hz

      CALL timeset(routineN, handle)

      CALL cite_reference(Genovese2006)
      CALL cite_reference(Genovese2007)

      IF (ASSOCIATED(wavelet)) THEN
         CALL ps_wavelet_release(wavelet)
         NULLIFY (wavelet)
      END IF

      ALLOCATE (wavelet)

      NULLIFY (wavelet%karray, wavelet%rho_z_sliced)

      wavelet%geocode = poisson_params%wavelet_geocode
      wavelet%method = poisson_params%wavelet_method
      wavelet%special_dimension = poisson_params%wavelet_special_dimension
      wavelet%itype_scf = poisson_params%wavelet_scf_type
      wavelet%datacode = "D"

      CALL set_wavelet_axis(wavelet)
      hx = pw_grid%dr(wavelet%axis(1))
      hy = pw_grid%dr(wavelet%axis(2))
      hz = pw_grid%dr(wavelet%axis(3))

      IF (poisson_params%wavelet_method == WAVELET0D) THEN
         IF (hx /= hy) &
            CPABORT("Poisson solver for non cubic cells not yet implemented")
         IF (hz /= hy) &
            CPABORT("Poisson solver for non cubic cells not yet implemented")
      END IF

      CALL RS_z_slice_distribution(wavelet, pw_grid)

      CALL timestop(handle)
   END SUBROUTINE ps_wavelet_create

! **************************************************************************************************
!> \brief ...
!> \param wavelet ...
!> \param pw_grid ...
! **************************************************************************************************
   SUBROUTINE RS_z_slice_distribution(wavelet, pw_grid)

      TYPE(ps_wavelet_type), POINTER                     :: wavelet
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      CHARACTER(len=*), PARAMETER :: routineN = 'RS_z_slice_distribution'

      CHARACTER(LEN=1)                                   :: geocode
      INTEGER                                            :: handle, iproc, m1, m2, m3, md1, md2, &
                                                            md3, n1, n2, n3, nd1, nd2, nd3, nproc, &
                                                            nx, ny, nz, z_dim
      REAL(KIND=dp)                                      :: hx, hy, hz

      CALL timeset(routineN, handle)
      nproc = PRODUCT(pw_grid%para%group%num_pe_cart)
      iproc = pw_grid%para%group%mepos
      geocode = wavelet%geocode
      CALL get_wavelet_grid(wavelet, pw_grid, nx, ny, nz, hx, hy, hz)

      !calculate Dimensions for the z-distributed density and for the kernel

      IF (geocode == 'P') THEN
         CALL P_FFT_dimensions(nx, ny, nz, m1, m2, m3, n1, n2, n3, md1, md2, md3, nd1, nd2, nd3, nproc)
      ELSE IF (geocode == 'S') THEN
         CALL S_FFT_dimensions(nx, ny, nz, m1, m2, m3, n1, n2, n3, md1, md2, md3, nd1, nd2, nd3, nproc)
      ELSE IF (geocode == 'F') THEN
         CALL F_FFT_dimensions(nx, ny, nz, m1, m2, m3, n1, n2, n3, md1, md2, md3, nd1, nd2, nd3, nproc)
      END IF

      wavelet%PS_grid(1) = md1
      wavelet%PS_grid(2) = md3
      wavelet%PS_grid(3) = md2
      z_dim = md2/nproc
      !!!!!!!!!      indices y and z are interchanged    !!!!!!!
      ALLOCATE (wavelet%rho_z_sliced(md1, md3, z_dim))

      CALL createKernel(geocode, nx, ny, nz, hx, hy, hz, wavelet%itype_scf, iproc, nproc, wavelet%karray, &
                        pw_grid%para%group)

      CALL timestop(handle)
   END SUBROUTINE RS_z_slice_distribution

! **************************************************************************************************
!> \brief ...
!> \param density ...
!> \param wavelet ...
!> \param pw_grid ...
! **************************************************************************************************
   SUBROUTINE cp2k_distribution_to_z_slices(density, wavelet, pw_grid)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: density
      TYPE(ps_wavelet_type), POINTER                     :: wavelet
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      CHARACTER(len=*), PARAMETER :: routineN = 'cp2k_distribution_to_z_slices'

      INTEGER                                            :: dest, handle, i, ii, iproc, j, k, l, &
                                                            local_z_dim, loz, m, m2, md2, nproc, &
                                                            should_warn
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rcount, rdispl, scount, sdispl, tmp
      INTEGER, DIMENSION(2)                              :: cart_pos, lox, loy
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: max_val_low, max_val_up
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rbuf, sbuf

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(wavelet))

      IF (.NOT. wavelet_axis_is_identity(wavelet)) THEN
         CALL warn_density_edges(density, wavelet, pw_grid)
         CALL cp2k_distribution_to_z_slices_permuted(density, wavelet, pw_grid)
         CALL timestop(handle)
         RETURN
      END IF

      nproc = PRODUCT(pw_grid%para%group%num_pe_cart)
      iproc = pw_grid%para%group%mepos
      md2 = wavelet%PS_grid(3)
      m2 = pw_grid%npts(3)
      lb(:) = pw_grid%bounds_local(1, :)
      ub(:) = pw_grid%bounds_local(2, :)
      local_z_dim = MAX((md2/nproc), 1)

      ALLOCATE (sbuf(PRODUCT(pw_grid%npts_local)))
      ALLOCATE (rbuf(PRODUCT(wavelet%PS_grid)/nproc))
      ALLOCATE (scount(nproc), sdispl(nproc), rcount(nproc), rdispl(nproc), tmp(nproc))

      rbuf = 0.0_dp
      ii = 1
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               sbuf(ii) = density%array(i, j, k)
               ii = ii + 1
            END DO
         END DO
      END DO

      should_warn = 0
      IF (wavelet%geocode == 'S' .OR. wavelet%geocode == 'F') THEN
         max_val_low = 0._dp
         max_val_up = 0._dp
         IF (lb(2) == pw_grid%bounds(1, 2)) max_val_low = MAXVAL(ABS(density%array(:, lb(2), :)))
         IF (ub(2) == pw_grid%bounds(2, 2)) max_val_up = MAXVAL(ABS(density%array(:, ub(2), :)))
         IF (max_val_low >= 0.0001_dp) should_warn = 1
         IF (max_val_up >= 0.0001_dp) should_warn = 1
         IF (wavelet%geocode == 'F') THEN
            max_val_low = 0._dp
            max_val_up = 0._dp
            IF (lb(1) == pw_grid%bounds(1, 1)) max_val_low = MAXVAL(ABS(density%array(lb(1), :, :)))
            IF (ub(1) == pw_grid%bounds(2, 1)) max_val_up = MAXVAL(ABS(density%array(ub(1), :, :)))
            IF (max_val_low >= 0.0001_dp) should_warn = 1
            IF (max_val_up >= 0.0001_dp) should_warn = 1
            max_val_low = 0._dp
            max_val_up = 0._dp
            IF (lb(3) == pw_grid%bounds(1, 3)) max_val_low = MAXVAL(ABS(density%array(:, :, lb(3))))
            IF (ub(3) == pw_grid%bounds(2, 3)) max_val_up = MAXVAL(ABS(density%array(:, :, ub(3))))
            IF (max_val_low >= 0.0001_dp) should_warn = 1
            IF (max_val_up >= 0.0001_dp) should_warn = 1
         END IF
      END IF

      CALL pw_grid%para%group%max(should_warn)
      IF (should_warn > 0 .AND. iproc == 0) THEN
         CPWARN("Density non-zero on the edges of the unit cell: wrong results in WAVELET solver")
      END IF
      DO i = 0, pw_grid%para%group%num_pe_cart(1) - 1
         DO j = 0, pw_grid%para%group%num_pe_cart(2) - 1
            cart_pos = [i, j]
            CALL pw_grid%para%group%rank_cart(cart_pos, dest)
            IF ((ub(1) >= lb(1)) .AND. (ub(2) >= lb(2))) THEN
               IF (dest*local_z_dim <= m2) THEN
                  IF ((dest + 1)*local_z_dim <= m2) THEN
                     scount(dest + 1) = ABS((ub(1) - lb(1) + 1)*(ub(2) - lb(2) + 1)*local_z_dim)
                  ELSE
                     scount(dest + 1) = ABS((ub(1) - lb(1) + 1)*(ub(2) - lb(2) + 1)*MOD(m2, local_z_dim))
                  END IF
               ELSE
                  scount(dest + 1) = 0
               END IF
            ELSE
               scount(dest + 1) = 0
            END IF
            lox = get_limit(pw_grid%npts(1), pw_grid%para%group%num_pe_cart(1), i)
            loy = get_limit(pw_grid%npts(2), pw_grid%para%group%num_pe_cart(2), j)
            IF ((lox(2) >= lox(1)) .AND. (loy(2) >= loy(1))) THEN
               IF (iproc*local_z_dim <= m2) THEN
                  IF ((iproc + 1)*local_z_dim <= m2) THEN
                     rcount(dest + 1) = ABS((lox(2) - lox(1) + 1)*(loy(2) - loy(1) + 1)*local_z_dim)
                  ELSE
                     rcount(dest + 1) = ABS((lox(2) - lox(1) + 1)*(loy(2) - loy(1) + 1)*MOD(m2, local_z_dim))
                  END IF
               ELSE
                  rcount(dest + 1) = 0
               END IF
            ELSE
               rcount(dest + 1) = 0
            END IF

         END DO
      END DO
      sdispl(1) = 0
      rdispl(1) = 0
      DO i = 2, nproc
         sdispl(i) = sdispl(i - 1) + scount(i - 1)
         rdispl(i) = rdispl(i - 1) + rcount(i - 1)
      END DO
      CALL pw_grid%para%group%alltoall(sbuf, scount, sdispl, rbuf, rcount, rdispl)
      !!!! and now, how to put the right cubes to the right position!!!!!!

      wavelet%rho_z_sliced = 0.0_dp

      DO i = 0, pw_grid%para%group%num_pe_cart(1) - 1
         DO j = 0, pw_grid%para%group%num_pe_cart(2) - 1
            cart_pos = [i, j]
            CALL pw_grid%para%group%rank_cart(cart_pos, dest)

            lox = get_limit(pw_grid%npts(1), pw_grid%para%group%num_pe_cart(1), i)
            loy = get_limit(pw_grid%npts(2), pw_grid%para%group%num_pe_cart(2), j)
            IF (iproc*local_z_dim <= m2) THEN
               IF ((iproc + 1)*local_z_dim <= m2) THEN
                  loz = local_z_dim
               ELSE
                  loz = MOD(m2, local_z_dim)
               END IF
               ii = 1
               DO k = 1, loz
                  DO l = loy(1), loy(2)
                     DO m = lox(1), lox(2)
                        wavelet%rho_z_sliced(m, l, k) = rbuf(ii + rdispl(dest + 1))
                        ii = ii + 1
                     END DO
                  END DO
               END DO
            END IF
         END DO
      END DO

      DEALLOCATE (sbuf, rbuf, scount, sdispl, rcount, rdispl, tmp)

      CALL timestop(handle)

   END SUBROUTINE cp2k_distribution_to_z_slices

! **************************************************************************************************
!> \brief ...
!> \param density ...
!> \param wavelet ...
!> \param pw_grid ...
! **************************************************************************************************
   SUBROUTINE z_slices_to_cp2k_distribution(density, wavelet, pw_grid)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: density
      TYPE(ps_wavelet_type), POINTER                     :: wavelet
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      INTEGER                                            :: dest, i, ii, iproc, j, k, l, &
                                                            local_z_dim, loz, m, m2, md2, nproc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rcount, rdispl, scount, sdispl, tmp
      INTEGER, DIMENSION(2)                              :: cart_pos, lox, loy, min_x, min_y
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rbuf, sbuf

      CPASSERT(ASSOCIATED(wavelet))

      IF (.NOT. wavelet_axis_is_identity(wavelet)) THEN
         CALL z_slices_to_cp2k_distribution_permuted(density, wavelet, pw_grid)
         RETURN
      END IF

      nproc = PRODUCT(pw_grid%para%group%num_pe_cart)
      iproc = pw_grid%para%group%mepos
      md2 = wavelet%PS_grid(3)
      m2 = pw_grid%npts(3)

      lb(:) = pw_grid%bounds_local(1, :)
      ub(:) = pw_grid%bounds_local(2, :)

      local_z_dim = MAX((md2/nproc), 1)

      ALLOCATE (rbuf(PRODUCT(pw_grid%npts_local)))
      ALLOCATE (sbuf(PRODUCT(wavelet%PS_grid)/nproc))
      ALLOCATE (scount(nproc), sdispl(nproc), rcount(nproc), rdispl(nproc), tmp(nproc))
      scount = 0
      rcount = 0
      rbuf = 0.0_dp
      ii = 1
      IF (iproc*local_z_dim <= m2) THEN
         IF ((iproc + 1)*local_z_dim <= m2) THEN
            loz = local_z_dim
         ELSE
            loz = MOD(m2, local_z_dim)
         END IF
      ELSE
         loz = 0
      END IF

      min_x = get_limit(pw_grid%npts(1), pw_grid%para%group%num_pe_cart(1), 0)
      min_y = get_limit(pw_grid%npts(2), pw_grid%para%group%num_pe_cart(2), 0)
      DO i = 0, pw_grid%para%group%num_pe_cart(1) - 1
         DO j = 0, pw_grid%para%group%num_pe_cart(2) - 1
            cart_pos = [i, j]
            CALL pw_grid%para%group%rank_cart(cart_pos, dest)
            IF ((ub(1) >= lb(1)) .AND. (ub(2) >= lb(2))) THEN
               IF (dest*local_z_dim <= m2) THEN
                  IF ((dest + 1)*local_z_dim <= m2) THEN
                     rcount(dest + 1) = ABS((ub(1) - lb(1) + 1)*(ub(2) - lb(2) + 1)*local_z_dim)
                  ELSE
                     rcount(dest + 1) = ABS((ub(1) - lb(1) + 1)*(ub(2) - lb(2) + 1)*MOD(m2, local_z_dim))
                  END IF
               ELSE
                  rcount(dest + 1) = 0
               END IF
            ELSE
               rcount(dest + 1) = 0
            END IF
            lox = get_limit(pw_grid%npts(1), pw_grid%para%group%num_pe_cart(1), i)
            loy = get_limit(pw_grid%npts(2), pw_grid%para%group%num_pe_cart(2), j)
            IF ((lox(2) >= lox(1)) .AND. (loy(2) >= loy(1))) THEN
               scount(dest + 1) = ABS((lox(2) - lox(1) + 1)*(loy(2) - loy(1) + 1)*loz)
               DO k = lox(1) - min_x(1) + 1, lox(2) - min_x(1) + 1
                  DO l = loy(1) - min_y(1) + 1, loy(2) - min_y(1) + 1
                     DO m = 1, loz
                        sbuf(ii) = wavelet%rho_z_sliced(k, l, m)
                        ii = ii + 1
                     END DO
                  END DO
               END DO
            ELSE
               scount(dest + 1) = 0
            END IF
         END DO
      END DO
      sdispl(1) = 0
      rdispl(1) = 0
      DO i = 2, nproc
         sdispl(i) = sdispl(i - 1) + scount(i - 1)
         rdispl(i) = rdispl(i - 1) + rcount(i - 1)
      END DO
      CALL pw_grid%para%group%alltoall(sbuf, scount, sdispl, rbuf, rcount, rdispl)

      !!!! and now, how to put the right cubes to the right position!!!!!!

      DO i = 0, pw_grid%para%group%num_pe_cart(1) - 1
         DO j = 0, pw_grid%para%group%num_pe_cart(2) - 1
            cart_pos = [i, j]
            CALL pw_grid%para%group%rank_cart(cart_pos, dest)
            IF (dest*local_z_dim <= m2) THEN
               IF ((dest + 1)*local_z_dim <= m2) THEN
                  loz = local_z_dim
               ELSE
                  loz = MOD(m2, local_z_dim)
               END IF
               ii = 1
               IF (lb(3) + (dest*local_z_dim) <= ub(3)) THEN
                  DO m = lb(1), ub(1)
                     DO l = lb(2), ub(2)
                        DO k = lb(3) + (dest*local_z_dim), lb(3) + (dest*local_z_dim) + loz - 1
                           density%array(m, l, k) = rbuf(ii + rdispl(dest + 1))
                           ii = ii + 1
                        END DO
                     END DO
                  END DO
               END IF
            END IF
         END DO
      END DO
      DEALLOCATE (sbuf, rbuf, scount, sdispl, rcount, rdispl, tmp)

   END SUBROUTINE z_slices_to_cp2k_distribution

! **************************************************************************************************
!> \brief Set the internal Wavelet solver axes for the requested boundary conditions.
!> \param wavelet ...
! **************************************************************************************************
   SUBROUTINE set_wavelet_axis(wavelet)

      TYPE(ps_wavelet_type), POINTER                     :: wavelet

      wavelet%axis = [1, 2, 3]

      IF (wavelet%method == WAVELET2D) THEN
         SELECT CASE (wavelet%special_dimension)
         CASE (1)
            wavelet%axis = [2, 1, 3]
         CASE (2)
            wavelet%axis = [1, 2, 3]
         CASE (3)
            wavelet%axis = [1, 3, 2]
         CASE DEFAULT
            CPABORT("Invalid isolated dimension for WAVELET 2D")
         END SELECT
      END IF

   END SUBROUTINE set_wavelet_axis

! **************************************************************************************************
!> \brief Return grid sizes and spacings in the internal Wavelet solver axes.
!> \param wavelet ...
!> \param pw_grid ...
!> \param nx ...
!> \param ny ...
!> \param nz ...
!> \param hx ...
!> \param hy ...
!> \param hz ...
! **************************************************************************************************
   SUBROUTINE get_wavelet_grid(wavelet, pw_grid, nx, ny, nz, hx, hy, hz)

      TYPE(ps_wavelet_type), POINTER                     :: wavelet
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      INTEGER, INTENT(OUT)                               :: nx, ny, nz
      REAL(KIND=dp), INTENT(OUT)                         :: hx, hy, hz

      nx = pw_grid%npts(wavelet%axis(1))
      ny = pw_grid%npts(wavelet%axis(2))
      nz = pw_grid%npts(wavelet%axis(3))
      hx = pw_grid%dr(wavelet%axis(1))
      hy = pw_grid%dr(wavelet%axis(2))
      hz = pw_grid%dr(wavelet%axis(3))

   END SUBROUTINE get_wavelet_grid

! **************************************************************************************************
!> \brief Return whether CP2K and internal Wavelet solver axes are identical.
!> \param wavelet ...
!> \return ...
! **************************************************************************************************
   FUNCTION wavelet_axis_is_identity(wavelet) RESULT(is_identity)

      TYPE(ps_wavelet_type), POINTER                     :: wavelet
      LOGICAL                                            :: is_identity

      is_identity = ALL(wavelet%axis == [1, 2, 3])

   END FUNCTION wavelet_axis_is_identity

! **************************************************************************************************
!> \brief Warn if the density is non-zero at isolated Wavelet cell edges.
!> \param density ...
!> \param wavelet ...
!> \param pw_grid ...
! **************************************************************************************************
   SUBROUTINE warn_density_edges(density, wavelet, pw_grid)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: density
      TYPE(ps_wavelet_type), POINTER                     :: wavelet
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      INTEGER                                            :: idir, iproc, should_warn

      should_warn = 0
      iproc = pw_grid%para%group%mepos

      IF (wavelet%geocode == 'S') THEN
         CALL update_edge_warning(density, pw_grid, wavelet%special_dimension, should_warn)
      ELSE IF (wavelet%geocode == 'F') THEN
         DO idir = 1, 3
            CALL update_edge_warning(density, pw_grid, idir, should_warn)
         END DO
      END IF

      CALL pw_grid%para%group%max(should_warn)
      IF (should_warn > 0 .AND. iproc == 0) THEN
         CPWARN("Density non-zero on the edges of the unit cell: wrong results in WAVELET solver")
      END IF

   END SUBROUTINE warn_density_edges

! **************************************************************************************************
!> \brief Update the density edge warning for one real-space direction.
!> \param density ...
!> \param pw_grid ...
!> \param direction ...
!> \param should_warn ...
! **************************************************************************************************
   SUBROUTINE update_edge_warning(density, pw_grid, direction, should_warn)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: density
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      INTEGER, INTENT(IN)                                :: direction
      INTEGER, INTENT(INOUT)                             :: should_warn

      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: max_val_low, max_val_up

      lb(:) = pw_grid%bounds_local(1, :)
      ub(:) = pw_grid%bounds_local(2, :)
      IF (.NOT. ALL(ub >= lb)) RETURN

      max_val_low = 0._dp
      max_val_up = 0._dp
      SELECT CASE (direction)
      CASE (1)
         IF (lb(1) == pw_grid%bounds(1, 1)) max_val_low = MAXVAL(ABS(density%array(lb(1), :, :)))
         IF (ub(1) == pw_grid%bounds(2, 1)) max_val_up = MAXVAL(ABS(density%array(ub(1), :, :)))
      CASE (2)
         IF (lb(2) == pw_grid%bounds(1, 2)) max_val_low = MAXVAL(ABS(density%array(:, lb(2), :)))
         IF (ub(2) == pw_grid%bounds(2, 2)) max_val_up = MAXVAL(ABS(density%array(:, ub(2), :)))
      CASE (3)
         IF (lb(3) == pw_grid%bounds(1, 3)) max_val_low = MAXVAL(ABS(density%array(:, :, lb(3))))
         IF (ub(3) == pw_grid%bounds(2, 3)) max_val_up = MAXVAL(ABS(density%array(:, :, ub(3))))
      CASE DEFAULT
         CPABORT("Invalid WAVELET isolated dimension")
      END SELECT

      IF (max_val_low >= 0.0001_dp) should_warn = 1
      IF (max_val_up >= 0.0001_dp) should_warn = 1

   END SUBROUTINE update_edge_warning

! **************************************************************************************************
!> \brief Convert counts into zero-based displacements for all-to-all communication.
!> \param counts ...
!> \param displacements ...
! **************************************************************************************************
   SUBROUTINE set_displacements(counts, displacements)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: counts
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: displacements

      INTEGER                                            :: i

      displacements(1) = 0
      DO i = 2, SIZE(counts)
         displacements(i) = displacements(i - 1) + counts(i - 1)
      END DO

   END SUBROUTINE set_displacements

! **************************************************************************************************
!> \brief Return the distributed grid owner of a compact one-based grid index.
!> \param index ...
!> \param npts ...
!> \param nparts ...
!> \return ...
! **************************************************************************************************
   FUNCTION grid_owner(index, npts, nparts) RESULT(owner)

      INTEGER, INTENT(IN)                                :: index, npts, nparts
      INTEGER                                            :: owner

      INTEGER                                            :: ipart
      INTEGER, DIMENSION(2)                              :: limits

      DO ipart = 0, nparts - 1
         limits = get_limit(npts, nparts, ipart)
         IF (index >= limits(1) .AND. index <= limits(2)) THEN
            owner = ipart
            RETURN
         END IF
      END DO
      CPABORT("Grid index outside distributed bounds")

   END FUNCTION grid_owner

! **************************************************************************************************
!> \brief Return the Wavelet z-slice owner of a compact one-based grid index.
!> \param index ...
!> \param local_z_dim ...
!> \param nproc ...
!> \return ...
! **************************************************************************************************
   FUNCTION z_slice_owner(index, local_z_dim, nproc) RESULT(owner)

      INTEGER, INTENT(IN)                                :: index, local_z_dim, nproc
      INTEGER                                            :: owner

      owner = MIN((index - 1)/local_z_dim, nproc - 1)

   END FUNCTION z_slice_owner

! **************************************************************************************************
!> \brief Return the CP2K real-space rank owning a real-space x/y grid point.
!> \param ix ...
!> \param iy ...
!> \param pw_grid ...
!> \return ...
! **************************************************************************************************
   FUNCTION cp2k_rank_owner(ix, iy, pw_grid) RESULT(rank)

      INTEGER, INTENT(IN)                                :: ix, iy
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      INTEGER                                            :: rank

      INTEGER, DIMENSION(2)                              :: cart_pos

      cart_pos(1) = grid_owner(ix - pw_grid%bounds(1, 1) + 1, pw_grid%npts(1), &
                               pw_grid%para%group%num_pe_cart(1))
      cart_pos(2) = grid_owner(iy - pw_grid%bounds(1, 2) + 1, pw_grid%npts(2), &
                               pw_grid%para%group%num_pe_cart(2))
      CALL pw_grid%para%group%rank_cart(cart_pos, rank)

   END FUNCTION cp2k_rank_owner

! **************************************************************************************************
!> \brief Transfer a CP2K real-space grid into a permuted Wavelet z-sliced layout.
!> \param density ...
!> \param wavelet ...
!> \param pw_grid ...
! **************************************************************************************************
   SUBROUTINE cp2k_distribution_to_z_slices_permuted(density, wavelet, pw_grid)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: density
      TYPE(ps_wavelet_type), POINTER                     :: wavelet
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      INTEGER                                            :: dest, i, idir, ii, j, k, local_z_dim, &
                                                            md2, nproc, nrecv, nsend
      INTEGER, ALLOCATABLE, DIMENSION(:) :: coord_rbuf, coord_rcount, coord_rdispl, coord_sbuf, &
         coord_scount, coord_sdispl, rcount, rdispl, scount, sdispl, send_pos
      INTEGER, DIMENSION(3)                              :: lb, q, u, ub
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: rbuf, sbuf

      nproc = PRODUCT(pw_grid%para%group%num_pe_cart)
      md2 = wavelet%PS_grid(3)
      local_z_dim = MAX(md2/nproc, 1)
      lb(:) = pw_grid%bounds_local(1, :)
      ub(:) = pw_grid%bounds_local(2, :)

      ALLOCATE (scount(nproc), sdispl(nproc), rcount(nproc), rdispl(nproc), send_pos(nproc))
      scount = 0
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               u = [i, j, k]
               DO idir = 1, 3
                  q(idir) = u(wavelet%axis(idir)) - pw_grid%bounds(1, wavelet%axis(idir)) + 1
               END DO
               dest = z_slice_owner(q(3), local_z_dim, nproc)
               scount(dest + 1) = scount(dest + 1) + 1
            END DO
         END DO
      END DO

      CALL pw_grid%para%group%alltoall(scount, rcount, 1)
      CALL set_displacements(scount, sdispl)
      CALL set_displacements(rcount, rdispl)
      nsend = SUM(scount)
      nrecv = SUM(rcount)

      ALLOCATE (sbuf(MAX(nsend, 1)), rbuf(MAX(nrecv, 1)))
      ALLOCATE (coord_sbuf(MAX(3*nsend, 1)), coord_rbuf(MAX(3*nrecv, 1)))
      ALLOCATE (coord_scount(nproc), coord_sdispl(nproc), coord_rcount(nproc), coord_rdispl(nproc))
      coord_scount(:) = 3*scount(:)
      coord_rcount(:) = 3*rcount(:)
      coord_sdispl(:) = 3*sdispl(:)
      coord_rdispl(:) = 3*rdispl(:)
      send_pos(:) = sdispl(:) + 1

      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               u = [i, j, k]
               DO idir = 1, 3
                  q(idir) = u(wavelet%axis(idir)) - pw_grid%bounds(1, wavelet%axis(idir)) + 1
               END DO
               dest = z_slice_owner(q(3), local_z_dim, nproc)
               ii = send_pos(dest + 1)
               sbuf(ii) = density%array(i, j, k)
               coord_sbuf(3*ii - 2) = q(1)
               coord_sbuf(3*ii - 1) = q(2)
               coord_sbuf(3*ii) = q(3)
               send_pos(dest + 1) = ii + 1
            END DO
         END DO
      END DO

      CALL pw_grid%para%group%alltoall(sbuf, scount, sdispl, rbuf, rcount, rdispl)
      CALL pw_grid%para%group%alltoall(coord_sbuf, coord_scount, coord_sdispl, &
                                       coord_rbuf, coord_rcount, coord_rdispl)

      wavelet%rho_z_sliced = 0.0_dp
      DO ii = 1, nrecv
         q(1) = coord_rbuf(3*ii - 2)
         q(2) = coord_rbuf(3*ii - 1)
         q(3) = coord_rbuf(3*ii)
         wavelet%rho_z_sliced(q(1), q(2), q(3) - pw_grid%para%group%mepos*local_z_dim) = rbuf(ii)
      END DO

      DEALLOCATE (sbuf, rbuf, coord_sbuf, coord_rbuf, coord_scount, coord_sdispl, coord_rcount, &
                  coord_rdispl, scount, sdispl, rcount, rdispl, send_pos)

   END SUBROUTINE cp2k_distribution_to_z_slices_permuted

! **************************************************************************************************
!> \brief Transfer a permuted Wavelet z-sliced layout back to a CP2K real-space grid.
!> \param density ...
!> \param wavelet ...
!> \param pw_grid ...
! **************************************************************************************************
   SUBROUTINE z_slices_to_cp2k_distribution_permuted(density, wavelet, pw_grid)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: density
      TYPE(ps_wavelet_type), POINTER                     :: wavelet
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      INTEGER                                            :: dest, i, idir, ii, j, k, local_z, &
                                                            local_z_dim, md2, n1, n2, n3, nproc, &
                                                            nrecv, nsend, z_end, z_start
      INTEGER, ALLOCATABLE, DIMENSION(:) :: coord_rbuf, coord_rcount, coord_rdispl, coord_sbuf, &
         coord_scount, coord_sdispl, rcount, rdispl, scount, sdispl, send_pos
      INTEGER, DIMENSION(3)                              :: lb, q, u, ub
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: rbuf, sbuf

      nproc = PRODUCT(pw_grid%para%group%num_pe_cart)
      md2 = wavelet%PS_grid(3)
      local_z_dim = MAX(md2/nproc, 1)
      n1 = pw_grid%npts(wavelet%axis(1))
      n2 = pw_grid%npts(wavelet%axis(2))
      n3 = pw_grid%npts(wavelet%axis(3))
      z_start = pw_grid%para%group%mepos*local_z_dim + 1
      z_end = MIN((pw_grid%para%group%mepos + 1)*local_z_dim, n3)
      lb(:) = pw_grid%bounds_local(1, :)
      ub(:) = pw_grid%bounds_local(2, :)

      ALLOCATE (scount(nproc), sdispl(nproc), rcount(nproc), rdispl(nproc), send_pos(nproc))
      scount = 0
      DO k = z_start, z_end
         DO j = 1, n2
            DO i = 1, n1
               q = [i, j, k]
               DO idir = 1, 3
                  u(wavelet%axis(idir)) = q(idir) + pw_grid%bounds(1, wavelet%axis(idir)) - 1
               END DO
               dest = cp2k_rank_owner(u(1), u(2), pw_grid)
               scount(dest + 1) = scount(dest + 1) + 1
            END DO
         END DO
      END DO

      CALL pw_grid%para%group%alltoall(scount, rcount, 1)
      CALL set_displacements(scount, sdispl)
      CALL set_displacements(rcount, rdispl)
      nsend = SUM(scount)
      nrecv = SUM(rcount)

      ALLOCATE (sbuf(MAX(nsend, 1)), rbuf(MAX(nrecv, 1)))
      ALLOCATE (coord_sbuf(MAX(3*nsend, 1)), coord_rbuf(MAX(3*nrecv, 1)))
      ALLOCATE (coord_scount(nproc), coord_sdispl(nproc), coord_rcount(nproc), coord_rdispl(nproc))
      coord_scount(:) = 3*scount(:)
      coord_rcount(:) = 3*rcount(:)
      coord_sdispl(:) = 3*sdispl(:)
      coord_rdispl(:) = 3*rdispl(:)
      send_pos(:) = sdispl(:) + 1

      DO k = z_start, z_end
         DO j = 1, n2
            DO i = 1, n1
               q = [i, j, k]
               DO idir = 1, 3
                  u(wavelet%axis(idir)) = q(idir) + pw_grid%bounds(1, wavelet%axis(idir)) - 1
               END DO
               dest = cp2k_rank_owner(u(1), u(2), pw_grid)
               ii = send_pos(dest + 1)
               local_z = k - pw_grid%para%group%mepos*local_z_dim
               sbuf(ii) = wavelet%rho_z_sliced(i, j, local_z)
               coord_sbuf(3*ii - 2) = u(1)
               coord_sbuf(3*ii - 1) = u(2)
               coord_sbuf(3*ii) = u(3)
               send_pos(dest + 1) = ii + 1
            END DO
         END DO
      END DO

      CALL pw_grid%para%group%alltoall(sbuf, scount, sdispl, rbuf, rcount, rdispl)
      CALL pw_grid%para%group%alltoall(coord_sbuf, coord_scount, coord_sdispl, &
                                       coord_rbuf, coord_rcount, coord_rdispl)

      DO ii = 1, nrecv
         u(1) = coord_rbuf(3*ii - 2)
         u(2) = coord_rbuf(3*ii - 1)
         u(3) = coord_rbuf(3*ii)
         IF (u(1) >= lb(1) .AND. u(1) <= ub(1) .AND. u(2) >= lb(2) .AND. u(2) <= ub(2) .AND. &
             u(3) >= lb(3) .AND. u(3) <= ub(3)) THEN
            density%array(u(1), u(2), u(3)) = rbuf(ii)
         END IF
      END DO

      DEALLOCATE (sbuf, rbuf, coord_sbuf, coord_rbuf, coord_scount, coord_sdispl, coord_rcount, &
                  coord_rdispl, scount, sdispl, rcount, rdispl, send_pos)

   END SUBROUTINE z_slices_to_cp2k_distribution_permuted

! **************************************************************************************************
!> \brief ...
!> \param wavelet ...
!> \param pw_grid ...
! **************************************************************************************************
   SUBROUTINE ps_wavelet_solve(wavelet, pw_grid)

      TYPE(ps_wavelet_type), POINTER                     :: wavelet
      TYPE(pw_grid_type), POINTER                        :: pw_grid

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ps_wavelet_solve'

      CHARACTER(LEN=1)                                   :: geocode
      INTEGER                                            :: handle, iproc, nproc, nx, ny, nz
      REAL(KIND=dp)                                      :: hx, hy, hz

      CALL timeset(routineN, handle)
      nproc = PRODUCT(pw_grid%para%group%num_pe_cart)
      iproc = pw_grid%para%group%mepos
      geocode = wavelet%geocode
      CALL get_wavelet_grid(wavelet, pw_grid, nx, ny, nz, hx, hy, hz)

      CALL PSolver(geocode, iproc, nproc, nx, ny, nz, hx, hy, hz, &
                   wavelet%rho_z_sliced, wavelet%karray, pw_grid)
      CALL timestop(handle)
   END SUBROUTINE ps_wavelet_solve

END MODULE ps_wavelet_methods
