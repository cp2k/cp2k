!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Transfers densities from PW to RS grids and potentials from PW to RS
!> \par History
!>      - Copied from qs_coolocate_Density and qs_integrate_potenntial
!> \author JGH (04.2014)
! **************************************************************************************************
MODULE rs_pw_interface
   USE cp_log_handling,                 ONLY: cp_to_string
   USE cp_spline_utils,                 ONLY: pw_interp,&
                                              pw_prolongate_s3,&
                                              pw_restrict_s3,&
                                              spline3_pbc_interp
   USE gaussian_gridlevels,             ONLY: gridlevel_info_type
   USE input_section_types,             ONLY: section_vals_val_get
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: FULLSPACE,&
                                              pw_grid_type,&
                                              pw_nyquist_map_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                              pw_pools_create_pws,&
                                              pw_pools_give_back_pws
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                              realspace_grid_type,&
                                              transfer_pw2rs,&
                                              transfer_rs2pw
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rs_pw_interface'
! *** Public subroutines ***

   PUBLIC :: pw_real_transfer
   PUBLIC :: density_rs2pw, &
             potential_pw2rs

CONTAINS

! **************************************************************************************************
!> \brief Transfer a real Fourier field between compatible multigrids, or apply its transpose.
!>        An even coarse-grid Nyquist coefficient represents both signs independently in each
!>        split direction. Ordinary complex coefficient copying does not supply that interpolation.
!>        Use the same real weights in potential restriction to retain the density adjoint.
!>        Half-space and spherical grids retain their existing frequency selection.
!> \param source Input coefficients of a real field.
!> \param destination Output coefficients, on the same grid or a compatible finer/coarser grid.
!> \param additive Add rather than overwrite the destination.
! **************************************************************************************************
   SUBROUTINE pw_real_transfer(source, destination, additive)
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: source
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: destination
      LOGICAL, INTENT(IN), OPTIONAL                      :: additive

      COMPLEX(dp), ALLOCATABLE                           :: incoming(:), outgoing(:)
      INTEGER                                            :: i
      LOGICAL                                            :: prolongate, split(3)
      TYPE(pw_grid_type), POINTER                        :: coarse, fine

      IF (PRESENT(additive)) THEN
         IF (additive) THEN
            CALL pw_axpy(source, destination)
         ELSE
            CALL pw_copy(source, destination)
         END IF
      ELSE
         CALL pw_copy(source, destination)
      END IF
      IF (ASSOCIATED(source%pw_grid, destination%pw_grid)) RETURN
      IF (source%pw_grid%grid_span /= FULLSPACE .OR. destination%pw_grid%grid_span /= FULLSPACE) RETURN
      IF (source%pw_grid%spherical .OR. destination%pw_grid%spherical) RETURN
      prolongate = ALL(source%pw_grid%npts <= destination%pw_grid%npts)
      IF (prolongate) THEN
         coarse => source%pw_grid
         fine => destination%pw_grid
      ELSE
         coarse => destination%pw_grid
         fine => source%pw_grid
      END IF
      CPASSERT(ALL(coarse%npts <= fine%npts))
      split = MOD(coarse%npts, 2) == 0 .AND. coarse%npts < fine%npts
      IF (.NOT. ANY(split)) RETURN
      IF (coarse%real_interpolation%reference /= fine%id_nr) CALL prepare_real_interpolation(coarse, fine, split)
      ASSOCIATE (map => coarse%real_interpolation, group => coarse%para%group)
         ALLOCATE (outgoing(SIZE(map%coarse_index)), incoming(SIZE(map%fine_index)))
         IF (prolongate) THEN
            outgoing(:) = map%weight*source%array(map%coarse_index)
            CALL group%alltoall(outgoing, map%send_count, map%send_offset, &
                                incoming, map%receive_count, map%receive_offset)
            DO i = 1, SIZE(incoming)
               destination%array(map%fine_index(i)) = destination%array(map%fine_index(i)) + incoming(i)
            END DO
         ELSE
            incoming(:) = source%array(map%fine_index)
            CALL group%alltoall(incoming, map%receive_count, map%receive_offset, &
                                outgoing, map%send_count, map%send_offset)
            DO i = 1, SIZE(outgoing)
               destination%array(map%coarse_index(i)) = destination%array(map%coarse_index(i)) + map%weight(i)*outgoing(i)
            END DO
         END IF
      END ASSOCIATE
   END SUBROUTINE pw_real_transfer

! **************************************************************************************************
!> \brief Cache only boundary-plane communication for a compatible full-box real interpolation.
!> \param coarse Grid owning the plan.
!> \param fine Reference grid receiving the split frequencies.
!> \param split Directions with an even, strictly smaller coarse extent.
! **************************************************************************************************
   SUBROUTINE prepare_real_interpolation(coarse, fine, split)
      TYPE(pw_grid_type), INTENT(INOUT)                  :: coarse
      TYPE(pw_grid_type), INTENT(IN)                     :: fine
      LOGICAL, INTENT(IN)                                :: split(3)

      INTEGER                                            :: axis, destination, i, ig, mask, &
                                                            nboundary, nproc, pass, slot, &
                                                            table_size, v(3)
      INTEGER(int_8)                                     :: key
      INTEGER(int_8), ALLOCATABLE                        :: keys(:)
      INTEGER, ALLOCATABLE                               :: coordinates(:), cursor(:), lookup(:), &
                                                            owner(:, :), received(:)
      LOGICAL                                            :: boundary(3)
      REAL(dp)                                           :: weight

      CPASSERT(coarse%para%group%num_pe == fine%para%group%num_pe)
      CPASSERT(coarse%para%group%mepos == fine%para%group%mepos)
      coarse%real_interpolation = pw_nyquist_map_type()
      nproc = coarse%para%group%num_pe
      ASSOCIATE (map => coarse%real_interpolation)
         ALLOCATE (map%send_count(0:nproc - 1), map%send_offset(0:nproc - 1), &
                   map%receive_count(0:nproc - 1), map%receive_offset(0:nproc - 1), cursor(0:nproc - 1))
         ALLOCATE (owner(fine%npts(2), fine%npts(3)))
         owner = 0
         ! Full-box reciprocal distributions own entire x rays. Infer their owners
         ! from actual Miller indices, independently of FFT remapping and blocking.
         DO ig = 1, SIZE(fine%g_hat, 2)
            v = MODULO(fine%g_hat(:, ig), fine%npts) + 1
            owner(v(2), v(3)) = fine%para%group%mepos + 1
         END DO
         CALL fine%para%group%max(owner)
         map%send_count = 0
         DO pass = 1, 2
            DO ig = 1, SIZE(coarse%g_hat, 2)
               boundary = split .AND. coarse%g_hat(:, ig) == -coarse%npts/2
               IF (.NOT. ANY(boundary)) CYCLE
               weight = 0.5_dp**COUNT(boundary)
               masks: DO mask = 0, 7
                  v = coarse%g_hat(:, ig)
                  DO axis = 1, 3
                     IF (.NOT. BTEST(mask, axis - 1)) CYCLE
                     IF (.NOT. boundary(axis)) CYCLE masks
                     v(axis) = -v(axis)
                  END DO
                  destination = owner(MODULO(v(2), fine%npts(2)) + 1, MODULO(v(3), fine%npts(3)) + 1) - 1
                  CPASSERT(destination >= 0 .AND. destination < nproc)
                  IF (pass == 1) THEN
                     map%send_count(destination) = map%send_count(destination) + 1
                  ELSE
                     slot = cursor(destination) + 1
                     cursor(destination) = slot
                     map%coarse_index(slot) = ig
                     map%weight(slot) = weight
                     IF (mask == 0) map%weight(slot) = weight - 1.0_dp
                     coordinates(3*slot - 2:3*slot) = v
                  END IF
               END DO masks
            END DO
            IF (pass == 1) THEN
               CALL coarse%para%group%alltoall(map%send_count, map%receive_count, 1)
               map%send_offset(0) = 0
               map%receive_offset(0) = 0
               DO i = 1, nproc - 1
                  map%send_offset(i) = map%send_offset(i - 1) + map%send_count(i - 1)
                  map%receive_offset(i) = map%receive_offset(i - 1) + map%receive_count(i - 1)
               END DO
               ALLOCATE (map%coarse_index(SUM(map%send_count)), map%weight(SUM(map%send_count)), &
                         coordinates(3*SUM(map%send_count)), received(3*SUM(map%receive_count)), &
                         map%fine_index(SUM(map%receive_count)))
               cursor(:) = map%send_offset
            END IF
         END DO
         CALL coarse%para%group%alltoall(coordinates, 3*map%send_count, 3*map%send_offset, &
                                         received, 3*map%receive_count, 3*map%receive_offset)
         nboundary = 0
         DO ig = 1, SIZE(fine%g_hat, 2)
            IF (ANY(split .AND. ABS(fine%g_hat(:, ig)) == coarse%npts/2)) nboundary = nboundary + 1
         END DO
         table_size = 2*nboundary + 1
         ALLOCATE (keys(table_size), lookup(table_size))
         keys = -1_int_8
         lookup = 0
         DO ig = 1, SIZE(fine%g_hat, 2)
            IF (.NOT. ANY(split .AND. ABS(fine%g_hat(:, ig)) == coarse%npts/2)) CYCLE
            key = real_interpolation_key(fine%g_hat(:, ig), fine)
            slot = INT(MODULO(key, INT(table_size, int_8))) + 1
            DO WHILE (keys(slot) /= -1_int_8)
               slot = MODULO(slot, table_size) + 1
            END DO
            keys(slot) = key
            lookup(slot) = ig
         END DO
         DO i = 1, SIZE(map%fine_index)
            key = real_interpolation_key(received(3*i - 2:3*i), fine)
            slot = INT(MODULO(key, INT(table_size, int_8))) + 1
            DO WHILE (keys(slot) /= key)
               CPASSERT(keys(slot) /= -1_int_8)
               slot = MODULO(slot, table_size) + 1
            END DO
            map%fine_index(i) = lookup(slot)
         END DO
         map%reference = fine%id_nr
      END ASSOCIATE
   END SUBROUTINE prepare_real_interpolation

! **************************************************************************************************
!> \brief Collision-free integer key for a full-box Miller index.
!> \param v Miller index.
!> \param grid Full reciprocal box.
!> \return Nonnegative lexicographic index.
! **************************************************************************************************
   PURE FUNCTION real_interpolation_key(v, grid) RESULT(key)
      INTEGER, INTENT(IN)                                :: v(3)
      TYPE(pw_grid_type), INTENT(IN)                     :: grid
      INTEGER(int_8)                                     :: key

      key = INT(v(1) - grid%bounds(1, 1), int_8) + INT(grid%npts(1), int_8)* &
            (INT(v(2) - grid%bounds(1, 2), int_8) + &
             INT(grid%npts(2), int_8)*INT(v(3) - grid%bounds(1, 3), int_8))
   END FUNCTION real_interpolation_key

! **************************************************************************************************
!> \brief given partial densities on the realspace multigrids,
!>      computes the full density on the plane wave grids, both in real and
!>      gspace
!> \param pw_env ...
!> \param rs_rho ...
!> \param rho ...
!> \param rho_gspace ...
!> \note
!>      should contain all communication in the collocation of the density
!>      in the case of replicated grids
! **************************************************************************************************
   SUBROUTINE density_rs2pw(pw_env, rs_rho, rho, rho_gspace)

      TYPE(pw_env_type), INTENT(IN)                      :: pw_env
      TYPE(realspace_grid_type), DIMENSION(:), &
         INTENT(IN)                                      :: rs_rho
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_gspace

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'density_rs2pw'

      INTEGER                                            :: handle, igrid_level, interp_kind
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(pw_c1d_gs_type), ALLOCATABLE, DIMENSION(:)    :: mgrid_gspace
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_r3d_rs_type), ALLOCATABLE, DIMENSION(:)    :: mgrid_rspace
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         POINTER                                         :: rs_descs

      CALL timeset(routineN, handle)
      NULLIFY (gridlevel_info, rs_descs, pw_pools)
      CALL pw_env_get(pw_env, rs_descs=rs_descs, pw_pools=pw_pools)

      gridlevel_info => pw_env%gridlevel_info

      CALL section_vals_val_get(pw_env%interp_section, "KIND", i_val=interp_kind)

      CALL pw_pools_create_pws(pw_pools, mgrid_rspace)

      CALL pw_pools_create_pws(pw_pools, mgrid_gspace)

      IF (gridlevel_info%ngrid_levels == 1) THEN
         CALL transfer_rs2pw(rs_rho(1), rho)
         CALL pw_transfer(rho, rho_gspace)
         IF (rho%pw_grid%spherical) THEN ! rho_gspace = rho
            CALL pw_transfer(rho_gspace, rho)
         END IF
      ELSE
         DO igrid_level = 1, gridlevel_info%ngrid_levels
            CALL transfer_rs2pw(rs_rho(igrid_level), &
                                mgrid_rspace(igrid_level))
         END DO

         ! we want both rho and rho_gspace, the latter for Hartree and co-workers.
         SELECT CASE (interp_kind)
         CASE (pw_interp)
            CALL pw_zero(rho_gspace)
            DO igrid_level = 1, gridlevel_info%ngrid_levels
               CALL pw_transfer(mgrid_rspace(igrid_level), &
                                mgrid_gspace(igrid_level))
               IF (pw_env%real_grid_interpolation) THEN
                  CALL pw_real_transfer(mgrid_gspace(igrid_level), rho_gspace, additive=.TRUE.)
               ELSE
                  CALL pw_axpy(mgrid_gspace(igrid_level), rho_gspace)
               END IF
            END DO
            CALL pw_transfer(rho_gspace, rho)
         CASE (spline3_pbc_interp)
            DO igrid_level = gridlevel_info%ngrid_levels, 2, -1
               CALL pw_prolongate_s3(mgrid_rspace(igrid_level), &
                                     mgrid_rspace(igrid_level - 1), pw_pools(igrid_level)%pool, &
                                     pw_env%interp_section)
            END DO
            CALL pw_copy(mgrid_rspace(1), rho)
            CALL pw_transfer(rho, rho_gspace)
         CASE default
            CALL cp_abort(__LOCATION__, &
                          "interpolator "// &
                          cp_to_string(interp_kind))
         END SELECT
      END IF

      ! *** give back the pw multi-grids
      CALL pw_pools_give_back_pws(pw_pools, mgrid_gspace)
      CALL pw_pools_give_back_pws(pw_pools, mgrid_rspace)
      CALL timestop(handle)

   END SUBROUTINE density_rs2pw

! **************************************************************************************************
!> \brief transfers a potential from a pw_grid to a vector of
!>      realspace multigrids
!> \param rs_v OUTPUT: the potential on the realspace multigrids
!> \param v_rspace INPUT : the potential on a planewave grid in Rspace
!> \param pw_env ...
!> \par History
!>      09.2006 created [Joost VandeVondele]
!> \note
!>      extracted from integrate_v_rspace
!>      should contain all parallel communication of integrate_v_rspace in the
!>      case of replicated grids.
! **************************************************************************************************
   SUBROUTINE potential_pw2rs(rs_v, v_rspace, pw_env)

      TYPE(realspace_grid_type), DIMENSION(:), &
         INTENT(IN)                                      :: rs_v
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: v_rspace
      TYPE(pw_env_type), INTENT(IN)                      :: pw_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'potential_pw2rs'

      INTEGER                                            :: auxbas_grid, handle, igrid_level, &
                                                            interp_kind
      REAL(KIND=dp)                                      :: scale
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(pw_c1d_gs_type), ALLOCATABLE, DIMENSION(:)    :: mgrid_gspace
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_r3d_rs_type), ALLOCATABLE, DIMENSION(:)    :: mgrid_rspace

      CALL timeset(routineN, handle)

      ! *** set up of the potential on the multigrids
      CALL pw_env_get(pw_env, pw_pools=pw_pools, gridlevel_info=gridlevel_info, &
                      auxbas_grid=auxbas_grid)

      CALL pw_pools_create_pws(pw_pools, mgrid_rspace)

      ! use either realspace or fft techniques to get the potential on the rs multigrids
      CALL section_vals_val_get(pw_env%interp_section, "KIND", i_val=interp_kind)
      SELECT CASE (interp_kind)
      CASE (pw_interp)
         CALL pw_pools_create_pws(pw_pools, mgrid_gspace)
         CALL pw_transfer(v_rspace, mgrid_gspace(auxbas_grid))
         DO igrid_level = 1, gridlevel_info%ngrid_levels
            IF (igrid_level /= auxbas_grid) THEN
               IF (pw_env%real_grid_interpolation) THEN
                  CALL pw_real_transfer(mgrid_gspace(auxbas_grid), mgrid_gspace(igrid_level))
               ELSE
                  CALL pw_copy(mgrid_gspace(auxbas_grid), mgrid_gspace(igrid_level))
               END IF
               CALL pw_transfer(mgrid_gspace(igrid_level), mgrid_rspace(igrid_level))
            ELSE
               IF (mgrid_gspace(auxbas_grid)%pw_grid%spherical) THEN
                  CALL pw_transfer(mgrid_gspace(auxbas_grid), mgrid_rspace(auxbas_grid))
               ELSE ! fft forward + backward should be identical
                  CALL pw_copy(v_rspace, mgrid_rspace(auxbas_grid))
               END IF
            END IF
            ! *** Multiply by the grid volume element ratio ***
            IF (igrid_level /= auxbas_grid) THEN
               scale = mgrid_rspace(igrid_level)%pw_grid%dvol/ &
                       mgrid_rspace(auxbas_grid)%pw_grid%dvol
               mgrid_rspace(igrid_level)%array = &
                  scale*mgrid_rspace(igrid_level)%array
            END IF
         END DO
         CALL pw_pools_give_back_pws(pw_pools, mgrid_gspace)
      CASE (spline3_pbc_interp)
         CALL pw_copy(v_rspace, mgrid_rspace(1))
         DO igrid_level = 1, gridlevel_info%ngrid_levels - 1
            CALL pw_zero(mgrid_rspace(igrid_level + 1))
            CALL pw_restrict_s3(mgrid_rspace(igrid_level), &
                                mgrid_rspace(igrid_level + 1), pw_pools(igrid_level + 1)%pool, &
                                pw_env%interp_section)
            ! *** Multiply by the grid volume element ratio
            mgrid_rspace(igrid_level + 1)%array = &
               mgrid_rspace(igrid_level + 1)%array*8._dp
         END DO
      CASE default
         CALL cp_abort(__LOCATION__, &
                       "interpolation not supported "// &
                       cp_to_string(interp_kind))
      END SELECT

      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL transfer_pw2rs(rs_v(igrid_level), &
                             mgrid_rspace(igrid_level))
      END DO
      ! *** give back the pw multi-grids
      CALL pw_pools_give_back_pws(pw_pools, mgrid_rspace)

      CALL timestop(handle)

   END SUBROUTINE potential_pw2rs

END MODULE rs_pw_interface
