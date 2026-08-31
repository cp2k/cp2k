!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Routines to handle an external density
!>        The external density can be generic and is provided by user input
!> \author D. Varsano
! **************************************************************************************************
MODULE qs_external_density
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE gaussian_gridlevels,             ONLY: gridlevel_info_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_integrate_function
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                              realspace_grid_type,&
                                              rs_grid_create,&
                                              rs_grid_release,&
                                              rs_grid_zero
   USE rs_pw_interface,                 ONLY: density_rs2pw
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_external_density'

   PUBLIC :: external_read_density, read_cube_density

CONTAINS

! **************************************************************************************************
!> \brief  Computes the external density on the grid
!> \param qs_env ...
!> \date   03.2011
!> \author D. Varsano
! **************************************************************************************************
   SUBROUTINE external_read_density(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'external_read_density'

      CHARACTER(LEN=default_string_length)               :: filename
      INTEGER                                            :: handle
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_rho_type), POINTER                         :: rho_external
      TYPE(section_vals_type), POINTER                   :: ext_den_section, input

      CALL timeset(routineN, handle)
      NULLIFY (input, ext_den_section, dft_control, rho_external)

      CALL get_qs_env(qs_env, &
                      rho_external=rho_external, &
                      input=input, &
                      dft_control=dft_control)

      IF (dft_control%apply_external_density .AND. dft_control%read_external_density) THEN
         ext_den_section => section_vals_get_subs_vals(input, "DFT%EXTERNAL_DENSITY")
         CALL section_vals_val_get(ext_den_section, "FILE_DENSITY", c_val=filename)
         CALL read_cube_density(qs_env, rho_external, TRIM(filename), &
                                total_density_sign=1, source_label="ZMP")
      END IF

      CALL timestop(handle)

   END SUBROUTINE external_read_density

! **************************************************************************************************
!> \brief Read an electron density from a Gaussian cube file into a QS density grid
!> \param qs_env ...
!> \param rho_target target density structure
!> \param filename cube filename
!> \param total_density_sign sign used for rho_target%tot_rho_r
!> \param source_label label used in output
! **************************************************************************************************
   SUBROUTINE read_cube_density(qs_env, rho_target, filename, total_density_sign, source_label)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_target
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      INTEGER, INTENT(IN)                                :: total_density_sign
      CHARACTER(LEN=*), INTENT(IN)                       :: source_label

      CHARACTER(len=*), PARAMETER                        :: routineN = 'read_cube_density'

      INTEGER                                            :: extunit, handle, i, igrid_level, j, k, &
                                                            nat, ndum
      INTEGER, DIMENSION(3)                              :: lbounds, lbounds_local, npoints, &
                                                            ubounds, ubounds_local
      LOGICAL                                            :: grid_mismatch
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: buffer
      REAL(kind=dp), DIMENSION(3)                        :: cube_origin, voxel
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         POINTER                                         :: rs_descs
      TYPE(realspace_grid_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: rs_rho

      CALL timeset(routineN, handle)
      NULLIFY (pw_env, rho_r, rho_g, tot_rho_r, rs_descs)

      IF (total_density_sign /= -1 .AND. total_density_sign /= 1) THEN
         CPABORT("total_density_sign has to be -1 or 1")
      END IF

      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL qs_rho_get(rho_target, rho_r=rho_r, rho_g=rho_g, tot_rho_r=tot_rho_r)
      gridlevel_info => pw_env%gridlevel_info
      CALL pw_env_get(pw_env, rs_descs=rs_descs)

      ALLOCATE (rs_rho(gridlevel_info%ngrid_levels))
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL rs_grid_create(rs_rho(igrid_level), rs_descs(igrid_level)%rs_desc)
         CALL rs_grid_zero(rs_rho(igrid_level))
      END DO
      igrid_level = gridlevel_info%ngrid_levels

      npoints = rs_descs(igrid_level)%rs_desc%npts
      lbounds = rs_descs(igrid_level)%rs_desc%lb
      ubounds = rs_descs(igrid_level)%rs_desc%ub
      lbounds_local = rho_r(1)%pw_grid%bounds_local(1, :)
      ubounds_local = rho_r(1)%pw_grid%bounds_local(2, :)
      ALLOCATE (buffer(lbounds(3):ubounds(3)))

      ASSOCIATE (gid => rho_r(1)%pw_grid%para%group, &
                 my_rank => rho_r(1)%pw_grid%para%group%mepos)
         grid_mismatch = .FALSE.
         IF (my_rank == 0) THEN
            WRITE (*, FMT="(/,T3,A,A)") TRIM(source_label)//"| Reading electron density: ", TRIM(filename)
            CALL open_file(file_name=filename, file_status="OLD", file_form="FORMATTED", &
                           file_action="READ", unit_number=extunit)

            READ (extunit, *)
            READ (extunit, *)
            READ (extunit, *) nat, cube_origin
            IF (MAXVAL(ABS(cube_origin)) > 1.0E-4_dp) THEN
               grid_mismatch = .TRUE.
               WRITE (*, FMT="(T3,A,3ES16.8)") TRIM(source_label)// &
                  "| Cube origin is not the CP2K grid origin: ", cube_origin
            END IF
            IF (nat < 0) THEN
               grid_mismatch = .TRUE.
               WRITE (*, FMT="(T3,A)") TRIM(source_label)// &
                  "| Multi-orbital cube files are not supported"
            END IF
            DO i = 1, 3
               READ (extunit, *) ndum, voxel
               IF (ndum /= npoints(i) .OR. &
                   MAXVAL(ABS(voxel - rs_descs(igrid_level)%rs_desc%dh(:, i))) > 1.0E-4_dp) THEN
                  grid_mismatch = .TRUE.
                  WRITE (*, FMT="(T3,A,I0)") TRIM(source_label)// &
                     "| Cube grid does not coincide with CP2K grid along axis ", i
                  WRITE (*, FMT="(T3,A,I0,A,I0)") TRIM(source_label)//"| Grid points: ", &
                     ndum, " instead of ", npoints(i)
                  WRITE (*, FMT="(T3,A,3ES16.8)") TRIM(source_label)//"| Cube vector: ", voxel
                  WRITE (*, FMT="(T3,A,3ES16.8)") TRIM(source_label)//"| CP2K vector: ", &
                     rs_descs(igrid_level)%rs_desc%dh(:, i)
               END IF
            END DO
            DO i = 1, ABS(nat)
               READ (extunit, *)
            END DO
         END IF

         CALL gid%bcast(grid_mismatch, 0)
         IF (grid_mismatch) THEN
            IF (my_rank == 0) CALL close_file(unit_number=extunit)
            CPABORT("Cube density grid is incompatible with the CP2K real-space grid")
         END IF

         DO i = lbounds(1), ubounds(1)
            DO j = lbounds(2), ubounds(2)
               IF (my_rank == 0) THEN
                  READ (extunit, *) (buffer(k), k=lbounds(3), ubounds(3))
               END IF
               CALL gid%bcast(buffer, 0)
               IF ((lbounds_local(1) <= i) .AND. (i <= ubounds_local(1)) .AND. &
                   (lbounds_local(2) <= j) .AND. (j <= ubounds_local(2))) THEN
                  rs_rho(igrid_level)%r(i, j, lbounds(3):ubounds(3)) = buffer
               END IF
            END DO
         END DO
         IF (my_rank == 0) CALL close_file(unit_number=extunit)

         CALL density_rs2pw(pw_env, rs_rho, rho=rho_r(1), rho_gspace=rho_g(1))
         tot_rho_r(1) = pw_integrate_function(rho_r(1), isign=total_density_sign)
         IF (my_rank == 0) THEN
            WRITE (*, FMT="(T3,A,T61,F20.10)") TRIM(source_label)// &
               "| Integrated electron density:", REAL(total_density_sign, dp)*tot_rho_r(1)
         END IF
         CALL gid%sync()
      END ASSOCIATE

      DO igrid_level = 1, SIZE(rs_rho)
         CALL rs_grid_release(rs_rho(igrid_level))
      END DO
      DEALLOCATE (buffer, rs_rho)

      CALL timestop(handle)

   END SUBROUTINE read_cube_density

END MODULE qs_external_density
