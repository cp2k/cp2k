!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to Maxwell equation solver
!> \par History
!>      11/2020 created [mbrehm]
!> \author Martin Brehm
! **************************************************************************************************
MODULE maxwell_solver_interface
   USE cp_control_types, ONLY: maxwell_control_type
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_get_default_io_unit, &
                              cp_logger_type
   USE kinds, ONLY: dp
   USE pw_types, ONLY: pw_p_type
   USE message_passing, ONLY: &
      file_amode_rdonly, file_offset, mp_bcast, mp_file_close, mp_file_descriptor_type, &
      mp_file_get_position, mp_file_open, mp_file_read_all_chv, mp_file_type_free, &
      mp_file_type_hindexed_make_chv, mp_file_type_set_view_chv, mp_file_write_all_chv, &
      mp_file_write_at, mp_maxloc, mp_recv, mp_send, mp_sum, mp_sync, mpi_character_size
   USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE

#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'maxwell_solver_interface'

! *** Public subroutines ***
   PUBLIC :: maxwell_solver

#if defined(__LIBMAXWELL)

   INTERFACE

      INTEGER(C_INT) FUNCTION libcp2kmw_setgrid(rx, ry, rz, ax, ay, az, bx, by, bz, cx, cy, cz) BIND(C, NAME='libcp2kmw_setgrid')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
         IMPLICIT NONE
         INTEGER(C_INT) :: rx, ry, rz
         REAL(C_DOUBLE) :: ax, ay, az, bx, by, bz, cx, cy, cz
      END FUNCTION libcp2kmw_setgrid

      INTEGER(C_INT) FUNCTION libcp2kmw_step(step, t) BIND(C, NAME='libcp2kmw_step')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
         IMPLICIT NONE
         INTEGER(C_INT) :: step
         REAL(C_DOUBLE) :: t
      END FUNCTION libcp2kmw_step

      INTEGER(C_INT) FUNCTION libcp2kmw_getzrow(buf, px, py, zmin, zmax) BIND(C, NAME='libcp2kmw_getzrow')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
         IMPLICIT NONE
         REAL(C_DOUBLE) :: buf(*)
         INTEGER(C_INT) :: px, py, zmin, zmax
      END FUNCTION libcp2kmw_getzrow

   END INTERFACE

#endif

CONTAINS

! **************************************************************************************************
!> \brief  Computes the external potential on the grid
!> \param maxwell_control the Maxwell control section
!> \param v_ee the realspace grid with the potential
!> \param sim_step current simulation step
!> \param sim_time current physical simulation time
!> \param scaling_factor a factor to scale the potential with
!> \date   12/2020
!> \author Martin Brehm
! **************************************************************************************************
   SUBROUTINE maxwell_solver(maxwell_control, v_ee, sim_step, sim_time, scaling_factor)
      TYPE(maxwell_control_type), INTENT(IN)             :: maxwell_control
      TYPE(pw_p_type), POINTER                           :: v_ee
      INTEGER, INTENT(IN)                                :: sim_step
      REAL(KIND=dp), INTENT(IN)                          :: sim_time
      REAL(KIND=dp), INTENT(IN)                          :: scaling_factor

#if defined(__LIBMAXWELL)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'maxwell_solver'

      INTEGER                                            :: handle, iounit, res, my_rank, num_pe, &
                                                            gid, master, tag, i, j, ip
      TYPE(cp_logger_type), POINTER                      :: logger

      INTEGER, DIMENSION(3)                              :: lbounds, lbounds_local, npoints, &
                                                            npoints_local, ubounds, ubounds_local
      REAL(C_DOUBLE), ALLOCATABLE, DIMENSION(:)          :: buffer

      MARK_USED(maxwell_control)
      MARK_USED(v_ee)
      MARK_USED(sim_step)
      MARK_USED(sim_time)

      CALL timeset(routineN, handle)
      NULLIFY (logger)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      my_rank = v_ee%pw%pw_grid%para%my_pos
      num_pe = v_ee%pw%pw_grid%para%group_size
      gid = v_ee%pw%pw_grid%para%group
      tag = 1

      lbounds = v_ee%pw%pw_grid%bounds(1, :)
      ubounds = v_ee%pw%pw_grid%bounds(2, :)
      npoints = v_ee%pw%pw_grid%npts

      lbounds_local = v_ee%pw%pw_grid%bounds_local(1, :)
      ubounds_local = v_ee%pw%pw_grid%bounds_local(2, :)
      npoints_local = v_ee%pw%pw_grid%npts_local

      ALLOCATE (buffer(lbounds(3):ubounds(3)))

      IF (my_rank == 0) THEN

         IF (iounit > 0) THEN
            WRITE (iounit, *) ""
            WRITE (iounit, *) "MAXWELL| Called, step = ", sim_step, " time = ", sim_time
         ENDIF

         res = libcp2kmw_setgrid( &
               ubounds(1) - lbounds(1) + 1, &
               ubounds(2) - lbounds(2) + 1, &
               ubounds(3) - lbounds(3) + 1, &
               v_ee%pw%pw_grid%dh(1, 1)*(ubounds(1) - lbounds(1) + 1), &
               v_ee%pw%pw_grid%dh(2, 1)*(ubounds(1) - lbounds(1) + 1), &
               v_ee%pw%pw_grid%dh(3, 1)*(ubounds(1) - lbounds(1) + 1), &
               v_ee%pw%pw_grid%dh(1, 2)*(ubounds(2) - lbounds(2) + 1), &
               v_ee%pw%pw_grid%dh(2, 2)*(ubounds(2) - lbounds(2) + 1), &
               v_ee%pw%pw_grid%dh(3, 2)*(ubounds(2) - lbounds(2) + 1), &
               v_ee%pw%pw_grid%dh(1, 3)*(ubounds(3) - lbounds(3) + 1), &
               v_ee%pw%pw_grid%dh(2, 3)*(ubounds(3) - lbounds(3) + 1), &
               v_ee%pw%pw_grid%dh(3, 3)*(ubounds(3) - lbounds(3) + 1) &
               )

         res = libcp2kmw_step(sim_step, sim_time)

         IF (iounit > 0) THEN
            WRITE (iounit, *) "MAXWELL| Returned with value ", res
            WRITE (iounit, *) "MAXWELL| Distributing potential to MPI processes..."
         ENDIF

      ENDIF

      ! The following code block is copied from src/pw/realspace_grid_cube.F

      !master sends all data to everyone
      DO i = lbounds(1), ubounds(1)
         DO j = lbounds(2), ubounds(2)

            IF (my_rank == 0) THEN
               !buffer(lbounds(3)) = 23.0;
               !buffer(lbounds(3)+1) = 42.0;
               !buffer(lbounds(3)+2) = 69.0;
               res = libcp2kmw_getzrow(buffer, i - lbounds(1) + 1, j - lbounds(2) + 1, 1, ubounds(3) - lbounds(3) + 1)
               IF (num_pe > 1) THEN
                  DO ip = 1, num_pe - 1
                     CALL mp_send(buffer(lbounds(3):ubounds(3)), ip, tag, gid)
                  ENDDO
               ENDIF
            ELSE
               master = 0
               CALL mp_recv(buffer(lbounds(3):ubounds(3)), master, tag, gid)
            ENDIF

            CALL mp_sync(gid)

            !only use data that is local to me - i.e. in slice of pencil I own
            IF ((lbounds_local(1) <= i) .AND. (i <= ubounds_local(1)) .AND. &
                (lbounds_local(2) <= j) .AND. (j <= ubounds_local(2))) THEN
               !allow scaling of external potential values by factor 'scaling' (SCALING_FACTOR in input file)
               v_ee%pw%cr3d(i, j, lbounds(3):ubounds(3)) = buffer(lbounds(3):ubounds(3))*scaling_factor
            ENDIF

         ENDDO
      ENDDO

      IF (iounit > 0) THEN
         WRITE (iounit, *) "MAXWELL| All done."
      ENDIF

      CALL timestop(handle)

#else

      MARK_USED(maxwell_control)
      MARK_USED(v_ee)
      MARK_USED(sim_step)
      MARK_USED(sim_time)
      MARK_USED(scaling_factor)

      CALL cp_abort(__LOCATION__, &
                    "The Maxwell solver interface requires CP2k to be compiled &
                     &with the -D__LIBMAXWELL preprocessor option.")

#endif

   END SUBROUTINE maxwell_solver

END MODULE maxwell_solver_interface

