!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface for Voronoi Integration and output of BQB files
!> \par History
!>      11/2020 created [mbrehm]
!> \author Martin Brehm
! **************************************************************************************************
MODULE voronoi_interface
   USE input_section_types, ONLY: section_vals_type, &
                                  section_vals_val_get
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_get_default_io_unit, &
                              cp_logger_type
   USE bibliography, ONLY: Rycroft2009, Thomas2015, Brehm2018, Brehm2020, &
                           Brehm2021, cite_reference
   USE kinds, ONLY: dp, default_path_length
   USE cell_types, ONLY: cell_type
   USE pw_types, ONLY: pw_p_type, pw_type
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE qs_rho_types, ONLY: qs_rho_get, &
                           qs_rho_type
   USE atomic_kind_types, ONLY: atomic_kind_type, &
                                get_atomic_kind
   USE particle_list_types, ONLY: particle_list_type
   USE particle_types, ONLY: particle_type
   USE qs_kind_types, ONLY: get_qs_kind, &
                            qs_kind_type
   USE cp_para_types, ONLY: cp_para_env_type
   USE qs_subsys_types, ONLY: qs_subsys_get, &
                              qs_subsys_type
   USE cp_control_types, ONLY: dft_control_type
   USE pw_grid_types, ONLY: PW_MODE_LOCAL
   USE physcon, ONLY: angstrom, femtoseconds
   USE message_passing, ONLY: &
      file_amode_rdonly, file_offset, mp_bcast, mp_file_close, mp_file_descriptor_type, &
      mp_file_get_position, mp_file_open, mp_file_read_all_chv, mp_file_type_free, &
      mp_file_type_hindexed_make_chv, mp_file_type_set_view_chv, mp_file_write_all_chv, &
      mp_file_write_at, mp_maxloc, mp_recv, mp_send, mp_sum, mp_sync, mpi_character_size
   USE input_constants, ONLY: &
      voro_radii_unity, voro_radii_vdw, voro_radii_cov, voro_radii_user, &
      bqb_opt_off, bqb_opt_quick, bqb_opt_normal, bqb_opt_patient, bqb_opt_exhaustive
   USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_CHAR
#if defined(__HAS_IEEE_EXCEPTIONS)
   USE ieee_exceptions, ONLY: ieee_get_halting_mode, &
                              ieee_set_halting_mode, &
                              ieee_all
#endif

#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Global parameters
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'voronoi_interface'
   PUBLIC :: entry_voronoi_or_bqb, finalize_libvori

#if defined(__LIBVORI)

   ! The C interface to libvori

   INTERFACE

      INTEGER(C_INT) FUNCTION libvori_setBQBSkipFirst(i) BIND(C, NAME='libvori_setBQBSkipFirst')
         USE ISO_C_BINDING, ONLY: C_INT
      INTEGER(C_INT), VALUE                              :: i

      END FUNCTION libvori_setBQBSkipFirst

      INTEGER(C_INT) FUNCTION libvori_setBQBStoreStep(i) BIND(C, NAME='libvori_setBQBStoreStep')
         USE ISO_C_BINDING, ONLY: C_INT
      INTEGER(C_INT), VALUE                              :: i

      END FUNCTION libvori_setBQBStoreStep

      INTEGER(C_INT) FUNCTION libvori_setVoronoiSkipFirst(i) BIND(C, NAME='libvori_setVoronoiSkipFirst')
         USE ISO_C_BINDING, ONLY: C_INT
      INTEGER(C_INT), VALUE                              :: i

      END FUNCTION libvori_setVoronoiSkipFirst

      INTEGER(C_INT) FUNCTION libvori_setBQBCheck(i) BIND(C, NAME='libvori_setBQBCheck')
         USE ISO_C_BINDING, ONLY: C_INT
      INTEGER(C_INT), VALUE                              :: i

      END FUNCTION libvori_setBQBCheck

      INTEGER(C_INT) FUNCTION libvori_setBQBFilename(len, s) BIND(C, NAME='libvori_setBQBFilename')
         USE ISO_C_BINDING, ONLY: C_INT, C_CHAR
      INTEGER(C_INT), VALUE                              :: len
      CHARACTER(C_CHAR)                                  :: s(*)

      END FUNCTION libvori_setBQBFilename

      INTEGER(C_INT) FUNCTION libvori_setBQBParmString(len, s) BIND(C, NAME='libvori_setBQBParmString')
         USE ISO_C_BINDING, ONLY: C_INT, C_CHAR
      INTEGER(C_INT), VALUE                              :: len
      CHARACTER(C_CHAR)                                  :: s(*)

      END FUNCTION libvori_setBQBParmString

      INTEGER(C_INT) FUNCTION libvori_setBQBHistory(i) BIND(C, NAME='libvori_setBQBHistory')
         USE ISO_C_BINDING, ONLY: C_INT
      INTEGER(C_INT), VALUE                              :: i

      END FUNCTION libvori_setBQBHistory

      INTEGER(C_INT) FUNCTION libvori_setBQBOptimization(i) BIND(C, NAME='libvori_setBQBOptimization')
         USE ISO_C_BINDING, ONLY: C_INT
      INTEGER(C_INT), VALUE                              :: i

      END FUNCTION libvori_setBQBOptimization

      INTEGER(C_INT) FUNCTION libvori_processBQBFrame(step, t) BIND(C, NAME='libvori_processBQBFrame')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), VALUE                              :: step
      REAL(C_DOUBLE), VALUE                              :: t

      END FUNCTION libvori_processBQBFrame

      INTEGER(C_INT) FUNCTION libvori_setPrefix_Voronoi() BIND(C, NAME='libvori_setPrefix_Voronoi')
         USE ISO_C_BINDING, ONLY: C_INT
      END FUNCTION libvori_setPrefix_Voronoi

      INTEGER(C_INT) FUNCTION libvori_setPrefix_BQB() BIND(C, NAME='libvori_setPrefix_BQB')
         USE ISO_C_BINDING, ONLY: C_INT
      END FUNCTION libvori_setPrefix_BQB

      INTEGER(C_INT) FUNCTION libvori_setRefinementFactor(i) BIND(C, NAME='libvori_setRefinementFactor')
         USE ISO_C_BINDING, ONLY: C_INT
      INTEGER(C_INT), VALUE                              :: i

      END FUNCTION libvori_setRefinementFactor

      INTEGER(C_INT) FUNCTION libvori_setEMPOutput(i) BIND(C, NAME='libvori_setEMPOutput')
         USE ISO_C_BINDING, ONLY: C_INT
      INTEGER(C_INT), VALUE                              :: i

      END FUNCTION libvori_setEMPOutput

      INTEGER(C_INT) FUNCTION libvori_setPrintLevel_Verbose() BIND(C, NAME='libvori_setPrintLevel_Verbose')
         USE ISO_C_BINDING, ONLY: C_INT
      END FUNCTION libvori_setPrintLevel_Verbose

      INTEGER(C_INT) FUNCTION libvori_setRadii_Unity() BIND(C, NAME='libvori_setRadii_Unity')
         USE ISO_C_BINDING, ONLY: C_INT
      END FUNCTION libvori_setRadii_Unity

      INTEGER(C_INT) FUNCTION libvori_setRadii_Covalent() BIND(C, NAME='libvori_setRadii_Covalent')
         USE ISO_C_BINDING, ONLY: C_INT
      END FUNCTION libvori_setRadii_Covalent

      INTEGER(C_INT) FUNCTION libvori_setRadii_User(factor, rad) BIND(C, NAME='libvori_setRadii_User')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      REAL(C_DOUBLE), VALUE                              :: factor
      REAL(C_DOUBLE)                                     :: rad(*)

      END FUNCTION libvori_setRadii_User

      INTEGER(C_INT) FUNCTION libvori_step(step, t) BIND(C, NAME='libvori_step')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), VALUE                              :: step
      REAL(C_DOUBLE), VALUE                              :: t

      END FUNCTION libvori_step

      INTEGER(C_INT) FUNCTION libvori_sanitycheck(step, t) BIND(C, NAME='libvori_sanitycheck')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), VALUE                              :: step
      REAL(C_DOUBLE), VALUE                              :: t

      END FUNCTION libvori_sanitycheck

      INTEGER(C_INT) FUNCTION libvori_setGrid(rx, ry, rz, ax, ay, az, bx, by, bz, cx, cy, cz, tax, tay, taz, tbx, tby, tbz, &
                                              tcx, tcy, tcz) BIND(C, NAME='libvori_setGrid')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), VALUE                              :: rx, ry, rz
      REAL(C_DOUBLE), VALUE                              :: ax, ay, az, bx, by, bz, cx, cy, cz, tax, &
                                                            tay, taz, tbx, tby, tbz, tcx, tcy, tcz

      END FUNCTION libvori_setGrid

      INTEGER(C_INT) FUNCTION libvori_pushAtoms(n, pord, pchg, posx, posy, posz) BIND(C, NAME='libvori_pushAtoms')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), VALUE                              :: n
      INTEGER(C_INT)                                     :: pord(*)
      REAL(C_DOUBLE)                                     :: pchg(*), posx(*), posy(*), posz(*)

      END FUNCTION libvori_pushAtoms

      INTEGER(C_INT) FUNCTION libvori_push_rho_zrow(ix, iy, length, buf) BIND(C, NAME='libvori_push_rho_zrow')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), VALUE                              :: ix, iy, length
      REAL(C_DOUBLE)                                     :: buf(*)

      END FUNCTION libvori_push_rho_zrow

      INTEGER(C_INT) FUNCTION libvori_setBQBOverwrite(i) BIND(C, NAME='libvori_setBQBOverwrite')
         USE ISO_C_BINDING, ONLY: C_INT
      INTEGER(C_INT), VALUE                              :: i

      END FUNCTION libvori_setBQBOverwrite

      INTEGER(C_INT) FUNCTION libvori_setVoriOverwrite(i) BIND(C, NAME='libvori_setVoriOverwrite')
         USE ISO_C_BINDING, ONLY: C_INT
      INTEGER(C_INT), VALUE                              :: i

      END FUNCTION libvori_setVoriOverwrite

      INTEGER(C_INT) FUNCTION libvori_get_radius(length, buf) BIND(C, NAME='libvori_get_radius')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), VALUE                              :: length
      REAL(C_DOUBLE)                                     :: buf(*)

      END FUNCTION libvori_get_radius

      INTEGER(C_INT) FUNCTION libvori_get_volume(length, buf) BIND(C, NAME='libvori_get_volume')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), VALUE                              :: length
      REAL(C_DOUBLE)                                     :: buf(*)

      END FUNCTION libvori_get_volume

      INTEGER(C_INT) FUNCTION libvori_get_charge(length, buf) BIND(C, NAME='libvori_get_charge')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), VALUE                              :: length
      REAL(C_DOUBLE)                                     :: buf(*)

      END FUNCTION libvori_get_charge

      INTEGER(C_INT) FUNCTION libvori_get_dipole(component, length, buf) BIND(C, NAME='libvori_get_dipole')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), VALUE                              :: component, length
      REAL(C_DOUBLE)                                     :: buf(*)

      END FUNCTION libvori_get_dipole

      INTEGER(C_INT) FUNCTION libvori_get_quadrupole(component, length, buf) BIND(C, NAME='libvori_get_quadrupole')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), VALUE                              :: component, length
      REAL(C_DOUBLE)                                     :: buf(*)

      END FUNCTION libvori_get_quadrupole

      INTEGER(C_INT) FUNCTION libvori_get_wrapped_pos(component, length, buf) BIND(C, NAME='libvori_get_wrapped_pos')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), VALUE                              :: component, length
      REAL(C_DOUBLE)                                     :: buf(*)

      END FUNCTION libvori_get_wrapped_pos

      INTEGER(C_INT) FUNCTION libvori_get_charge_center(component, length, buf) BIND(C, NAME='libvori_get_charge_center')
         USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      INTEGER(C_INT), VALUE                              :: component, length
      REAL(C_DOUBLE)                                     :: buf(*)

      END FUNCTION libvori_get_charge_center

      INTEGER(C_INT) FUNCTION libvori_finalize() BIND(C, NAME='libvori_finalize')
         USE ISO_C_BINDING, ONLY: C_INT
      END FUNCTION libvori_finalize
   END INTERFACE

#endif

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Does a Voronoi integration of density or stores the density to compressed BQB format
!> \param do_voro whether a Voronoi integration shall be performed
!> \param do_bqb whether the density shall be written to compressed BQB format
!> \param input_voro the input section for Voronoi integration
!> \param input_bqb the input section for the BQB compression
!> \param unit_voro the output unit number for the Voronoi integration results
!> \param qs_env the qs_env where to calculate the charges
!> \param rspace_pw the grid with the real-space electron density to integrate/compress
!> \author Martin Brehm
! **************************************************************************************************
   SUBROUTINE entry_voronoi_or_bqb(do_voro, do_bqb, input_voro, input_bqb, unit_voro, qs_env, rspace_pw)
      INTEGER                                            :: do_voro, do_bqb
      TYPE(section_vals_type), POINTER                   :: input_voro, input_bqb
      INTEGER, INTENT(IN)                                :: unit_voro
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_type), POINTER                             :: rspace_pw

#if defined(__LIBVORI)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'entry_voronoi_or_bqb'
      INTEGER                                            :: handle, iounit, L1, L2, L3, U1, U2, &
                                                            U3, ret, i, my_rank, num_pe, gid, tag, &
                                                            nkind, natom, ikind, iat, ord, source, dest, &
                                                            ip, i1, i2, reffac, radius_type, bqb_optimize, &
                                                            bqb_history
      LOGICAL                                            :: outemp, bqb_skip_first, voro_skip_first, &
                                                            bqb_store_step, bqb_check, voro_sanity, &
                                                            bqb_overwrite, vori_overwrite
      REAL(KIND=dp)                                      :: zeff
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(cp_logger_type), POINTER                      :: logger
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: buf
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: particles_z
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: particles_r
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: particles_c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: particles_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER                                            :: sim_step
      REAL(KIND=dp)                                      :: r, sim_time
      TYPE(qs_subsys_type), POINTER                      :: subsys
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: voro_radii
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: voro_charge
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: voro_volume
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: voro_dipole
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: voro_quadrupole
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: voro_buffer
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: voro_wrapped_pos
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: voro_charge_center
      REAL(KIND=dp), DIMENSION(:), POINTER               :: user_radii
      CHARACTER(len=default_path_length)                 :: bqb_file_name
      CHARACTER(len=128)                                 :: bqb_parm_string
#if defined(__HAS_IEEE_EXCEPTIONS)
      LOGICAL, DIMENSION(5)                              :: halt
#endif

      CALL timeset(routineN, handle)
      NULLIFY (logger)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env=qs_env, rho=rho, qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, para_env=para_env, &
                      nkind=nkind, natom=natom, subsys=subsys, dft_control=dft_control, &
                      cell=cell)

      my_rank = para_env%mepos
      num_pe = para_env%num_pe
      gid = para_env%group
      tag = 1

      sim_step = qs_env%sim_step
      sim_time = qs_env%sim_time

      IF (do_voro /= 0) THEN
         CALL section_vals_val_get(input_voro, "REFINEMENT_FACTOR", i_val=reffac)
         CALL section_vals_val_get(input_voro, "OUTPUT_EMP", l_val=outemp)
         CALL section_vals_val_get(input_voro, "VORONOI_RADII", i_val=radius_type)
         CALL section_vals_val_get(input_voro, "SKIP_FIRST", l_val=voro_skip_first)
         CALL section_vals_val_get(input_voro, "SANITY_CHECK", l_val=voro_sanity)
         CALL section_vals_val_get(input_voro, "OVERWRITE", l_val=vori_overwrite)
         IF (radius_type == voro_radii_user) THEN
            CALL section_vals_val_get(input_voro, "USER_RADII", r_vals=user_radii)
         END IF
         IF (qs_env%single_point_run) THEN
            voro_skip_first = .FALSE.
         END IF
         CALL cite_reference(Rycroft2009)
         CALL cite_reference(Thomas2015)
         CALL cite_reference(Brehm2018)
         CALL cite_reference(Brehm2020)
         CALL cite_reference(Brehm2021)
      END IF

      IF (do_bqb /= 0) THEN
         CALL section_vals_val_get(input_bqb, "HISTORY", i_val=bqb_history)
         CALL section_vals_val_get(input_bqb, "OPTIMIZE", i_val=bqb_optimize)
         CALL section_vals_val_get(input_bqb, "FILENAME", c_val=bqb_file_name)
         CALL section_vals_val_get(input_bqb, "SKIP_FIRST", l_val=bqb_skip_first)
         CALL section_vals_val_get(input_bqb, "STORE_STEP_NUMBER", l_val=bqb_store_step)
         CALL section_vals_val_get(input_bqb, "CHECK", l_val=bqb_check)
         CALL section_vals_val_get(input_bqb, "OVERWRITE", l_val=bqb_overwrite)
         CALL section_vals_val_get(input_bqb, "PARAMETER_KEY", c_val=bqb_parm_string)
         IF (qs_env%single_point_run) THEN
            bqb_skip_first = .FALSE.
            bqb_history = 1
         END IF
         IF (bqb_history < 1) THEN
            bqb_history = 1
         END IF
         CALL cite_reference(Brehm2018)
      END IF

      CALL qs_subsys_get(subsys, particles=particles)

      ! The boundaries for the grid index
      L1 = rspace_pw%pw_grid%bounds(1, 1)
      L2 = rspace_pw%pw_grid%bounds(1, 2)
      L3 = rspace_pw%pw_grid%bounds(1, 3)
      U1 = rspace_pw%pw_grid%bounds(2, 1)
      U2 = rspace_pw%pw_grid%bounds(2, 2)
      U3 = rspace_pw%pw_grid%bounds(2, 3)

      ! Temporarily disable floating point traps because libvori raise IEEE754 exceptions.
#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_get_halting_mode(IEEE_ALL, halt)
      CALL ieee_set_halting_mode(IEEE_ALL, .FALSE.)
#endif

      IF (my_rank == 0) THEN

         IF (iounit > 0) THEN
            WRITE (iounit, *) ""
         END IF

         IF (do_voro /= 0) THEN
            ret = libvori_setPrefix_Voronoi()
            ret = libvori_setRefinementFactor(reffac)
            ret = libvori_setVoronoiSkipFirst(MERGE(1, 0, voro_skip_first))
            ret = libvori_setVoriOverwrite(MERGE(1, 0, vori_overwrite))
            ret = libvori_setEMPOutput(MERGE(1, 0, outemp))
         ELSE
            ret = libvori_setPrefix_BQB()
         END IF

         IF (do_bqb /= 0) THEN
            ret = libvori_setBQBFilename(default_path_length, bqb_file_name)
            ret = libvori_setBQBParmString(128, bqb_parm_string)
            SELECT CASE (bqb_optimize)
            CASE (bqb_opt_off)
               bqb_optimize = 0
            CASE (bqb_opt_quick)
               bqb_optimize = 1
            CASE (bqb_opt_normal)
               bqb_optimize = 2
            CASE (bqb_opt_patient)
               bqb_optimize = 3
            CASE (bqb_opt_exhaustive)
               bqb_optimize = 4
            END SELECT
            ret = libvori_setBQBOptimization(bqb_optimize)
            ret = libvori_setBQBHistory(bqb_history)
            ret = libvori_setBQBSkipFirst(MERGE(1, 0, bqb_skip_first))
            ret = libvori_setBQBCheck(MERGE(1, 0, bqb_check))
            ret = libvori_setBQBOverwrite(MERGE(1, 0, bqb_overwrite))
            ret = libvori_setBQBStoreStep(MERGE(1, 0, bqb_store_step))
         END IF

         ret = libvori_setgrid( &
               U1 - L1 + 1, &
               U2 - L2 + 1, &
               U3 - L3 + 1, &
               rspace_pw%pw_grid%dh(1, 1)*(U1 - L1 + 1), &
               rspace_pw%pw_grid%dh(2, 1)*(U1 - L1 + 1), &
               rspace_pw%pw_grid%dh(3, 1)*(U1 - L1 + 1), &
               rspace_pw%pw_grid%dh(1, 2)*(U2 - L2 + 1), &
               rspace_pw%pw_grid%dh(2, 2)*(U2 - L2 + 1), &
               rspace_pw%pw_grid%dh(3, 2)*(U2 - L2 + 1), &
               rspace_pw%pw_grid%dh(1, 3)*(U3 - L3 + 1), &
               rspace_pw%pw_grid%dh(2, 3)*(U3 - L3 + 1), &
               rspace_pw%pw_grid%dh(3, 3)*(U3 - L3 + 1), &
               cell%hmat(1, 1), &
               cell%hmat(2, 1), &
               cell%hmat(3, 1), &
               cell%hmat(1, 2), &
               cell%hmat(2, 2), &
               cell%hmat(3, 2), &
               cell%hmat(1, 3), &
               cell%hmat(2, 3), &
               cell%hmat(3, 3) &
               )

         IF (ret /= 0) THEN
            CPABORT("The library returned an error. Aborting.")
         END IF

         ALLOCATE (particles_z(natom))
         ALLOCATE (particles_c(natom))
         ALLOCATE (particles_r(3, natom))
         ALLOCATE (particles_radius(natom))

         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff, vdw_radius=r)
            r = r*angstrom
            atomic_kind => atomic_kind_set(ikind)
            CALL get_atomic_kind(atomic_kind, atom_list=atom_list, z=ord)
            DO iat = 1, SIZE(atom_list)
               i = atom_list(iat)
               particles_c(i) = zeff
               particles_z(i) = ord
               particles_r(:, i) = particles%els(i)%r(:)
               particles_radius(i) = r
            END DO
         END DO

         ret = libvori_pushatoms(natom, particles_z, particles_c, particles_r(1, :), particles_r(2, :), particles_r(3, :))

         IF (ret /= 0) THEN
            CPABORT("The library returned an error. Aborting.")
         END IF

      END IF

      IF (iounit > 0) THEN
         IF (do_voro /= 0) THEN
            WRITE (iounit, *) "VORONOI| Collecting electron density from MPI ranks and sending to library..."
         ELSE
            WRITE (iounit, *) "BQB| Collecting electron density from MPI ranks and sending to library..."
         END IF
      END IF

      ALLOCATE (buf(L3:U3))

      dest = 0

      DO I1 = L1, U1
         DO I2 = L2, U2

            ! cycling through the CPUs, check if the current ray (I1,I2) is local to that CPU
            IF (rspace_pw%pw_grid%para%mode .NE. PW_MODE_LOCAL) THEN
               DO ip = 0, num_pe - 1
                  IF (rspace_pw%pw_grid%para%bo(1, 1, ip, 1) <= I1 - L1 + 1 &
                      .AND. rspace_pw%pw_grid%para%bo(2, 1, ip, 1) >= I1 - L1 + 1 .AND. &
                      rspace_pw%pw_grid%para%bo(1, 2, ip, 1) <= I2 - L2 + 1 &
                      .AND. rspace_pw%pw_grid%para%bo(2, 2, ip, 1) >= I2 - L2 + 1) THEN
                     source = ip
                  END IF
               END DO
            ELSE
               source = dest
            END IF

            IF (source == dest) THEN
               IF (my_rank == source) THEN
                  buf(:) = rspace_pw%cr3d(I1, I2, :)
               END IF
            ELSE
               IF (my_rank == source) THEN
                  buf(:) = rspace_pw%cr3d(I1, I2, :)
                  CALL mp_send(buf, dest, tag, gid)
               END IF
               IF (my_rank == dest) THEN
                  CALL mp_recv(buf, source, tag, gid)
               END IF
            END IF

            IF (my_rank == dest) THEN
               ret = libvori_push_rho_zrow(I1 - L1, I2 - L2, U3 - L3 + 1, buf)
               IF (ret /= 0) THEN
                  CPABORT("The library returned an error. Aborting.")
               END IF
            END IF

            ! this double loop generates so many messages that it can overload
            ! the message passing system, e.g. on XT3
            ! we therefore put a barrier here that limits the amount of message
            ! that flies around at any given time.
            ! if ever this routine becomes a bottleneck, we should go for a
            ! more complicated rewrite
            CALL mp_sync(gid)

         END DO
      END DO

      DEALLOCATE (buf)

      IF (my_rank == 0) THEN

         IF (do_voro /= 0) THEN

            IF (radius_type == voro_radii_unity) THEN
               ret = libvori_setRadii_Unity()
            ELSE IF (radius_type == voro_radii_cov) THEN
               ! Use the covalent radii from LIBVORI
               ret = libvori_setRadii_Covalent()
            ELSE IF (radius_type == voro_radii_vdw) THEN
               ! Use the van der Waals radii from CP2K
               ret = libvori_setRadii_User(100.0_dp, particles_radius)
            ELSE IF (radius_type == voro_radii_user) THEN
               ! Use the user defined atomic radii
               IF (natom /= SIZE(user_radii)) THEN
                  CALL cp_abort(__LOCATION__, &
                                "Length of keyword VORONOI\USER_RADII does not "// &
                                "match number of atoms in the input coordinate file.")
               END IF
               ret = libvori_setRadii_User(100.0_dp, user_radii)
            ELSE
               CPABORT("No valid radius type was specified for VORONOI")
            END IF

            IF (voro_sanity) THEN

               ret = libvori_sanitycheck(sim_step, sim_time)

               IF (ret /= 0) THEN
                  CPABORT("The library returned an error. Aborting.")
               END IF

            END IF

            ret = libvori_step(sim_step, sim_time)

            IF (ret /= 0) THEN
               CPABORT("The library returned an error. Aborting.")
            END IF

            ALLOCATE (voro_radii(natom))
            ALLOCATE (voro_charge(natom))
            ALLOCATE (voro_volume(natom))
            ALLOCATE (voro_dipole(natom, 3))
            ALLOCATE (voro_quadrupole(natom, 9))
            ALLOCATE (voro_buffer(natom))
            ALLOCATE (voro_wrapped_pos(natom, 3))
            ALLOCATE (voro_charge_center(natom, 3))

            ret = libvori_get_radius(natom, voro_radii)

            ret = libvori_get_charge(natom, voro_charge)

            ret = libvori_get_volume(natom, voro_volume)

            DO i1 = 1, 3
               ret = libvori_get_dipole(i1, natom, voro_buffer)
               voro_dipole(:, i1) = voro_buffer(:)
            END DO

            DO i1 = 1, 9
               ret = libvori_get_quadrupole(i1, natom, voro_buffer)
               voro_quadrupole(:, i1) = voro_buffer(:)
            END DO

            DO i1 = 1, 3
               ret = libvori_get_wrapped_pos(i1, natom, voro_buffer)
               voro_wrapped_pos(:, i1) = voro_buffer(:)
            END DO

            DO i1 = 1, 3
               ret = libvori_get_charge_center(i1, natom, voro_buffer)
               voro_charge_center(:, i1) = voro_buffer(:)
            END DO

            IF (unit_voro > 0) THEN
               WRITE (unit_voro, FMT="(I0)") natom
               WRITE (unit_voro, FMT="(A,I8,A,F12.4,A)") "# Simulation step ", sim_step, ",  Simulation time ", &
                  sim_time*femtoseconds, " fs"
               WRITE (unit_voro, FMT="(A,22A20)") "# Atom     Z", &
                  "Radius", "Position(X)", "Position(Y)", "Position(Z)", &
                  "Voronoi_Volume", "Z(eff)", "Charge", "Dipole(X)", "Dipole(Y)", "Dipole(Z)", &
                  "ChargeCenter(X)", "ChargeCenter(Y)", "ChargeCenter(Z)", &
                  "Quadrupole(XX)", "Quadrupole(XY)", "Quadrupole(XZ)", &
                  "Quadrupole(YX)", "Quadrupole(YY)", "Quadrupole(YZ)", &
                  "Quadrupole(ZX)", "Quadrupole(ZY)", "Quadrupole(ZZ)"
               DO i1 = 1, natom
                  WRITE (unit_voro, FMT="(2I6,22F20.10)") &
                     i1, &
                     particles_z(i1), &
                     voro_radii(i1)/100.0_dp, &
                     voro_wrapped_pos(i1, 1:3)/100.0_dp, &
                     voro_volume(i1)/1000000.0_dp, &
                     particles_c(i1), &
                     voro_charge(i1), &
                     voro_dipole(i1, 1:3), &
                     voro_charge_center(i1, 1:3)/100.0_dp, &
                     voro_quadrupole(i1, 1:9)
               END DO
            END IF

            DEALLOCATE (voro_radii)
            DEALLOCATE (voro_charge)
            DEALLOCATE (voro_volume)
            DEALLOCATE (voro_dipole)
            DEALLOCATE (voro_quadrupole)
            DEALLOCATE (voro_buffer)
            DEALLOCATE (voro_wrapped_pos)
            DEALLOCATE (voro_charge_center)

            IF (iounit > 0) THEN
               WRITE (iounit, *) "VORONOI| Voronoi integration finished."
            END IF

         END IF ! do_voro

         IF (do_bqb /= 0) THEN

            ret = libvori_processBQBFrame(sim_step, sim_time*femtoseconds)

            IF (ret /= 0) THEN
               CPABORT("The library returned an error. Aborting.")
            END IF

            IF (do_voro /= 0) THEN
               IF (iounit > 0) THEN
                  WRITE (iounit, *) "VORONOI| BQB compression finished."
               END IF
            ELSE
               IF (iounit > 0) THEN
                  WRITE (iounit, *) "BQB| BQB compression finished."
               END IF
            END IF

         END IF  ! do_bqb

      END IF

      IF (my_rank == 0) THEN
         DEALLOCATE (particles_z)
         DEALLOCATE (particles_c)
         DEALLOCATE (particles_r)
         DEALLOCATE (particles_radius)
      END IF

#if defined(__HAS_IEEE_EXCEPTIONS)
      CALL ieee_set_halting_mode(IEEE_ALL, halt)
#endif

      CALL timestop(handle)

#else

      MARK_USED(do_voro)
      MARK_USED(do_bqb)
      MARK_USED(input_voro)
      MARK_USED(input_bqb)
      MARK_USED(unit_voro)
      MARK_USED(qs_env)
      MARK_USED(rspace_pw)

      CALL cp_warn(__LOCATION__, &
                   "Voronoi integration and BQB output require CP2k to be compiled"// &
                   " with the -D__LIBVORI preprocessor option.")

#endif

   END SUBROUTINE entry_voronoi_or_bqb

! **************************************************************************************************
!> \brief Call libvori's finalize if support is compiled in
! **************************************************************************************************
   SUBROUTINE finalize_libvori()
#if defined(__LIBVORI)
      INTEGER(KIND=C_INT) :: ret
      ret = libvori_finalize()
#endif
   END SUBROUTINE
END MODULE voronoi_interface


