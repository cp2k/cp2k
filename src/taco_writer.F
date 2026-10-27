!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Bridge between CP2K trajectory output and the TACO binary trajectory format.
!> \par History
!>      * 2026-10 Nils Holle: created
! **************************************************************************************************
MODULE taco_writer

#ifdef __TACO
   USE kinds, ONLY: dp
   USE particle_types, ONLY: particle_type
   USE cell_types, ONLY: cell_type
   USE ISO_C_BINDING, ONLY: c_associated, c_char, c_double, c_float, c_int, c_int32_t, &
                            c_int64_t, c_loc, c_null_char, c_null_ptr, c_ptr
   USE physcon, ONLY: angstrom, evolt, femtoseconds
#endif
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: taco_close_file
#ifdef __TACO
   PUBLIC :: taco_write_frame

   INTEGER, PARAMETER, PRIVATE :: symbol_width = 4

   TYPE(c_ptr), SAVE           :: taco_handle = c_null_ptr
   INTEGER, SAVE               :: taco_natoms = -1
   LOGICAL, SAVE               :: taco_open = .FALSE.

   INTERFACE
      FUNCTION taco_cp2k_create(path, natoms, time_step, software, temperature, pressure, &
                                timestep_fs, symbols, masses) RESULT(w) &
         BIND(C, name="taco_cp2k_create")
         IMPORT :: c_char, c_double, c_int32_t, c_ptr
         TYPE(c_ptr), VALUE               :: path
         INTEGER(c_int32_t), VALUE        :: natoms
         REAL(c_double), VALUE            :: time_step
         TYPE(c_ptr), VALUE               :: software
         REAL(c_double), VALUE            :: temperature, pressure, timestep_fs
         TYPE(c_ptr), VALUE               :: symbols, masses
         TYPE(c_ptr)                      :: w
      END FUNCTION taco_cp2k_create

      FUNCTION taco_cp2k_write_frame(writer, frame_number, time, natoms, pos, vel, force, box, &
                                     potential_energy) RESULT(ierr) &
         BIND(C, name="taco_cp2k_write_frame")
         IMPORT :: c_double, c_float, c_int, c_int32_t, c_int64_t, c_ptr
         TYPE(c_ptr), VALUE               :: writer
         INTEGER(c_int64_t), VALUE        :: frame_number
         REAL(c_double), VALUE            :: time
         INTEGER(c_int32_t), VALUE        :: natoms
         TYPE(c_ptr), VALUE               :: pos, vel, force, box
         REAL(c_double), VALUE            :: potential_energy
         INTEGER(c_int)                   :: ierr
      END FUNCTION taco_cp2k_write_frame

      FUNCTION taco_cp2k_finish(writer) RESULT(ierr) BIND(C, name="taco_cp2k_finish")
         IMPORT :: c_int, c_ptr
         TYPE(c_ptr), VALUE               :: writer
         INTEGER(c_int)                   :: ierr
      END FUNCTION taco_cp2k_finish

      SUBROUTINE taco_cp2k_free(writer) BIND(C, name="taco_cp2k_free")
         IMPORT :: c_ptr
         TYPE(c_ptr), VALUE               :: writer
      END SUBROUTINE taco_cp2k_free
   END INTERFACE
#endif

CONTAINS

! **************************************************************************************************
!> \brief Write one frame to the TACO trajectory file, creating it lazily on the first call.
!> \param filename   Path of the TACO file.
!> \param particle_set CP2K particle set (positions, velocities, forces, atomic kinds).
!> \param cell       Simulation cell.
!> \param it         Step number.
!> \param time       Time (CP2K internal units).
!> \param dtime      Time step (CP2K internal units).
!> \param etot       Total energy (CP2K internal units).
! **************************************************************************************************
#ifdef __TACO
   SUBROUTINE taco_write_frame(filename, particle_set, cell, it, time, dtime, etot)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: it
      REAL(KIND=dp), INTENT(IN)                          :: time, dtime, etot

      CALL taco_write_frame_impl(filename, particle_set, cell, it, time, dtime, etot)
   END SUBROUTINE taco_write_frame
#endif

! **************************************************************************************************
!> \brief Finalise and free the currently open TACO file.
! **************************************************************************************************
   SUBROUTINE taco_close_file()
#ifdef __TACO
      CALL taco_close_file_impl()
#endif
   END SUBROUTINE taco_close_file

#ifdef __TACO
! **************************************************************************************************
!> \brief Implementation of the TACO frame writer (only compiled with __TACO).
!> \param filename ...
!> \param particle_set ...
!> \param cell ...
!> \param it ...
!> \param time ...
!> \param dtime ...
!> \param etot ...
! **************************************************************************************************
   SUBROUTINE taco_write_frame_impl(filename, particle_set, cell, it, time, dtime, etot)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: it
      REAL(KIND=dp), INTENT(IN)                          :: time, dtime, etot

      CHARACTER(LEN=*), PARAMETER :: routineN = 'taco_write_frame_impl'

      INTEGER                                            :: handle, i, k, nat
      INTEGER(c_int)                                     :: ierr
      REAL(c_double)                                     :: tstep_ps
      REAL(c_float), ALLOCATABLE, TARGET                 :: force_c(:), pos_c(:), vel_c(:)
      REAL(c_float), TARGET                              :: box_c(3)

      CALL timeset(routineN, handle)

      nat = SIZE(particle_set)
      IF (.NOT. taco_open) THEN
         CALL taco_open_file(filename, particle_set, nat, dtime)
      END IF
      CPASSERT(taco_open)
      CPASSERT(nat == taco_natoms)

      ALLOCATE (pos_c(3*nat), vel_c(3*nat), force_c(3*nat))
      DO i = 1, nat
         k = 3*(i - 1)
         ! positions: [Bohr] -> [Angstrom]
         pos_c(k + 1) = REAL(particle_set(i)%r(1)*angstrom, KIND=c_float)
         pos_c(k + 2) = REAL(particle_set(i)%r(2)*angstrom, KIND=c_float)
         pos_c(k + 3) = REAL(particle_set(i)%r(3)*angstrom, KIND=c_float)
         ! velocities: [Bohr/a.u.] -> [Angstrom/ps]
         vel_c(k + 1) = REAL(particle_set(i)%v(1)*angstrom*1000.0_dp/femtoseconds, KIND=c_float)
         vel_c(k + 2) = REAL(particle_set(i)%v(2)*angstrom*1000.0_dp/femtoseconds, KIND=c_float)
         vel_c(k + 3) = REAL(particle_set(i)%v(3)*angstrom*1000.0_dp/femtoseconds, KIND=c_float)
         ! forces: [Hartree/Bohr] -> [eV/Angstrom]
         force_c(k + 1) = REAL(particle_set(i)%f(1)*evolt/angstrom, KIND=c_float)
         force_c(k + 2) = REAL(particle_set(i)%f(2)*evolt/angstrom, KIND=c_float)
         force_c(k + 3) = REAL(particle_set(i)%f(3)*evolt/angstrom, KIND=c_float)
      END DO

      ! box dimensions: lattice vector lengths [Bohr] -> [Angstrom]
      box_c(1) = REAL(SQRT(SUM(cell%hmat(:, 1)**2))*angstrom, KIND=c_float)
      box_c(2) = REAL(SQRT(SUM(cell%hmat(:, 2)**2))*angstrom, KIND=c_float)
      box_c(3) = REAL(SQRT(SUM(cell%hmat(:, 3)**2))*angstrom, KIND=c_float)

      tstep_ps = REAL(dtime*femtoseconds/1000.0_dp, KIND=c_double)

      ierr = taco_cp2k_write_frame(taco_handle, INT(it, c_int64_t), &
                                   REAL(time*femtoseconds/1000.0_dp, c_double), &
                                   INT(nat, c_int32_t), c_loc(pos_c(1)), c_loc(vel_c(1)), &
                                   c_loc(force_c(1)), c_loc(box_c(1)), &
                                   REAL(etot*evolt, c_double))
      IF (ierr /= 0) THEN
         CALL cp_warn(__LOCATION__, "Error writing TACO trajectory frame")
      END IF

      DEALLOCATE (pos_c, vel_c, force_c)

      CALL timestop(handle)

   END SUBROUTINE taco_write_frame_impl

! **************************************************************************************************
!> \brief Open (or reopen) the TACO file, writing the header with atom metadata.
!> \param filename ...
!> \param particle_set ...
!> \param nat ...
!> \param dtime ...
! **************************************************************************************************
   SUBROUTINE taco_open_file(filename, particle_set, nat, dtime)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, INTENT(IN)                                :: nat
      REAL(KIND=dp), INTENT(IN)                          :: dtime

      CHARACTER(KIND=c_char), ALLOCATABLE, TARGET        :: path_c(:), symbols(:)
      CHARACTER(KIND=c_char), TARGET                     :: software_c(5)
      INTEGER                                            :: i, k
      REAL(c_double)                                     :: timestep_fs, tstep_ps
      REAL(c_double), ALLOCATABLE, TARGET                :: masses(:)

      IF (taco_open) CALL taco_close_file_impl()

      ALLOCATE (path_c(LEN_TRIM(filename) + 1))
      DO k = 1, LEN_TRIM(filename)
         path_c(k) = filename(k:k)
      END DO
      path_c(LEN_TRIM(filename) + 1) = c_null_char

      software_c(1:4) = ['C', 'P', '2', 'K']
      software_c(5) = c_null_char

      ALLOCATE (masses(nat))
      ALLOCATE (symbols(nat*symbol_width))
      symbols = ' '
      DO i = 1, nat
         masses(i) = REAL(particle_set(i)%atomic_kind%mass, KIND=c_double)
         k = (i - 1)*symbol_width
         CALL taco_copy_symbol(particle_set(i)%atomic_kind%element_symbol, symbols(k + 1:k + symbol_width))
      END DO

      timestep_fs = REAL(dtime*femtoseconds, KIND=c_double)
      tstep_ps = REAL(dtime*femtoseconds/1000.0_dp, KIND=c_double)

      taco_handle = taco_cp2k_create(c_loc(path_c(1)), INT(nat, c_int32_t), tstep_ps, &
                                     c_loc(software_c(1)), 0.0_c_double, 0.0_c_double, &
                                     timestep_fs, c_loc(symbols(1)), c_loc(masses(1)))
      IF (.NOT. C_ASSOCIATED(taco_handle)) THEN
         CALL cp_abort(__LOCATION__, "Could not open TACO file "//TRIM(filename))
      END IF
      taco_open = .TRUE.
      taco_natoms = nat

      DEALLOCATE (path_c, masses, symbols)

   END SUBROUTINE taco_open_file

! **************************************************************************************************
!> \brief Copy a (short) element symbol into a fixed-width space-padded buffer.
!> \param symbol ...
!> \param buf ...
! **************************************************************************************************
   SUBROUTINE taco_copy_symbol(symbol, buf)
      CHARACTER(LEN=*), INTENT(IN)                       :: symbol
      CHARACTER(KIND=c_char), INTENT(OUT)                :: buf(:)

      INTEGER                                            :: i, n

      buf = ' '
      n = MIN(LEN_TRIM(symbol), SIZE(buf))
      DO i = 1, n
         buf(i) = symbol(i:i)
      END DO

   END SUBROUTINE taco_copy_symbol

! **************************************************************************************************
!> \brief Implementation of the close/free routine (only compiled with __TACO).
! **************************************************************************************************
   SUBROUTINE taco_close_file_impl()
      CHARACTER(LEN=*), PARAMETER :: routineN = 'taco_close_file_impl'

      INTEGER                                            :: handle
      INTEGER(c_int)                                     :: ierr

      CALL timeset(routineN, handle)

      IF (taco_open) THEN
         ierr = taco_cp2k_finish(taco_handle)
         CALL taco_cp2k_free(taco_handle)
         taco_open = .FALSE.
         taco_handle = c_null_ptr
         taco_natoms = -1
      END IF

      CALL timestop(handle)

   END SUBROUTINE taco_close_file_impl

#endif

END MODULE taco_writer
