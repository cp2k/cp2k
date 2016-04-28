!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief CP2K C/C++ interface
!> \par History
!>       12.2012 created [Hossein Bani-Hashemian]
!>       04.2016 restructured [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
MODULE cp2k_c_interface
   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_DOUBLE,&
                                              C_FUNPTR,&
                                              C_INT,&
                                              C_NULL_CHAR
   USE cp2k_info,                       ONLY: cp2k_version
   USE cp2k_runs,                       ONLY: run_input
   USE f77_interface,                   ONLY: &
        calc_energy, calc_energy_force, calc_force, create_force_env, destroy_force_env, do_shake, &
        f_env_add_defaults, f_env_rm_defaults, f_env_type, finalize_cp2k, get_energy, get_force, &
        get_natom, get_nparticle, get_pos, get_result_r1, init_cp2k, set_pos, set_vel
   USE force_env_types,                 ONLY: force_env_get,&
                                              use_qs_force
   USE input_cp2k,                      ONLY: create_cp2k_root_section
   USE input_section_types,             ONLY: section_release,&
                                              section_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

CONTAINS

! **************************************************************************************************
!> \brief Copy the content of a \0-terminated C-string to a finite-length Fortran string
!>
!> The content of the new string may be truncated if the number of characters before the '\0'
!> in the source string exceed the length of the destination string.
!> \param fstring destination string
!> \param cstring source string
!> \author Tiziano Mueller
! **************************************************************************************************
   SUBROUTINE strncpy_c2f(fstring, cstring)
      CHARACTER(LEN=*), INTENT(OUT)                      :: fstring
      CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN)          :: cstring(*)

      INTEGER                                            :: i

      DO i = 1, LEN(fstring)
         IF (cstring(i) == C_NULL_CHAR) EXIT
         fstring(i:i) = cstring(i)
      END DO
   END SUBROUTINE strncpy_c2f

! **************************************************************************************************
!> \brief ...
!> \param version_str ...
!> \param str_length ...
! **************************************************************************************************
   SUBROUTINE cp_c_get_version(version_str, str_length) BIND(C)
      CHARACTER(LEN=1, KIND=C_CHAR), INTENT(INOUT)       :: version_str(*)
      INTEGER(C_INT), VALUE                              :: str_length

      INTEGER                                            :: i, n

      n = LEN_TRIM(cp2k_version)
      CPASSERT(str_length > n+1)

      ! copy string
      DO i = 1, n
         version_str(i) = cp2k_version(i:i)
      ENDDO
      version_str(n+1) = C_NULL_CHAR
   END SUBROUTINE cp_c_get_version

! **************************************************************************************************
!> \brief ...
!> \param init_mpi ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE cp_c_init_cp2k(init_mpi, ierr) BIND(C)
      INTEGER(C_INT)                                     :: init_mpi, ierr

      LOGICAL                                            :: i_mpi

      i_mpi = init_mpi /= 0
      CALL init_cp2k(i_mpi, ierr)
   END SUBROUTINE cp_c_init_cp2k

! **************************************************************************************************
!> \brief ...
!> \param finalize_mpi ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE cp_c_finalize_cp2k(finalize_mpi, ierr) BIND(C)
      INTEGER(C_INT)                                     :: finalize_mpi, ierr

      LOGICAL                                            :: f_mpi

      f_mpi = finalize_mpi /= 0
      CALL finalize_cp2k(f_mpi, ierr)
   END SUBROUTINE cp_c_finalize_cp2k

! **************************************************************************************************
!> \brief ...
!> \param new_env_id ...
!> \param input_file_path ...
!> \param output_file_path ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE cp_c_create_fenv(new_env_id, input_file_path, output_file_path, ierr) BIND(C)
      INTEGER(C_INT)                                     :: new_env_id
      CHARACTER(LEN=1, KIND=C_CHAR)                      :: input_file_path(*), output_file_path(*)
      INTEGER(C_INT)                                     :: ierr

      CHARACTER(LEN=default_path_length)                 :: ifp, ofp
      TYPE(section_type), POINTER                        :: input_declaration

      ifp = " "; ofp = " "
      CALL strncpy_c2f(ifp, input_file_path)
      CALL strncpy_c2f(ofp, output_file_path)

      NULLIFY (input_declaration)
      CALL create_cp2k_root_section(input_declaration)
      CALL create_force_env(new_env_id, input_declaration, ifp, ofp, ierr=ierr)
      CALL section_release(input_declaration)
   END SUBROUTINE cp_c_create_fenv

! **************************************************************************************************
!> \brief ...
!> \param new_env_id ...
!> \param input_file_path ...
!> \param output_file_path ...
!> \param mpi_comm ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE cp_c_create_fenv_comm(new_env_id, input_file_path, output_file_path, mpi_comm, ierr) BIND(C)
      INTEGER(C_INT)                                     :: new_env_id
      CHARACTER(LEN=1, KIND=C_CHAR)                      :: input_file_path(*), output_file_path(*)
      INTEGER(C_INT)                                     :: mpi_comm, ierr

      CHARACTER(LEN=default_path_length)                 :: ifp, ofp
      TYPE(section_type), POINTER                        :: input_declaration

      ifp = " "; ofp = " "
      CALL strncpy_c2f(ifp, input_file_path)
      CALL strncpy_c2f(ofp, output_file_path)

      NULLIFY (input_declaration)
      CALL create_cp2k_root_section(input_declaration)
      CALL create_force_env(new_env_id, input_declaration, ifp, ofp, mpi_comm, ierr=ierr)
      CALL section_release(input_declaration)
   END SUBROUTINE cp_c_create_fenv_comm

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param ierr ...
! **************************************************************************************************
   RECURSIVE SUBROUTINE cp_c_destroy_fenv(env_id, ierr) BIND(C)
      INTEGER(C_INT)                                     :: env_id, ierr

      CALL destroy_force_env(env_id, ierr)
   END SUBROUTINE cp_c_destroy_fenv

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param new_pos ...
!> \param n_el ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE cp_c_set_pos(env_id, new_pos, n_el, ierr) BIND(C)
      INTEGER(C_INT)                                     :: env_id, n_el
      REAL(C_DOUBLE), DIMENSION(1:n_el)                  :: new_pos
      INTEGER(C_INT)                                     :: ierr

      CALL set_pos(env_id, new_pos, n_el, ierr)
   END SUBROUTINE cp_c_set_pos

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param new_vel ...
!> \param n_el ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE cp_c_set_vel(env_id, new_vel, n_el, ierr) BIND(C)
      INTEGER(C_INT)                                     :: env_id, n_el
      REAL(C_DOUBLE), DIMENSION(1:n_el)                  :: new_vel
      INTEGER(C_INT)                                     :: ierr

      CALL set_vel(env_id, new_vel, n_el, ierr)
   END SUBROUTINE cp_c_set_vel

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param description ...
!> \param N ...
!> \param RESULT ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE cp_c_get_result_r1(env_id, description, N, RESULT, ierr) BIND(C)
      INTEGER(C_INT)                                     :: env_id
      CHARACTER(LEN=1, KIND=C_CHAR)                      :: description(*)
      INTEGER(C_INT)                                     :: N
      REAL(C_DOUBLE), DIMENSION(1:N)                     :: RESULT
      INTEGER(C_INT)                                     :: ierr

      CHARACTER(LEN=default_string_length)               :: desc_low

      desc_low = " "
      CALL strncpy_c2f(desc_low, description)

      CALL get_result_r1(env_id, desc_low, N, RESULT, ierr=ierr)
   END SUBROUTINE cp_c_get_result_r1

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param natom ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE cp_c_get_natom(env_id, natom, ierr) BIND(C)
      INTEGER(C_INT)                                     :: env_id, natom, ierr

      CALL get_natom(env_id, natom, ierr)
   END SUBROUTINE cp_c_get_natom

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param nparticle ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE cp_c_get_nparticle(env_id, nparticle, ierr) BIND(C)
      INTEGER(C_INT)                                     :: env_id, nparticle, ierr

      CALL get_nparticle(env_id, nparticle, ierr)
   END SUBROUTINE cp_c_get_nparticle

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param pos ...
!> \param n_el ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE cp_c_get_pos(env_id, pos, n_el, ierr) BIND(C)
      INTEGER(C_INT)                                     :: env_id, n_el
      REAL(C_DOUBLE), DIMENSION(1:n_el)                  :: pos
      INTEGER(C_INT)                                     :: ierr

      CALL get_pos(env_id, pos, n_el, ierr)
   END SUBROUTINE cp_c_get_pos

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param force ...
!> \param n_el ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE cp_c_get_force(env_id, force, n_el, ierr) BIND(C)
      INTEGER(C_INT)                                     :: env_id, n_el
      REAL(C_DOUBLE), DIMENSION(1:n_el)                  :: force
      INTEGER(C_INT)                                     :: ierr

      CALL get_force(env_id, force, n_el, ierr)
   END SUBROUTINE cp_c_get_force

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param e_pot ...
!> \param ierr ...
! **************************************************************************************************
   RECURSIVE SUBROUTINE cp_c_get_energy(env_id, e_pot, ierr) BIND(C)
      INTEGER(C_INT)                                     :: env_id
      REAL(C_DOUBLE)                                     :: e_pot
      INTEGER(C_INT)                                     :: ierr

      CALL get_energy(env_id, e_pot, ierr)
   END SUBROUTINE cp_c_get_energy

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param calc_force ...
!> \param ierr ...
! **************************************************************************************************
   RECURSIVE SUBROUTINE cp_c_calc_energy_force(env_id, calc_force, ierr) BIND(C)
      INTEGER(C_INT)                                     :: env_id, calc_force, ierr

      LOGICAL                                            :: cf

      cf = (calc_force /= 0)
      CALL calc_energy_force(env_id, cf, ierr)
   END SUBROUTINE cp_c_calc_energy_force

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param pos ...
!> \param n_el ...
!> \param e_pot ...
!> \param ierr ...
! **************************************************************************************************
   RECURSIVE SUBROUTINE cp_c_calc_energy(env_id, pos, n_el, e_pot, ierr) BIND(C)
      INTEGER(C_INT)                                     :: env_id, n_el
      REAL(C_DOUBLE), DIMENSION(1:n_el)                  :: pos
      REAL(C_DOUBLE)                                     :: e_pot
      INTEGER(C_INT)                                     :: ierr

      CALL calc_energy(env_id, pos, n_el, e_pot, ierr)
   END SUBROUTINE cp_c_calc_energy

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param pos ...
!> \param n_el_pos ...
!> \param e_pot ...
!> \param force ...
!> \param n_el_force ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE cp_c_calc_force(env_id, pos, n_el_pos, e_pot, force, n_el_force, ierr) BIND(C)
      INTEGER(C_INT)                                     :: env_id, n_el_pos
      REAL(C_DOUBLE), DIMENSION(1:n_el_pos)              :: pos
      REAL(C_DOUBLE)                                     :: e_pot
      INTEGER(C_INT)                                     :: n_el_force
      REAL(C_DOUBLE), DIMENSION(1:n_el_force)            :: force
      INTEGER(C_INT)                                     :: ierr

      CALL calc_force(env_id, pos, n_el_pos, e_pot, force, n_el_force, ierr)
   END SUBROUTINE cp_c_calc_force

! **************************************************************************************************
!> \brief ...
!> \param input_file_path ...
!> \param output_file_path ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE cp_c_run_input(input_file_path, output_file_path, ierr) BIND(C)
      CHARACTER(LEN=1, KIND=C_CHAR)                      :: input_file_path(*), output_file_path(*)
      INTEGER(C_INT)                                     :: ierr

      CHARACTER(LEN=default_path_length)                 :: ifp, ofp
      TYPE(section_type), POINTER                        :: input_declaration

      ifp = " "; ofp = " "
      CALL strncpy_c2f(ifp, input_file_path)
      CALL strncpy_c2f(ofp, output_file_path)

      NULLIFY (input_declaration)
      CALL create_cp2k_root_section(input_declaration)
      CALL run_input(input_declaration, ifp, ofp, ierr=ierr)
      CALL section_release(input_declaration)
   END SUBROUTINE cp_c_run_input

! **************************************************************************************************
!> \brief ...
!> \param input_file_path ...
!> \param output_file_path ...
!> \param mpi_comm ...
!> \param ierr ...
! **************************************************************************************************
   RECURSIVE SUBROUTINE cp_c_run_input_comm(input_file_path, output_file_path, mpi_comm, ierr) BIND(C)
      CHARACTER(LEN=1, KIND=C_CHAR)                      :: input_file_path(*), output_file_path(*)
      INTEGER(C_INT)                                     :: mpi_comm, ierr

      CHARACTER(LEN=default_path_length)                 :: ifp, ofp
      TYPE(section_type), POINTER                        :: input_declaration

      ifp = " "; ofp = " "
      CALL strncpy_c2f(ifp, input_file_path)
      CALL strncpy_c2f(ofp, output_file_path)

      NULLIFY (input_declaration)
      CALL create_cp2k_root_section(input_declaration)
      CALL run_input(input_declaration, ifp, ofp, ierr, mpi_comm)
      CALL section_release(input_declaration)
   END SUBROUTINE cp_c_run_input_comm

! **************************************************************************************************
!> \brief ...
!> \param f_env_id ...
!> \param dt ...
!> \param shake_tol ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE cp_c_do_shake(f_env_id, dt, shake_tol, ierr) BIND(C)
      INTEGER(C_INT)                                     :: f_env_id
      REAL(C_DOUBLE)                                     :: dt, shake_tol
      INTEGER(C_INT)                                     :: ierr

      CALL do_shake(f_env_id, dt, shake_tol, ierr)
   END SUBROUTINE cp_c_do_shake

! **************************************************************************************************
!> \brief Gets a function pointer pointing to a routine defined in C/C++ and
!>        passes it to the transport environment in force environment
!> \param f_env_id  the force env id
!> \param func      the function pointer
!> \param ierr      will return a number different from 0 if there is an error
!> \par History
!>      12.2012 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE cp_c_ext_method_set_ptr(f_env_id, func, ierr) BIND(C)
      INTEGER(C_INT)                                     :: f_env_id
      TYPE(C_FUNPTR), VALUE                              :: func
      INTEGER(C_INT)                                     :: ierr

      INTEGER                                            :: in_use
      TYPE(f_env_type), POINTER                          :: f_env

      NULLIFY (f_env)
      CALL f_env_add_defaults(f_env_id, f_env)
      CALL force_env_get(f_env%force_env, in_use=in_use)
      IF (in_use .EQ. use_qs_force) THEN
         f_env%force_env%qs_env%transport_env%ext_c_method_ptr = func
      END IF
      CALL f_env_rm_defaults(f_env, ierr)
   END SUBROUTINE cp_c_ext_method_set_ptr

END MODULE cp2k_c_interface

