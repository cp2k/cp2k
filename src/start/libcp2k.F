!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief CP2K C/C++ interface
!> \par History
!>       12.2012 created [Hossein Bani-Hashemian]
!>       04.2016 restructured [Hossein Bani-Hashemian, Ole Schuett]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
MODULE libcp2k
   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_DOUBLE,&
                                              C_FUNPTR,&
                                              C_INT,&
                                              C_NULL_CHAR
   USE cp2k_info,                       ONLY: cp2k_version
   USE cp2k_runs,                       ONLY: run_input
   USE f77_interface,                   ONLY: &
        calc_energy_force, create_force_env, destroy_force_env, f_env_add_defaults, &
        f_env_rm_defaults, f_env_type, finalize_cp2k, get_energy, get_force, get_natom, &
        get_nparticle, get_pos, get_result_r1, init_cp2k, set_pos, set_vel
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
!> \brief ...
!> \param version_str ...
!> \param str_length ...
! **************************************************************************************************
   SUBROUTINE cp2k_get_version(version_str, str_length) BIND(C)
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
   END SUBROUTINE cp2k_get_version

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE cp2k_init() BIND(C)
      INTEGER                                            :: ierr

      CALL init_cp2k(.TRUE., ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_init

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE cp2k_init_without_mpi() BIND(C)
      INTEGER                                            :: ierr

      CALL init_cp2k(.FALSE., ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_init_without_mpi

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE cp2k_finalize() BIND(C)
      INTEGER                                            :: ierr

      CALL finalize_cp2k(.TRUE., ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_finalize

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE cp2k_finalize_without_mpi() BIND(C)
      INTEGER                                            :: ierr

      CALL finalize_cp2k(.FALSE., ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_finalize_without_mpi

! **************************************************************************************************
!> \brief ...
!> \param new_env_id ...
!> \param input_file_path ...
!> \param output_file_path ...
! **************************************************************************************************
   SUBROUTINE cp2k_create_force_env(new_env_id, input_file_path, output_file_path) BIND(C)
      INTEGER(C_INT)                                     :: new_env_id
      CHARACTER(LEN=1, KIND=C_CHAR)                      :: input_file_path(*), output_file_path(*)

      CHARACTER(LEN=default_path_length)                 :: ifp, ofp
      INTEGER                                            :: ierr
      TYPE(section_type), POINTER                        :: input_declaration

      ifp = " "; ofp = " "
      CALL strncpy_c2f(ifp, input_file_path)
      CALL strncpy_c2f(ofp, output_file_path)

      NULLIFY (input_declaration)
      CALL create_cp2k_root_section(input_declaration)
      CALL create_force_env(new_env_id, input_declaration, ifp, ofp, ierr=ierr)
      CALL section_release(input_declaration)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_create_force_env

! **************************************************************************************************
!> \brief ...
!> \param new_env_id ...
!> \param input_file_path ...
!> \param output_file_path ...
!> \param mpi_comm ...
! **************************************************************************************************
   SUBROUTINE cp2k_create_force_env_comm(new_env_id, input_file_path, output_file_path, mpi_comm) BIND(C)
      INTEGER(C_INT)                                     :: new_env_id
      CHARACTER(LEN=1, KIND=C_CHAR)                      :: input_file_path(*), output_file_path(*)
      INTEGER(C_INT), VALUE                              :: mpi_comm

      CHARACTER(LEN=default_path_length)                 :: ifp, ofp
      INTEGER                                            :: ierr
      TYPE(section_type), POINTER                        :: input_declaration

      ifp = " "; ofp = " "
      CALL strncpy_c2f(ifp, input_file_path)
      CALL strncpy_c2f(ofp, output_file_path)

      NULLIFY (input_declaration)
      CALL create_cp2k_root_section(input_declaration)
      CALL create_force_env(new_env_id, input_declaration, ifp, ofp, mpi_comm, ierr=ierr)
      CALL section_release(input_declaration)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_create_force_env_comm

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
! **************************************************************************************************
   SUBROUTINE cp2k_destroy_force_env(env_id) BIND(C)
      INTEGER(C_INT), VALUE                              :: env_id

      INTEGER                                            :: ierr

      CALL destroy_force_env(env_id, ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_destroy_force_env

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param new_pos ...
!> \param n_el ...
! **************************************************************************************************
   SUBROUTINE cp2k_set_positions(env_id, new_pos, n_el) BIND(C)
      INTEGER(C_INT), VALUE                              :: env_id, n_el
      REAL(C_DOUBLE), DIMENSION(1:n_el)                  :: new_pos

      INTEGER                                            :: ierr

      CALL set_pos(env_id, new_pos, n_el, ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_set_positions

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param new_vel ...
!> \param n_el ...
! **************************************************************************************************
   SUBROUTINE cp2k_set_velocities(env_id, new_vel, n_el) BIND(C)
      INTEGER(C_INT), VALUE                              :: env_id, n_el
      REAL(C_DOUBLE), DIMENSION(1:n_el)                  :: new_vel

      INTEGER                                            :: ierr

      CALL set_vel(env_id, new_vel, n_el, ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_set_velocities

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param description ...
!> \param RESULT ...
!> \param n_el ...
! **************************************************************************************************
   SUBROUTINE cp2k_get_result(env_id, description, RESULT, n_el) BIND(C)
      INTEGER(C_INT), VALUE                              :: env_id
      CHARACTER(LEN=1, KIND=C_CHAR)                      :: description(*)
      INTEGER(C_INT), VALUE                              :: n_el
      REAL(C_DOUBLE), DIMENSION(1:n_el)                  :: RESULT

      CHARACTER(LEN=default_string_length)               :: desc_low
      INTEGER                                            :: ierr

      desc_low = " "
      CALL strncpy_c2f(desc_low, description)

      CALL get_result_r1(env_id, desc_low, n_el, RESULT, ierr=ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_get_result

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param natom ...
! **************************************************************************************************
   SUBROUTINE cp2k_get_natom(env_id, natom) BIND(C)
      INTEGER(C_INT), VALUE                              :: env_id
      INTEGER(C_INT)                                     :: natom

      INTEGER                                            :: ierr

      CALL get_natom(env_id, natom, ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_get_natom

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param nparticle ...
! **************************************************************************************************
   SUBROUTINE cp2k_get_nparticle(env_id, nparticle) BIND(C)
      INTEGER(C_INT), VALUE                              :: env_id
      INTEGER(C_INT)                                     :: nparticle

      INTEGER                                            :: ierr

      CALL get_nparticle(env_id, nparticle, ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_get_nparticle

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param pos ...
!> \param n_el ...
! **************************************************************************************************
   SUBROUTINE cp2k_get_positions(env_id, pos, n_el) BIND(C)
      INTEGER(C_INT), VALUE                              :: env_id, n_el
      REAL(C_DOUBLE), DIMENSION(1:n_el)                  :: pos

      INTEGER                                            :: ierr

      CALL get_pos(env_id, pos, n_el, ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_get_positions

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param force ...
!> \param n_el ...
! **************************************************************************************************
   SUBROUTINE cp2k_get_forces(env_id, force, n_el) BIND(C)
      INTEGER(C_INT), VALUE                              :: env_id, n_el
      REAL(C_DOUBLE), DIMENSION(1:n_el)                  :: force

      INTEGER                                            :: ierr

      CALL get_force(env_id, force, n_el, ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_get_forces

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param e_pot ...
! **************************************************************************************************
   SUBROUTINE cp2k_get_potential_energy(env_id, e_pot) BIND(C)
      INTEGER(C_INT), VALUE                              :: env_id
      REAL(C_DOUBLE)                                     :: e_pot

      INTEGER                                            :: ierr

      CALL get_energy(env_id, e_pot, ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_get_potential_energy

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
! **************************************************************************************************
   SUBROUTINE cp2k_calc_energy_force(env_id) BIND(C)
      INTEGER(C_INT), VALUE                              :: env_id

      INTEGER                                            :: ierr

      CALL calc_energy_force(env_id, .TRUE., ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_calc_energy_force

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
! **************************************************************************************************
   SUBROUTINE cp2k_calc_energy(env_id) BIND(C)
      INTEGER(C_INT), VALUE                              :: env_id

      INTEGER                                            :: ierr

      CALL calc_energy_force(env_id, .FALSE., ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_calc_energy

! **************************************************************************************************
!> \brief ...
!> \param input_file_path ...
!> \param output_file_path ...
! **************************************************************************************************
   SUBROUTINE cp2k_run_input(input_file_path, output_file_path) BIND(C)
      CHARACTER(LEN=1, KIND=C_CHAR)                      :: input_file_path(*), output_file_path(*)

      CHARACTER(LEN=default_path_length)                 :: ifp, ofp
      INTEGER                                            :: ierr
      TYPE(section_type), POINTER                        :: input_declaration

      ifp = " "; ofp = " "
      CALL strncpy_c2f(ifp, input_file_path)
      CALL strncpy_c2f(ofp, output_file_path)

      NULLIFY (input_declaration)
      CALL create_cp2k_root_section(input_declaration)
      CALL run_input(input_declaration, ifp, ofp, ierr=ierr)
      CALL section_release(input_declaration)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_run_input

! **************************************************************************************************
!> \brief ...
!> \param input_file_path ...
!> \param output_file_path ...
!> \param mpi_comm ...
! **************************************************************************************************
   SUBROUTINE cp2k_run_input_comm(input_file_path, output_file_path, mpi_comm) BIND(C)
      CHARACTER(LEN=1, KIND=C_CHAR)                      :: input_file_path(*), output_file_path(*)
      INTEGER(C_INT), VALUE                              :: mpi_comm

      CHARACTER(LEN=default_path_length)                 :: ifp, ofp
      INTEGER                                            :: ierr
      TYPE(section_type), POINTER                        :: input_declaration

      ifp = " "; ofp = " "
      CALL strncpy_c2f(ifp, input_file_path)
      CALL strncpy_c2f(ofp, output_file_path)

      NULLIFY (input_declaration)
      CALL create_cp2k_root_section(input_declaration)
      CALL run_input(input_declaration, ifp, ofp, ierr, mpi_comm)
      CALL section_release(input_declaration)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_run_input_comm

! **************************************************************************************************
!> \brief Gets a function pointer pointing to a routine defined in C/C++ and
!>        passes it to the transport environment in force environment
!> \param f_env_id  the force env id
!> \param func_ptr the function pointer
!> \par History
!>      12.2012 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE cp2k_transport_set_callback(f_env_id, func_ptr) BIND(C)
      INTEGER(C_INT), VALUE                              :: f_env_id
      TYPE(C_FUNPTR), VALUE                              :: func_ptr

      INTEGER                                            :: ierr, in_use
      TYPE(f_env_type), POINTER                          :: f_env

      NULLIFY (f_env)
      CALL f_env_add_defaults(f_env_id, f_env)
      CALL force_env_get(f_env%force_env, in_use=in_use)
      IF (in_use .EQ. use_qs_force) THEN
         f_env%force_env%qs_env%transport_env%ext_c_method_ptr = func_ptr
      END IF
      CALL f_env_rm_defaults(f_env, ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_transport_set_callback

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

END MODULE libcp2k
