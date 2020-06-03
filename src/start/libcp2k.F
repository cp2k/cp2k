!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! IMPORTANT: Update libcp2k.h when you add, remove or change a function in this file.              !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief CP2K C/C++ interface
!> \par History
!>       12.2012 created [Hossein Bani-Hashemian]
!>       04.2016 restructured [Hossein Bani-Hashemian, Ole Schuett]
!>       03.2018 added Active Space functions [Tiziano Mueller]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
MODULE libcp2k
   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_DOUBLE,&
                                              C_FUNPTR,&
                                              C_INT,&
                                              C_LONG,&
                                              C_NULL_CHAR
   USE cp2k_info,                       ONLY: cp2k_version
   USE cp2k_runs,                       ONLY: run_input
   USE cp_fm_types,                     ONLY: cp_fm_get_element
   USE f77_interface,                   ONLY: &
        calc_energy_force, create_force_env, destroy_force_env, f_env_add_defaults, &
        f_env_rm_defaults, f_env_type, finalize_cp2k, get_energy, get_force, get_natom, &
        get_nparticle, get_pos, get_result_r1, init_cp2k, set_cell, set_pos, set_vel
   USE force_env_types,                 ONLY: force_env_get,&
                                              use_qs_force
   USE input_cp2k,                      ONLY: create_cp2k_root_section
   USE input_cp2k_read,                 ONLY: empty_initial_variables
   USE input_section_types,             ONLY: section_release,&
                                              section_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE qs_active_space_types,           ONLY: eri_type_eri_element_func
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   TYPE, EXTENDS(eri_type_eri_element_func) :: eri2array
      INTEGER(C_INT), POINTER :: coords(:)
      REAL(C_DOUBLE), POINTER :: values(:)
      INTEGER                 :: idx = 1
   CONTAINS
      PROCEDURE :: func => eri2array_func
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param version_str ...
!> \param str_length ...
! **************************************************************************************************
   SUBROUTINE cp2k_get_version(version_str, str_length) BIND(C)
      CHARACTER(LEN=1, KIND=C_CHAR), INTENT(OUT)         :: version_str(*)
      INTEGER(C_INT), VALUE                              :: str_length

      INTEGER                                            :: i, n

      n = LEN_TRIM(cp2k_version)
      CPASSERT(str_length >= n + 1)
      MARK_USED(str_length)

      ! copy string
      DO i = 1, n
         version_str(i) = cp2k_version(i:i)
      ENDDO
      version_str(n + 1) = C_NULL_CHAR
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
      INTEGER(C_INT), INTENT(OUT)                        :: new_env_id
      CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN)          :: input_file_path(*), output_file_path(*)

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
      INTEGER(C_INT), INTENT(OUT)                        :: new_env_id
      CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN)          :: input_file_path(*), output_file_path(*)
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
      REAL(C_DOUBLE), DIMENSION(1:n_el), INTENT(IN)      :: new_pos

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
      REAL(C_DOUBLE), DIMENSION(1:n_el), INTENT(IN)      :: new_vel

      INTEGER                                            :: ierr

      CALL set_vel(env_id, new_vel, n_el, ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_set_velocities

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param new_cell ...
! **************************************************************************************************
   SUBROUTINE cp2k_set_cell(env_id, new_cell) BIND(C)
      INTEGER(C_INT), VALUE                              :: env_id
      REAL(C_DOUBLE), DIMENSION(3, 3), INTENT(IN)        :: new_cell

      INTEGER                                            :: ierr

      CALL set_cell(env_id, new_cell, ierr)
      CPASSERT(ierr == 0)
   END SUBROUTINE cp2k_set_cell

! **************************************************************************************************
!> \brief ...
!> \param env_id ...
!> \param description ...
!> \param RESULT ...
!> \param n_el ...
! **************************************************************************************************
   SUBROUTINE cp2k_get_result(env_id, description, RESULT, n_el) BIND(C)
      INTEGER(C_INT), VALUE                              :: env_id
      CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN)          :: description(*)
      INTEGER(C_INT), VALUE                              :: n_el
      REAL(C_DOUBLE), DIMENSION(1:n_el), INTENT(OUT)     :: RESULT

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
      INTEGER(C_INT), INTENT(OUT)                        :: natom

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
      INTEGER(C_INT), INTENT(OUT)                        :: nparticle

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
      REAL(C_DOUBLE), DIMENSION(1:n_el), INTENT(OUT)     :: pos

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
      REAL(C_DOUBLE), DIMENSION(1:n_el), INTENT(OUT)     :: force

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
      REAL(C_DOUBLE), INTENT(OUT)                        :: e_pot

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
      CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN)          :: input_file_path(*), output_file_path(*)

      CHARACTER(LEN=default_path_length)                 :: ifp, ofp
      TYPE(section_type), POINTER                        :: input_declaration

      ifp = " "; ofp = " "
      CALL strncpy_c2f(ifp, input_file_path)
      CALL strncpy_c2f(ofp, output_file_path)

      NULLIFY (input_declaration)
      CALL create_cp2k_root_section(input_declaration)
      CALL run_input(input_declaration, ifp, ofp, empty_initial_variables)
      CALL section_release(input_declaration)
   END SUBROUTINE cp2k_run_input

! **************************************************************************************************
!> \brief ...
!> \param input_file_path ...
!> \param output_file_path ...
!> \param mpi_comm ...
! **************************************************************************************************
   SUBROUTINE cp2k_run_input_comm(input_file_path, output_file_path, mpi_comm) BIND(C)
      CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN)          :: input_file_path(*), output_file_path(*)
      INTEGER(C_INT), VALUE                              :: mpi_comm

      CHARACTER(LEN=default_path_length)                 :: ifp, ofp
      TYPE(section_type), POINTER                        :: input_declaration

      ifp = " "; ofp = " "
      CALL strncpy_c2f(ifp, input_file_path)
      CALL strncpy_c2f(ofp, output_file_path)

      NULLIFY (input_declaration)
      CALL create_cp2k_root_section(input_declaration)
      CALL run_input(input_declaration, ifp, ofp, empty_initial_variables, mpi_comm)
      CALL section_release(input_declaration)
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
!> \brief Get the number of molecular orbitals
!> \param f_env_id  the force env id
!> \return The number of elements or -1 if unavailable
!> \author Tiziano Mueller
! **************************************************************************************************
   INTEGER(C_INT) FUNCTION cp2k_active_space_get_mo_count(f_env_id) RESULT(nmo) BIND(C)
      USE qs_active_space_types, ONLY: active_space_type
      USE qs_mo_types, ONLY: get_mo_set
      USE qs_environment_types, ONLY: get_qs_env
      INTEGER(C_INT), VALUE                              :: f_env_id

      INTEGER                                            :: ierr
      TYPE(active_space_type), POINTER                   :: active_space_env
      TYPE(f_env_type), POINTER                          :: f_env

      nmo = -1
      NULLIFY (f_env)

      CALL f_env_add_defaults(f_env_id, f_env)

      try: BLOCK
         CALL get_qs_env(f_env%force_env%qs_env, active_space=active_space_env)

         IF (.NOT. ASSOCIATED(active_space_env)) &
            EXIT try

         CALL get_mo_set(active_space_env%mos_active(1)%mo_set, nmo=nmo)
      END BLOCK try

      CALL f_env_rm_defaults(f_env, ierr)
      CPASSERT(ierr == 0)
   END FUNCTION cp2k_active_space_get_mo_count

! **************************************************************************************************
!> \brief Get the active space Fock sub-matrix (as a full matrix)
!> \param f_env_id the force env id
!> \param buf C array to write the data to
!> \param buf_len The length of the C array to write the data to (must be at least mo_count^2)
!> \return The number of elements written or -1 if unavailable or buffer too small
!> \author Tiziano Mueller
! **************************************************************************************************
   INTEGER(C_LONG) FUNCTION cp2k_active_space_get_fock_sub(f_env_id, buf, buf_len) RESULT(nelem) BIND(C)
      USE qs_active_space_types, ONLY: active_space_type
      USE qs_mo_types, ONLY: get_mo_set
      USE qs_environment_types, ONLY: get_qs_env
      INTEGER(C_INT), VALUE                              :: f_env_id
      INTEGER(C_LONG), VALUE                             :: buf_len
      REAL(C_DOUBLE), DIMENSION(0:buf_len-1), &
         INTENT(OUT)                                     :: buf

      INTEGER                                            :: i, ierr, j, norb
      REAL(C_DOUBLE)                                     :: mval
      TYPE(active_space_type), POINTER                   :: active_space_env
      TYPE(f_env_type), POINTER                          :: f_env

      nelem = -1
      NULLIFY (f_env)

      CALL f_env_add_defaults(f_env_id, f_env)

      try: BLOCK
         CALL get_qs_env(f_env%force_env%qs_env, active_space=active_space_env)

         IF (.NOT. ASSOCIATED(active_space_env)) &
            EXIT try

         CALL get_mo_set(active_space_env%mos_active(1)%mo_set, nmo=norb)

         IF (buf_len < norb*norb) &
            EXIT try

         DO i = 0, norb - 1
            DO j = 0, norb - 1
               CALL cp_fm_get_element(active_space_env%fock_sub(1)%matrix, i + 1, j + 1, mval)
               buf(norb*i + j) = mval
               buf(norb*j + i) = mval
            END DO
         END DO

         ! finished successfully, set number of written elements
         nelem = norb**norb
      END BLOCK try

      CALL f_env_rm_defaults(f_env, ierr)
      CPASSERT(ierr == 0)
   END FUNCTION cp2k_active_space_get_fock_sub

! **************************************************************************************************
!> \brief Get the number of non-zero elements of the ERI
!> \param f_env_id the force env id
!> \return The number of elements or -1 if unavailable
!> \author Tiziano Mueller
! **************************************************************************************************
   INTEGER(C_LONG) FUNCTION cp2k_active_space_get_eri_nze_count(f_env_id) RESULT(nze_count) BIND(C)
      USE qs_active_space_types, ONLY: active_space_type
      USE qs_environment_types, ONLY: get_qs_env
      INTEGER(C_INT), VALUE                              :: f_env_id

      INTEGER                                            :: ierr
      TYPE(active_space_type), POINTER                   :: active_space_env
      TYPE(f_env_type), POINTER                          :: f_env

      nze_count = -1
      NULLIFY (f_env)

      CALL f_env_add_defaults(f_env_id, f_env)

      try: BLOCK
         CALL get_qs_env(f_env%force_env%qs_env, active_space=active_space_env)

         IF (.NOT. ASSOCIATED(active_space_env)) &
            EXIT try

         nze_count = INT(active_space_env%eri%eri(1)%csr_mat%nze_total, KIND(nze_count))
      END BLOCK try

      CALL f_env_rm_defaults(f_env, ierr)
      CPASSERT(ierr == 0)
   END FUNCTION cp2k_active_space_get_eri_nze_count

! **************************************************************************************************
!> \brief Get the electron repulsion integrals (as a sparse tensor)
!> \param f_env_id the force env id
!> \param buf_coords C array to write the indizes (i,j,k,l) to
!> \param buf_coords_len size of the buffer, must be at least 4*nze_count
!> \param buf_values C array to write the values to
!> \param buf_values_len size of the buffer, must be at least nze_count
!> \return The number of elements written or -1 if unavailable or buffer too small
!> \author Tiziano Mueller
! **************************************************************************************************
   INTEGER(C_LONG) FUNCTION cp2k_active_space_get_eri(f_env_id, &
                                                      buf_coords, buf_coords_len, &
                                                      buf_values, buf_values_len) RESULT(nelem) BIND(C)
      USE qs_active_space_types, ONLY: active_space_type
      USE qs_mo_types, ONLY: get_mo_set
      USE qs_environment_types, ONLY: get_qs_env
      INTEGER(C_INT), INTENT(IN), VALUE                  :: f_env_id
      INTEGER(C_LONG), INTENT(IN), VALUE                 :: buf_coords_len
      INTEGER(C_INT), INTENT(OUT), TARGET                :: buf_coords(1:buf_coords_len)
      INTEGER(C_LONG), INTENT(IN), VALUE                 :: buf_values_len
      REAL(C_DOUBLE), INTENT(OUT), TARGET                :: buf_values(1:buf_values_len)

      INTEGER                                            :: ierr
      TYPE(active_space_type), POINTER                   :: active_space_env
      TYPE(f_env_type), POINTER                          :: f_env

      nelem = -1
      NULLIFY (f_env)

      CALL f_env_add_defaults(f_env_id, f_env)

      try: BLOCK
         CALL get_qs_env(f_env%force_env%qs_env, active_space=active_space_env)

         IF (.NOT. ASSOCIATED(active_space_env)) &
            EXIT try

         ASSOCIATE (nze=>active_space_env%eri%eri(1)%csr_mat%nze_total)
            IF (buf_coords_len < 4*nze .OR. buf_values_len < nze) &
               EXIT try

            CALL active_space_env%eri%eri_foreach(1, eri2array(buf_coords, buf_values))

            nelem = INT(nze, KIND(nelem))
         END ASSOCIATE
      END BLOCK try

      CALL f_env_rm_defaults(f_env, ierr)
      CPASSERT(ierr == 0)
   END FUNCTION cp2k_active_space_get_eri

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
!> \brief Copy the active space ERI to C buffers
!> \param this Class pointer
!> \param i The i index of the value `val`
!> \param j The j index of the value `val`
!> \param k The k index of the value `val`
!> \param l The l index of the value `val`
!> \param val The value at the given index
!> \return Always true to continue with the loop
!> \author Tiziano Mueller
! **************************************************************************************************
   LOGICAL FUNCTION eri2array_func(this, i, j, k, l, val) RESULT(cont)
      CLASS(eri2array), INTENT(inout) :: this
      INTEGER, INTENT(in)             :: i, j, k, l
      REAL(KIND=dp), INTENT(in)       :: val

      this%coords(4*(this%idx - 1) + 1) = i
      this%coords(4*(this%idx - 1) + 2) = j
      this%coords(4*(this%idx - 1) + 3) = k
      this%coords(4*(this%idx - 1) + 4) = l
      this%values(this%idx) = val

      this%idx = this%idx + 1

      cont = .TRUE.
   END FUNCTION eri2array_func

END MODULE libcp2k
