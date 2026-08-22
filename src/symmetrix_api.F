!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Fortran interface to libsymmetrix (bespoke C++ MACE kernels).
!>
!>        Thin ISO_C_BINDING layer over symmetrix_c_api.cpp. A MACE model is
!>        evaluated by handing symmetrix a CSR neighbour list (per-atom
!>        neighbour counts plus flattened neighbour indices/types, displacement
!>        vectors and distances); it returns the total energy, per-atom forces
!>        and the virial. symmetrix works in the model's native units (typically
!>        Angstrom / eV); the caller converts to/from atomic units.
!>
!>        Not thread-safe: symmetrix writes results into model-owned buffers, so
!>        a model handle must be driven from one thread at a time (the FIST
!>        manybody dispatch calls it serially per MPI rank).
!> \author Dhruv Sharma
! **************************************************************************************************
MODULE symmetrix_api
   USE ISO_C_BINDING,                   ONLY: C_ASSOCIATED,&
                                              C_CHAR,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_INT32_T,&
                                              C_NULL_CHAR,&
                                              C_NULL_PTR,&
                                              C_PTR
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: strlcpy_c2f
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: symmetrix_model_type
   PUBLIC :: symmetrix_model_load, symmetrix_model_release
   PUBLIC :: symmetrix_model_cutoff
   PUBLIC :: symmetrix_model_num_types, symmetrix_model_atomic_types
   PUBLIC :: symmetrix_eval
   PUBLIC :: symmetrix_is_kokkos

   INTEGER, PARAMETER :: symmetrix_err_len = 1024

   TYPE symmetrix_model_type
      PRIVATE
      TYPE(C_PTR) :: c_ptr = C_NULL_PTR
   END TYPE symmetrix_model_type

   INTERFACE
! **************************************************************************************************
!> \brief Loads a symmetrix MACE model (symmetrix_c_api.cpp).
!> \param path ...
!> \param model ...
!> \param status ...
!> \param err_msg ...
!> \param err_len ...
! **************************************************************************************************
      SUBROUTINE symmetrix_c_model_load(path, model, status, err_msg, err_len) &
         BIND(C, name="symmetrix_c_model_load")
         IMPORT :: C_PTR, C_CHAR, C_INT
         CHARACTER(kind=C_CHAR), DIMENSION(*)      :: path
         TYPE(C_PTR)                               :: model
         INTEGER(C_INT)                            :: status
         CHARACTER(kind=C_CHAR), DIMENSION(*)      :: err_msg
         INTEGER(C_INT), VALUE                     :: err_len
      END SUBROUTINE symmetrix_c_model_load

! **************************************************************************************************
!> \brief Neighbour cutoff of the model (native length unit).
!> \param model ...
!> \return ...
! **************************************************************************************************
      FUNCTION symmetrix_c_model_cutoff(model) &
         RESULT(r) BIND(C, name="symmetrix_c_model_cutoff")
         IMPORT :: C_PTR, C_DOUBLE
         TYPE(C_PTR), VALUE                        :: model
         REAL(C_DOUBLE)                            :: r
      END FUNCTION symmetrix_c_model_cutoff

! **************************************************************************************************
!> \brief Number of atomic types (elements) supported by the model.
!> \param model ...
!> \return ...
! **************************************************************************************************
      FUNCTION symmetrix_c_model_num_types(model) &
         RESULT(n) BIND(C, name="symmetrix_c_model_num_types")
         IMPORT :: C_PTR, C_INT
         TYPE(C_PTR), VALUE                        :: model
         INTEGER(C_INT)                            :: n
      END FUNCTION symmetrix_c_model_num_types

! **************************************************************************************************
!> \brief Copies the supported atomic numbers (Z) into out_types.
!> \param model ...
!> \param out_types ...
!> \param status ...
!> \param err_msg ...
!> \param err_len ...
! **************************************************************************************************
      SUBROUTINE symmetrix_c_model_types(model, out_types, status, err_msg, err_len) &
         BIND(C, name="symmetrix_c_model_types")
         IMPORT :: C_PTR, C_CHAR, C_INT, C_INT32_T
         TYPE(C_PTR), VALUE                        :: model
         INTEGER(C_INT32_T), DIMENSION(*)          :: out_types
         INTEGER(C_INT)                            :: status
         CHARACTER(kind=C_CHAR), DIMENSION(*)      :: err_msg
         INTEGER(C_INT), VALUE                     :: err_len
      END SUBROUTINE symmetrix_c_model_types

! **************************************************************************************************
!> \brief Evaluates over an owned-first node set: layer 0 over all n_total nodes,
!>        layer 1 and the readouts over the first n_local (owned) nodes only.
!> \param model ...
!> \param n_local ...
!> \param n_layer0 ...
!> \param n_total ...
!> \param node_types ...
!> \param n_edges ...
!> \param num_neigh ...
!> \param neigh_indices ...
!> \param neigh_types ...
!> \param edge_centre ...
!> \param xyz ...
!> \param r ...
!> \param do_virial ...
!> \param energy ...
!> \param forces ...
!> \param virial ...
!> \param status ...
!> \param err_msg ...
!> \param err_len ...
! **************************************************************************************************
      SUBROUTINE symmetrix_c_eval_mp(model, n_local, n_layer0, n_total, node_types, n_edges, num_neigh, &
                                     neigh_indices, neigh_types, edge_centre, xyz, r, &
                                     do_virial, energy, forces, virial, status, err_msg, err_len) &
         BIND(C, name="symmetrix_c_eval_mp")
         IMPORT :: C_PTR, C_CHAR, C_INT, C_INT32_T, C_DOUBLE
         TYPE(C_PTR), VALUE                        :: model
         INTEGER(C_INT), VALUE                     :: n_local
         INTEGER(C_INT), VALUE                     :: n_layer0
         INTEGER(C_INT), VALUE                     :: n_total
         INTEGER(C_INT32_T), DIMENSION(*)          :: node_types
         INTEGER(C_INT), VALUE                     :: n_edges
         INTEGER(C_INT32_T), DIMENSION(*)          :: num_neigh
         INTEGER(C_INT32_T), DIMENSION(*)          :: neigh_indices
         INTEGER(C_INT32_T), DIMENSION(*)          :: neigh_types
         INTEGER(C_INT32_T), DIMENSION(*)          :: edge_centre
         REAL(C_DOUBLE), DIMENSION(*)              :: xyz
         REAL(C_DOUBLE), DIMENSION(*)              :: r
         INTEGER(C_INT), VALUE                     :: do_virial
         REAL(C_DOUBLE)                            :: energy
         REAL(C_DOUBLE), DIMENSION(*)              :: forces
         REAL(C_DOUBLE), DIMENSION(*)              :: virial
         INTEGER(C_INT)                            :: status
         CHARACTER(kind=C_CHAR), DIMENSION(*)      :: err_msg
         INTEGER(C_INT), VALUE                     :: err_len
      END SUBROUTINE symmetrix_c_eval_mp

! **************************************************************************************************
!> \brief Frees a model handle.
!> \param model ...
! **************************************************************************************************
      SUBROUTINE symmetrix_c_model_release(model) &
         BIND(C, name="symmetrix_c_model_release")
         IMPORT :: C_PTR
         TYPE(C_PTR), VALUE                        :: model
      END SUBROUTINE symmetrix_c_model_release

! **************************************************************************************************
!> \brief Whether symmetrix was built with the Kokkos device backend.
!> \return ...
! **************************************************************************************************
      FUNCTION symmetrix_c_is_kokkos() &
         RESULT(k) BIND(C, name="symmetrix_c_is_kokkos")
         IMPORT :: C_INT
         INTEGER(C_INT)                            :: k
      END FUNCTION symmetrix_c_is_kokkos
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Convert a NUL-terminated C error buffer into a Fortran string.
!> \param buf the C character buffer
!> \return the Fortran string
! **************************************************************************************************
   FUNCTION symmetrix_err_string(buf) RESULT(msg)
      CHARACTER(kind=C_CHAR), DIMENSION(:), INTENT(IN)   :: buf
      CHARACTER(len=symmetrix_err_len)                   :: msg

      INTEGER                                            :: nchar

      nchar = strlcpy_c2f(msg, buf)
      IF (nchar < LEN(msg)) msg(nchar + 1:) = ""
   END FUNCTION symmetrix_err_string

! **************************************************************************************************
!> \brief Load a symmetrix MACE model from a .json file.
!> \param model handle, set on output
!> \param filename path to the symmetrix model file
! **************************************************************************************************
   SUBROUTINE symmetrix_model_load(model, filename)
      TYPE(symmetrix_model_type), INTENT(INOUT)          :: model
      CHARACTER(len=*), INTENT(IN)                       :: filename

#if defined(__SYMMETRIX)
      CHARACTER(kind=C_CHAR), DIMENSION(symmetrix_err_len) :: err_buf
      INTEGER(C_INT)                                     :: status

      CPASSERT(.NOT. C_ASSOCIATED(model%c_ptr))
      status = 0_C_INT
      CALL symmetrix_c_model_load(path=TRIM(filename)//C_NULL_CHAR, &
                                  model=model%c_ptr, status=status, &
                                  err_msg=err_buf, err_len=INT(symmetrix_err_len, C_INT))
      IF (status /= 0_C_INT) CPABORT(TRIM(symmetrix_err_string(err_buf)))
      CPASSERT(C_ASSOCIATED(model%c_ptr))
#else
      CPABORT("CP2K was compiled without symmetrix support.")
      MARK_USED(model)
      MARK_USED(filename)
#endif
   END SUBROUTINE symmetrix_model_load

! **************************************************************************************************
!> \brief Neighbour cutoff of the model (native length unit, usually Angstrom).
!> \param model loaded model handle
!> \return the cutoff
! **************************************************************************************************
   FUNCTION symmetrix_model_cutoff(model) RESULT(cutoff)
      TYPE(symmetrix_model_type), INTENT(IN)             :: model
      REAL(dp)                                           :: cutoff

#if defined(__SYMMETRIX)
      CPASSERT(C_ASSOCIATED(model%c_ptr))
      cutoff = symmetrix_c_model_cutoff(model%c_ptr)
#else
      cutoff = 0.0_dp
      CPABORT("CP2K was compiled without symmetrix support.")
      MARK_USED(model)
#endif
   END FUNCTION symmetrix_model_cutoff

! **************************************************************************************************
!> \brief Number of atomic types (elements) supported by the model.
!> \param model loaded model handle
!> \return the number of types
! **************************************************************************************************
   FUNCTION symmetrix_model_num_types(model) RESULT(n)
      TYPE(symmetrix_model_type), INTENT(IN)             :: model
      INTEGER                                            :: n

#if defined(__SYMMETRIX)
      CPASSERT(C_ASSOCIATED(model%c_ptr))
      n = INT(symmetrix_c_model_num_types(model%c_ptr))
#else
      n = 0
      CPABORT("CP2K was compiled without symmetrix support.")
      MARK_USED(model)
#endif
   END FUNCTION symmetrix_model_num_types

! **************************************************************************************************
!> \brief Atomic numbers supported by the model (size = symmetrix_model_num_types).
!> \param model loaded model handle
!> \param types the supported atomic numbers, filled on output
! **************************************************************************************************
   SUBROUTINE symmetrix_model_atomic_types(model, types)
      TYPE(symmetrix_model_type), INTENT(IN)             :: model
      INTEGER(C_INT32_T), DIMENSION(:), INTENT(OUT)      :: types

#if defined(__SYMMETRIX)
      CHARACTER(kind=C_CHAR), DIMENSION(symmetrix_err_len) :: err_buf
      INTEGER(C_INT)                                     :: status

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      status = 0_C_INT
      CALL symmetrix_c_model_types(model%c_ptr, types, status, err_buf, &
                                   INT(symmetrix_err_len, C_INT))
      IF (status /= 0_C_INT) CPABORT(TRIM(symmetrix_err_string(err_buf)))
#else
      types = 0_C_INT32_T
      CPABORT("CP2K was compiled without symmetrix support.")
      MARK_USED(model)
      MARK_USED(types)
#endif
   END SUBROUTINE symmetrix_model_atomic_types

! **************************************************************************************************
!> \brief Evaluate total energy, per-atom forces and the virial for a CSR neighbour list.
!> \param model loaded model handle
!> \param node_types (n_atoms) symmetrix 0-based type index per atom
!> \param num_neigh (n_atoms) number of edges whose central atom is this atom
!> \param neigh_indices (n_edges) neighbour (j) atom index per edge (0-based)
!> \param neigh_types (n_edges) symmetrix type index of the neighbour
!> \param edge_centre (n_edges) central (i) atom index per edge (0-based)
!> \param xyz (3, n_edges) displacement vectors r_j - r_i (model length unit)
!> \param r (n_edges) distances (model length unit)
!> \param do_virial compute the virial as well
!> \param energy total energy (model energy unit)
!> \param forces (3, n_atoms) forces (model energy/length unit)
!> \param virial (3, 3) virial (model energy unit); only meaningful when do_virial is .TRUE.
!> \param n_local ...
!> \param n_layer0 ...
! **************************************************************************************************
   SUBROUTINE symmetrix_eval(model, node_types, num_neigh, neigh_indices, neigh_types, &
                             edge_centre, xyz, r, do_virial, energy, forces, virial, &
                             n_local, n_layer0)
      TYPE(symmetrix_model_type), INTENT(IN)             :: model
      INTEGER(C_INT32_T), DIMENSION(:), INTENT(IN)       :: node_types, num_neigh
      INTEGER(C_INT32_T), DIMENSION(:), INTENT(IN)       :: neigh_indices, neigh_types, edge_centre
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: xyz
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: r
      LOGICAL, INTENT(IN)                                :: do_virial
      REAL(dp), INTENT(OUT)                              :: energy
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: forces
      REAL(dp), DIMENSION(3, 3), INTENT(OUT)             :: virial
      INTEGER, INTENT(IN), OPTIONAL                      :: n_local, n_layer0

#if defined(__SYMMETRIX)
      CHARACTER(kind=C_CHAR), DIMENSION(symmetrix_err_len) :: err_buf
      INTEGER                                            :: my_n_layer0, my_n_local
      INTEGER(C_INT)                                     :: status

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      ! the C interface carries 32-bit counts
      CPASSERT(SIZE(neigh_indices) <= HUGE(1_C_INT))
      my_n_local = SIZE(node_types)
      IF (PRESENT(n_local)) my_n_local = n_local
      my_n_layer0 = SIZE(node_types)
      IF (PRESENT(n_layer0)) my_n_layer0 = n_layer0
      CPASSERT(my_n_local >= 0 .AND. my_n_local <= my_n_layer0)
      CPASSERT(my_n_layer0 <= SIZE(node_types))
      energy = 0.0_dp
      forces = 0.0_dp
      virial = 0.0_dp
      status = 0_C_INT
      CALL symmetrix_c_eval_mp(model=model%c_ptr, &
                               n_local=INT(my_n_local, C_INT), &
                               n_layer0=INT(my_n_layer0, C_INT), &
                               n_total=INT(SIZE(node_types), C_INT), &
                               node_types=node_types, &
                               n_edges=INT(SIZE(neigh_indices), C_INT), &
                               num_neigh=num_neigh, neigh_indices=neigh_indices, &
                               neigh_types=neigh_types, edge_centre=edge_centre, &
                               xyz=xyz, r=r, &
                               do_virial=MERGE(1_C_INT, 0_C_INT, do_virial), &
                               energy=energy, forces=forces, virial=virial, &
                               status=status, err_msg=err_buf, &
                               err_len=INT(symmetrix_err_len, C_INT))
      IF (status /= 0_C_INT) CPABORT(TRIM(symmetrix_err_string(err_buf)))
#else
      energy = 0.0_dp
      forces = 0.0_dp
      virial = 0.0_dp
      CPABORT("CP2K was compiled without symmetrix support.")
      MARK_USED(model)
      MARK_USED(node_types)
      MARK_USED(num_neigh)
      MARK_USED(neigh_indices)
      MARK_USED(neigh_types)
      MARK_USED(edge_centre)
      MARK_USED(xyz)
      MARK_USED(r)
      MARK_USED(do_virial)
      MARK_USED(n_local)
      MARK_USED(n_layer0)
#endif
   END SUBROUTINE symmetrix_eval

! **************************************************************************************************
!> \brief Release a model handle.
!> \param model the model handle to release
! **************************************************************************************************
   SUBROUTINE symmetrix_model_release(model)
      TYPE(symmetrix_model_type), INTENT(INOUT)          :: model

#if defined(__SYMMETRIX)
      IF (C_ASSOCIATED(model%c_ptr)) THEN
         CALL symmetrix_c_model_release(model%c_ptr)
         model%c_ptr = C_NULL_PTR
      END IF
#else
      CPABORT("CP2K was compiled without symmetrix support.")
      MARK_USED(model)
#endif
   END SUBROUTINE symmetrix_model_release

! **************************************************************************************************
!> \brief Whether this binary evaluates symmetrix on a Kokkos device backend.
!>        The device is a build-time choice (CP2K_SYMMETRIX_KOKKOS), not a
!>        runtime one; input-level DEVICE requests are validated against it.
!> \return .TRUE. on a Kokkos (GPU) build, .FALSE. on a host build
! **************************************************************************************************
   FUNCTION symmetrix_is_kokkos() RESULT(kokkos)
      LOGICAL                                            :: kokkos

#if defined(__SYMMETRIX)
      kokkos = (symmetrix_c_is_kokkos() /= 0_C_INT)
#else
      kokkos = .FALSE.
      CPABORT("CP2K was compiled without symmetrix support.")
#endif
   END FUNCTION symmetrix_is_kokkos

END MODULE symmetrix_api
