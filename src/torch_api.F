!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE torch_api
   USE ISO_C_BINDING, ONLY: C_ASSOCIATED, &
                            C_BOOL, &
                            C_CHAR, &
                            C_FLOAT, &
                            C_DOUBLE, &
                            C_F_POINTER, &
                            C_INT, &
                            C_NULL_CHAR, &
                            C_NULL_PTR, &
                            C_PTR, &
                            C_INT32_T, &
                            C_INT64_T

   USE kinds, ONLY: sp, int_4, int_8, dp, default_string_length

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   TYPE torch_tensor_type
      PRIVATE
      TYPE(C_PTR)                          :: c_ptr = C_NULL_PTR
   END TYPE torch_tensor_type

   TYPE torch_dict_type
      PRIVATE
      TYPE(C_PTR)                          :: c_ptr = C_NULL_PTR
   END TYPE torch_dict_type

   TYPE torch_model_type
      PRIVATE
      TYPE(C_PTR)                          :: c_ptr = C_NULL_PTR
   END TYPE torch_model_type

   #:set max_dim = 3
   INTERFACE torch_tensor_from_array
      #:for ndims  in range(1, max_dim+1)
         MODULE PROCEDURE torch_tensor_from_array_int32_${ndims}$d
         MODULE PROCEDURE torch_tensor_from_array_float_${ndims}$d
         MODULE PROCEDURE torch_tensor_from_array_int64_${ndims}$d
         MODULE PROCEDURE torch_tensor_from_array_double_${ndims}$d
      #:endfor
   END INTERFACE torch_tensor_from_array

   INTERFACE torch_tensor_upload_pinned
      #:for ndims  in range(1, max_dim+1)
         MODULE PROCEDURE torch_tensor_upload_pinned_double_${ndims}$d
      #:endfor
   END INTERFACE torch_tensor_upload_pinned

   INTERFACE staged_upload
      MODULE PROCEDURE staged_upload_1d, staged_upload_2d, staged_upload_3d
   END INTERFACE staged_upload

   INTERFACE torch_tensor_data_ptr
      #:for ndims  in range(1, max_dim+1)
         MODULE PROCEDURE torch_tensor_data_ptr_int32_${ndims}$d
         MODULE PROCEDURE torch_tensor_data_ptr_float_${ndims}$d
         MODULE PROCEDURE torch_tensor_data_ptr_int64_${ndims}$d
         MODULE PROCEDURE torch_tensor_data_ptr_double_${ndims}$d
      #:endfor
   END INTERFACE torch_tensor_data_ptr

   INTERFACE torch_model_get_attr
      MODULE PROCEDURE torch_model_get_attr_string
      MODULE PROCEDURE torch_model_get_attr_double
      MODULE PROCEDURE torch_model_get_attr_int64
      MODULE PROCEDURE torch_model_get_attr_int32
      MODULE PROCEDURE torch_model_get_attr_strlist
   END INTERFACE torch_model_get_attr

   PUBLIC :: torch_tensor_type, torch_tensor_expand_dim, torch_tensor_from_array, &
             torch_tensor_clone, torch_tensor_narrow, &
             torch_tensor_release
   PUBLIC :: staged_upload, torch_pinned_reserve, torch_stream_synchronize
   PUBLIC :: torch_file_stamp, torch_grad_pool_clear, torch_model_metadata_clear, &
             torch_pinned_finish, torch_pinned_release
   PUBLIC :: torch_tensor_data_ptr, torch_tensor_backward, torch_tensor_backward_scalar
   PUBLIC :: torch_tensor_grad, torch_tensor_grad_batch, torch_tensor_grad_batch3
   PUBLIC :: torch_tensor_hessian_vector
   PUBLIC :: torch_tensor_grad_batch_deferred, torch_tensor_grad_batch3_deferred
   PUBLIC :: torch_tensor_to_device_leaf, torch_tensor_set_requires_grad
   PUBLIC :: torch_tensor_accumulate, torch_tensor_item_double, torch_tensor_weighted_sum
   PUBLIC :: torch_tensor_weighted_sum_segments_total
   PUBLIC :: torch_dict_type, torch_dict_clone, torch_dict_create, torch_dict_insert
   PUBLIC :: torch_dict_get, torch_dict_release, torch_dict_set
   PUBLIC :: torch_get_num_threads, torch_set_default_threads
   PUBLIC :: torch_model_type, torch_model_load, torch_model_load_with_metadata, &
             torch_model_forward, torch_model_release
   PUBLIC :: torch_model_disable_parameter_gradients, torch_model_forward_mol_tensor, &
             torch_model_remap_device_constants
   PUBLIC :: torch_model_promote_float32
   PUBLIC :: torch_model_prepare_higher_derivatives
   PUBLIC :: torch_model_get_attr, torch_model_read_metadata
   PUBLIC :: torch_cuda_device_count, torch_cuda_is_available
   PUBLIC :: torch_allow_tf32, torch_model_freeze, torch_model_freeze_preserving_method, &
             torch_use_cuda

   LOGICAL, SAVE :: torch_default_threads_set = .FALSE.

CONTAINS

   #:set typenames = ['int32', 'float', 'int64', 'double']
   #:set types_f = ['INTEGER(kind=int_4)', 'REAL(sp)', 'INTEGER(kind=int_8)', 'REAL(dp)']
   #:set types_c = ['INTEGER(kind=C_INT32_T)', 'REAL(kind=C_FLOAT)', 'INTEGER(kind=C_INT64_T)', 'REAL(kind=C_DOUBLE)']

   #:for ndims in range(1, max_dim+1)
      #:for typename, type_f, type_c in zip(typenames, types_f, types_c)

! **************************************************************************************************
!> \brief Creates a Torch tensor from an array. The passed array has to outlive the tensor!
!>        The source must be an ALLOCATABLE to prevent passing a temporary array.
!> \author Ole Schuett
! **************************************************************************************************
         SUBROUTINE torch_tensor_from_array_${typename}$_${ndims}$d(tensor, source, requires_grad)
            TYPE(torch_tensor_type), INTENT(INOUT)             :: tensor
            #:set arraydims = ", ".join(":" for i in range(ndims))
            ${type_f}$, DIMENSION(${arraydims}$), ALLOCATABLE, INTENT(IN)  :: source
            LOGICAL, OPTIONAL, INTENT(IN)                      :: requires_grad

#if defined(__LIBTORCH)
            INTEGER(kind=int_8), DIMENSION(${ndims}$)          :: sizes_c
            LOGICAL                                            :: my_req_grad

            INTERFACE
               SUBROUTINE torch_c_tensor_from_array_${typename}$ (tensor, req_grad, ndims, sizes, source) &
                  BIND(C, name="torch_c_tensor_from_array_${typename}$")
                  IMPORT :: C_PTR, C_INT, C_INT32_T, C_INT64_T, C_FLOAT, C_DOUBLE, C_BOOL
                  TYPE(C_PTR)                                  :: tensor
                  LOGICAL(kind=C_BOOL), VALUE                  :: req_grad
                  INTEGER(kind=C_INT), VALUE                   :: ndims
                  INTEGER(kind=C_INT64_T), DIMENSION(*)        :: sizes
                  ${type_c}$, DIMENSION(*)                     :: source
               END SUBROUTINE torch_c_tensor_from_array_${typename}$
            END INTERFACE

            my_req_grad = .FALSE.
            IF (PRESENT(requires_grad)) my_req_grad = requires_grad

            #:for axis in range(ndims)
               sizes_c(${axis + 1}$) = SIZE(source, ${ndims - axis}$) ! C arrays are stored row-major.
            #:endfor

            CPASSERT(.NOT. C_ASSOCIATED(tensor%c_ptr))
            CALL torch_c_tensor_from_array_${typename}$ (tensor=tensor%c_ptr, &
                                                         req_grad=LOGICAL(my_req_grad, C_BOOL), &
                                                         ndims=${ndims}$, &
                                                         sizes=sizes_c, &
                                                         source=source)
            CPASSERT(C_ASSOCIATED(tensor%c_ptr))
#else
            CPABORT("CP2K compiled without the Torch library.")
            MARK_USED(tensor)
            MARK_USED(source)
            MARK_USED(requires_grad)
#endif
         END SUBROUTINE torch_tensor_from_array_${typename}$_${ndims}$d

! **************************************************************************************************
!> \brief Copies data from a Torch tensor to an array.
!>        The returned pointer is only valide during the tensor's lifetime!
!> \author Ole Schuett
! **************************************************************************************************
         SUBROUTINE torch_tensor_data_ptr_${typename}$_${ndims}$d(tensor, data_ptr)
            TYPE(torch_tensor_type), INTENT(IN)                :: tensor
            #:set arraydims = ", ".join(":" for i in range(ndims))
            ${type_f}$, DIMENSION(${arraydims}$), POINTER      :: data_ptr

#if defined(__LIBTORCH)
            INTEGER(kind=int_8), DIMENSION(${ndims}$)          :: sizes_f, sizes_c
            TYPE(C_PTR)                                        :: data_ptr_c

            INTERFACE
               SUBROUTINE torch_c_tensor_data_ptr_${typename}$ (tensor, ndims, sizes, data_ptr) &
                  BIND(C, name="torch_c_tensor_data_ptr_${typename}$")
                  IMPORT :: C_CHAR, C_PTR, C_INT, C_INT32_T, C_INT64_T
                  TYPE(C_PTR), VALUE                           :: tensor
                  INTEGER(kind=C_INT), VALUE                   :: ndims
                  INTEGER(kind=C_INT64_T), DIMENSION(*)        :: sizes
                  TYPE(C_PTR)                                  :: data_ptr
               END SUBROUTINE torch_c_tensor_data_ptr_${typename}$
            END INTERFACE

            sizes_c(:) = -1
            data_ptr_c = C_NULL_PTR
            CPASSERT(C_ASSOCIATED(tensor%c_ptr))
            CPASSERT(.NOT. ASSOCIATED(data_ptr))
            CALL torch_c_tensor_data_ptr_${typename}$ (tensor=tensor%c_ptr, &
                                                       ndims=${ndims}$, &
                                                       sizes=sizes_c, &
                                                       data_ptr=data_ptr_c)

            #:for axis in range(ndims)
               sizes_f(${axis + 1}$) = sizes_c(${ndims - axis}$) ! C arrays are stored row-major.
            #:endfor

            IF (ALL(sizes_f /= 0)) THEN  ! Torch returns null pointer for zero-sized tensors.
               CPASSERT(C_ASSOCIATED(data_ptr_c))
               CALL C_F_POINTER(data_ptr_c, data_ptr, shape=sizes_f)
            END IF
#else
            CPABORT("CP2K compiled without the Torch library.")
            MARK_USED(tensor)
            MARK_USED(data_ptr)
#endif
         END SUBROUTINE torch_tensor_data_ptr_${typename}$_${ndims}$d

      #:endfor
   #:endfor

   #:for ndims in range(1, max_dim+1)

! **************************************************************************************************
!> \brief Reuses or creates a device leaf tensor and uploads data through the shared
!>        pinned staging buffer. The device copy is asynchronous; call
!>        torch_stream_synchronize() before consuming the data.
! **************************************************************************************************
      SUBROUTINE torch_tensor_upload_pinned_double_${ndims}$d(tensor, source, offset, &
                                                              requires_grad, reallocated)
         TYPE(torch_tensor_type), INTENT(INOUT)             :: tensor
         #:set arraydims = ", ".join(":" for i in range(ndims))
         REAL(dp), DIMENSION(${arraydims}$), ALLOCATABLE, INTENT(IN)  :: source
         INTEGER(kind=int_8), INTENT(IN)                    :: offset
         LOGICAL, OPTIONAL, INTENT(IN)                      :: requires_grad
         LOGICAL, OPTIONAL, INTENT(OUT)                     :: reallocated

#if defined(__LIBTORCH)
         INTEGER(kind=int_8), DIMENSION(${ndims}$)          :: sizes_c
         LOGICAL                                            :: my_req_grad
         LOGICAL(kind=C_BOOL)                               :: reallocated_c

         INTERFACE
            FUNCTION torch_c_tensor_upload_pinned(tensor, offset, req_grad, ndims, sizes, source) &
               RESULT(reallocated_flag) BIND(C, name="torch_c_tensor_upload_pinned")
               IMPORT :: C_PTR, C_INT, C_INT64_T, C_DOUBLE, C_BOOL
               TYPE(C_PTR)                                  :: tensor
               INTEGER(kind=C_INT64_T), VALUE               :: offset
               LOGICAL(kind=C_BOOL), VALUE                  :: req_grad
               INTEGER(kind=C_INT), VALUE                   :: ndims
               INTEGER(kind=C_INT64_T), DIMENSION(*)        :: sizes
               REAL(kind=C_DOUBLE), DIMENSION(*)            :: source
               LOGICAL(kind=C_BOOL)                         :: reallocated_flag
            END FUNCTION torch_c_tensor_upload_pinned
         END INTERFACE

         my_req_grad = .FALSE.
         IF (PRESENT(requires_grad)) my_req_grad = requires_grad

         #:for axis in range(ndims)
            sizes_c(${axis + 1}$) = SIZE(source, ${ndims - axis}$) ! C arrays are stored row-major.
         #:endfor

         reallocated_c = torch_c_tensor_upload_pinned( &
            tensor=tensor%c_ptr, &
            offset=offset, &
            req_grad=LOGICAL(my_req_grad, C_BOOL), &
            ndims=${ndims}$, &
            sizes=sizes_c, &
            source=source)
         IF (PRESENT(reallocated)) reallocated = reallocated_c
         CPASSERT(C_ASSOCIATED(tensor%c_ptr))
#else
         CPABORT("CP2K compiled without the Torch library.")
         MARK_USED(tensor)
         MARK_USED(source)
         MARK_USED(offset)
         MARK_USED(requires_grad)
         MARK_USED(reallocated)
#endif
      END SUBROUTINE torch_tensor_upload_pinned_double_${ndims}$d

   #:endfor

! **************************************************************************************************
!> \brief Upload a 1D array through the staged path and return the next offset.
! **************************************************************************************************
   FUNCTION staged_upload_1d(tensor, array, offset, requires_grad, reallocated) RESULT(next_offset)
      TYPE(torch_tensor_type), INTENT(INOUT)             :: tensor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: array
      INTEGER(KIND=int_8), INTENT(IN)                    :: offset
      LOGICAL, INTENT(IN)                                :: requires_grad
      LOGICAL, OPTIONAL, INTENT(OUT)                     :: reallocated
      INTEGER(KIND=int_8)                                :: next_offset

      CALL torch_tensor_upload_pinned(tensor, array, offset, requires_grad, reallocated)
      next_offset = offset + INT(SIZE(array), KIND=int_8)

   END FUNCTION staged_upload_1d

! **************************************************************************************************
!> \brief Upload a 2D array through the staged path and return the next offset.
! **************************************************************************************************
   FUNCTION staged_upload_2d(tensor, array, offset, requires_grad, reallocated) RESULT(next_offset)
      TYPE(torch_tensor_type), INTENT(INOUT)             :: tensor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(IN) :: array
      INTEGER(KIND=int_8), INTENT(IN)                    :: offset
      LOGICAL, INTENT(IN)                                :: requires_grad
      LOGICAL, OPTIONAL, INTENT(OUT)                     :: reallocated
      INTEGER(KIND=int_8)                                :: next_offset

      CALL torch_tensor_upload_pinned(tensor, array, offset, requires_grad, reallocated)
      next_offset = offset + INT(SIZE(array), KIND=int_8)

   END FUNCTION staged_upload_2d

! **************************************************************************************************
!> \brief Upload a 3D array through the staged path and return the next offset.
! **************************************************************************************************
   FUNCTION staged_upload_3d(tensor, array, offset, requires_grad, reallocated) RESULT(next_offset)
      TYPE(torch_tensor_type), INTENT(INOUT)             :: tensor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), INTENT(IN) :: array
      INTEGER(KIND=int_8), INTENT(IN)                    :: offset
      LOGICAL, INTENT(IN)                                :: requires_grad
      LOGICAL, OPTIONAL, INTENT(OUT)                     :: reallocated
      INTEGER(KIND=int_8)                                :: next_offset

      CALL torch_tensor_upload_pinned(tensor, array, offset, requires_grad, reallocated)
      next_offset = offset + INT(SIZE(array), KIND=int_8)

   END FUNCTION staged_upload_3d

! **************************************************************************************************
!> \brief Creates an expanded tensor view along one singleton dimension.
! **************************************************************************************************
   SUBROUTINE torch_tensor_expand_dim(tensor, dim, extent, result)
      TYPE(torch_tensor_type), INTENT(IN)                :: tensor
      INTEGER, INTENT(IN)                                :: dim, extent
      TYPE(torch_tensor_type), INTENT(INOUT)             :: result

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_tensor_expand_dim(tensor, dim, extent, result) &
            BIND(C, name="torch_c_tensor_expand_dim")
            IMPORT :: C_INT64_T, C_PTR
            TYPE(C_PTR), VALUE                           :: tensor
            INTEGER(kind=C_INT64_T), VALUE               :: dim, extent
            TYPE(C_PTR)                                  :: result
         END SUBROUTINE torch_c_tensor_expand_dim
      END INTERFACE

      CPASSERT(C_ASSOCIATED(tensor%c_ptr))
      CPASSERT(.NOT. C_ASSOCIATED(result%c_ptr))
      CPASSERT(dim >= 0)
      CPASSERT(extent >= 0)
      CALL torch_c_tensor_expand_dim(tensor=tensor%c_ptr, &
                                     dim=INT(dim, C_INT64_T), &
                                     extent=INT(extent, C_INT64_T), &
                                     result=result%c_ptr)
      CPASSERT(C_ASSOCIATED(result%c_ptr))
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(tensor)
      MARK_USED(dim)
      MARK_USED(extent)
      MARK_USED(result)
#endif
   END SUBROUTINE torch_tensor_expand_dim

! **************************************************************************************************
!> \brief Creates a view of a contiguous tensor slice.
! **************************************************************************************************
   SUBROUTINE torch_tensor_narrow(tensor, dim, start_index, length, result)
      TYPE(torch_tensor_type), INTENT(IN)                :: tensor
      INTEGER, INTENT(IN)                                :: dim, start_index, length
      TYPE(torch_tensor_type), INTENT(INOUT)             :: result

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_tensor_narrow(tensor, dim, start_index, length, result) &
            BIND(C, name="torch_c_tensor_narrow")
            IMPORT :: C_INT64_T, C_PTR
            TYPE(C_PTR), VALUE                           :: tensor
            INTEGER(kind=C_INT64_T), VALUE               :: dim, start_index, length
            TYPE(C_PTR)                                  :: result
         END SUBROUTINE torch_c_tensor_narrow
      END INTERFACE

      CPASSERT(C_ASSOCIATED(tensor%c_ptr))
      CPASSERT(.NOT. C_ASSOCIATED(result%c_ptr))
      CPASSERT(dim >= 0)
      CPASSERT(start_index >= 0)
      CPASSERT(length >= 0)
      CALL torch_c_tensor_narrow(tensor=tensor%c_ptr, &
                                 dim=INT(dim, C_INT64_T), &
                                 start_index=INT(start_index, C_INT64_T), &
                                 length=INT(length, C_INT64_T), &
                                 result=result%c_ptr)
      CPASSERT(C_ASSOCIATED(result%c_ptr))
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(tensor)
      MARK_USED(dim)
      MARK_USED(start_index)
      MARK_USED(length)
      MARK_USED(result)
#endif
   END SUBROUTINE torch_tensor_narrow

! **************************************************************************************************
!> \brief Deep-copies a tensor into an independent tensor on the current device.
! **************************************************************************************************
   SUBROUTINE torch_tensor_clone(tensor, result)
      TYPE(torch_tensor_type), INTENT(IN)                :: tensor
      TYPE(torch_tensor_type), INTENT(INOUT)             :: result

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_tensor_clone(tensor, result) BIND(C, name="torch_c_tensor_clone")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                           :: tensor
            TYPE(C_PTR)                                  :: result
         END SUBROUTINE torch_c_tensor_clone
      END INTERFACE

      CPASSERT(C_ASSOCIATED(tensor%c_ptr))
      CPASSERT(.NOT. C_ASSOCIATED(result%c_ptr))
      CALL torch_c_tensor_clone(tensor%c_ptr, result%c_ptr)
      CPASSERT(C_ASSOCIATED(result%c_ptr))
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(tensor)
      MARK_USED(result)
#endif
   END SUBROUTINE torch_tensor_clone

! **************************************************************************************************
!> \brief Runs autograd on a Torch tensor.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_tensor_backward(tensor, outer_grad)
      TYPE(torch_tensor_type), INTENT(IN)                :: tensor
      TYPE(torch_tensor_type), INTENT(IN)                :: outer_grad

#if defined(__LIBTORCH)
      CHARACTER(len=*), PARAMETER                        :: routineN = 'torch_tensor_backward'
      INTEGER                                            :: handle

      INTERFACE
         SUBROUTINE torch_c_tensor_backward(tensor, outer_grad) &
            BIND(C, name="torch_c_tensor_backward")
            IMPORT :: C_CHAR, C_PTR
            TYPE(C_PTR), VALUE                           :: tensor
            TYPE(C_PTR), VALUE                           :: outer_grad
         END SUBROUTINE torch_c_tensor_backward
      END INTERFACE

      CALL timeset(routineN, handle)
      CPASSERT(C_ASSOCIATED(tensor%c_ptr))
      CPASSERT(C_ASSOCIATED(outer_grad%c_ptr))
      CALL torch_c_tensor_backward(tensor=tensor%c_ptr, outer_grad=outer_grad%c_ptr)
      CALL timestop(handle)
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(tensor)
      MARK_USED(outer_grad)
#endif
   END SUBROUTINE torch_tensor_backward

! **************************************************************************************************
!> \brief Runs autograd on a scalar Torch tensor.
! **************************************************************************************************
   SUBROUTINE torch_tensor_backward_scalar(tensor)
      TYPE(torch_tensor_type), INTENT(IN)                :: tensor

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_tensor_backward_scalar(tensor) &
            BIND(C, name="torch_c_tensor_backward_scalar")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                           :: tensor
         END SUBROUTINE torch_c_tensor_backward_scalar
      END INTERFACE

      CPASSERT(C_ASSOCIATED(tensor%c_ptr))
      CALL torch_c_tensor_backward_scalar(tensor=tensor%c_ptr)
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(tensor)
#endif
   END SUBROUTINE torch_tensor_backward_scalar

! **************************************************************************************************
!> \brief Applies the Hessian of a scalar tensor to the given direction tensors.
! **************************************************************************************************
   SUBROUTINE torch_tensor_hessian_vector(tensor, inputs, directions, responses, status, to_host)
      TYPE(torch_tensor_type), INTENT(IN)                :: tensor, inputs(:), directions(:)
      TYPE(torch_tensor_type), INTENT(INOUT)             :: responses(:)
      INTEGER, OPTIONAL, INTENT(OUT)                     :: status
      LOGICAL, OPTIONAL, INTENT(IN)                      :: to_host
#if defined(__LIBTORCH)
      TYPE(C_PTR)                                        :: input_ptrs(SIZE(inputs)), &
                                                            direction_ptrs(SIZE(inputs)), response_ptrs(SIZE(inputs))
      INTEGER                                            :: i
      INTEGER(C_INT)                                     :: c_status
      LOGICAL                                            :: host
      INTERFACE
         SUBROUTINE torch_c_tensor_hessian_vector(tensor, count, inputs, directions, responses, to_host, status) &
            BIND(C, name="torch_c_tensor_hessian_vector")
            IMPORT :: C_INT, C_PTR, C_BOOL
            TYPE(C_PTR), VALUE                           :: tensor
            INTEGER(C_INT), VALUE                        :: count
            TYPE(C_PTR), INTENT(IN)                      :: inputs(*), directions(*)
            TYPE(C_PTR), INTENT(OUT)                     :: responses(*)
            LOGICAL(C_BOOL), VALUE                       :: to_host
            INTEGER(C_INT), INTENT(OUT)                  :: status
         END SUBROUTINE torch_c_tensor_hessian_vector
      END INTERFACE
      CPASSERT(C_ASSOCIATED(tensor%c_ptr))
      CPASSERT(SIZE(inputs) > 0)
      CPASSERT(SIZE(directions) == SIZE(inputs))
      CPASSERT(SIZE(responses) == SIZE(inputs))
      DO i = 1, SIZE(inputs)
         CPASSERT(C_ASSOCIATED(inputs(i)%c_ptr))
         CPASSERT(C_ASSOCIATED(directions(i)%c_ptr))
         CPASSERT(.NOT. C_ASSOCIATED(responses(i)%c_ptr))
         input_ptrs(i) = inputs(i)%c_ptr
         direction_ptrs(i) = directions(i)%c_ptr
      END DO
      response_ptrs = C_NULL_PTR
      host = .TRUE.
      IF (PRESENT(to_host)) host = to_host
      CALL torch_c_tensor_hessian_vector(tensor%c_ptr, INT(SIZE(inputs), C_INT), &
                                         input_ptrs, direction_ptrs, response_ptrs, LOGICAL(host, C_BOOL), c_status)
      IF (PRESENT(status)) THEN
         status = c_status
      ELSE
         CPASSERT(c_status == 0)
      END IF
      DO i = 1, SIZE(inputs)
         responses(i)%c_ptr = response_ptrs(i)
         IF (c_status == 0) THEN
            CPASSERT(C_ASSOCIATED(responses(i)%c_ptr))
         END IF
      END DO
#else
      MARK_USED(tensor)
      MARK_USED(inputs)
      MARK_USED(directions)
      MARK_USED(responses)
      MARK_USED(status)
      MARK_USED(to_host)
      CPABORT("CP2K was compiled without LibTorch support.")
#endif
   END SUBROUTINE torch_tensor_hessian_vector

! **************************************************************************************************
!> \brief Moves a tensor to the active Torch device and makes it an autograd leaf.
! **************************************************************************************************
   SUBROUTINE torch_tensor_to_device_leaf(tensor, requires_grad)
      TYPE(torch_tensor_type), INTENT(INOUT)             :: tensor
      LOGICAL, INTENT(IN)                                :: requires_grad

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_tensor_to_device_leaf(tensor, req_grad) &
            BIND(C, name="torch_c_tensor_to_device_leaf")
            IMPORT :: C_BOOL, C_PTR
            TYPE(C_PTR)                                  :: tensor
            LOGICAL(kind=C_BOOL), VALUE                  :: req_grad
         END SUBROUTINE torch_c_tensor_to_device_leaf
      END INTERFACE

      CPASSERT(C_ASSOCIATED(tensor%c_ptr))
      CALL torch_c_tensor_to_device_leaf(tensor=tensor%c_ptr, &
                                         req_grad=LOGICAL(requires_grad, C_BOOL))
      CPASSERT(C_ASSOCIATED(tensor%c_ptr))
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(tensor)
      MARK_USED(requires_grad)
#endif
   END SUBROUTINE torch_tensor_to_device_leaf

! **************************************************************************************************
!> \brief Changes a tensor's requires_grad flag and clears any stale gradient.
! **************************************************************************************************
   SUBROUTINE torch_tensor_set_requires_grad(tensor, requires_grad)
      TYPE(torch_tensor_type), INTENT(INOUT)             :: tensor
      LOGICAL, INTENT(IN)                                :: requires_grad

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_tensor_set_requires_grad(tensor, req_grad) &
            BIND(C, name="torch_c_tensor_set_requires_grad")
            IMPORT :: C_BOOL, C_PTR
            TYPE(C_PTR), VALUE                           :: tensor
            LOGICAL(kind=C_BOOL), VALUE                  :: req_grad
         END SUBROUTINE torch_c_tensor_set_requires_grad
      END INTERFACE

      CPASSERT(C_ASSOCIATED(tensor%c_ptr))
      CALL torch_c_tensor_set_requires_grad(tensor=tensor%c_ptr, &
                                            req_grad=LOGICAL(requires_grad, C_BOOL))
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(tensor)
      MARK_USED(requires_grad)
#endif
   END SUBROUTINE torch_tensor_set_requires_grad

! **************************************************************************************************
!> \brief Select whether Torch wrappers should use CUDA when available.
! **************************************************************************************************
   SUBROUTINE torch_use_cuda(use_cuda)
      LOGICAL, INTENT(IN)                                :: use_cuda

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_use_cuda(use_cuda) BIND(C, name="torch_c_use_cuda")
            IMPORT :: C_BOOL
            LOGICAL(kind=C_BOOL), VALUE                  :: use_cuda
         END SUBROUTINE torch_c_use_cuda
      END INTERFACE

      CALL torch_c_use_cuda(use_cuda=LOGICAL(use_cuda, C_BOOL))
#else
      MARK_USED(use_cuda)
#endif
   END SUBROUTINE torch_use_cuda

! **************************************************************************************************
!> \brief Set the Torch CPU intra-op thread count (honours CP2K_TORCH_NUM_THREADS).
!> \param num_threads ...
! **************************************************************************************************
   SUBROUTINE torch_set_num_threads(num_threads)
      INTEGER, INTENT(IN)                                :: num_threads

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_set_num_threads(num_threads) BIND(C, name="torch_c_set_num_threads")
            IMPORT :: C_INT
            INTEGER(kind=C_INT), VALUE                   :: num_threads
         END SUBROUTINE torch_c_set_num_threads
      END INTERFACE

      CALL torch_c_set_num_threads(num_threads=INT(num_threads, C_INT))
#else
      MARK_USED(num_threads)
#endif
   END SUBROUTINE torch_set_num_threads

! **************************************************************************************************
!> \brief Return the current Torch CPU intra-op thread count.
!> \return ...
! **************************************************************************************************
   FUNCTION torch_get_num_threads() RESULT(num_threads)
      INTEGER                                            :: num_threads

#if defined(__LIBTORCH)
      INTERFACE
         FUNCTION torch_c_get_num_threads() BIND(C, name="torch_c_get_num_threads")
            IMPORT :: C_INT
            INTEGER(kind=C_INT)                          :: torch_c_get_num_threads
         END FUNCTION torch_c_get_num_threads
      END INTERFACE

      num_threads = INT(torch_c_get_num_threads())
#else
      num_threads = 1
#endif
   END FUNCTION torch_get_num_threads

! **************************************************************************************************
!> \brief Set the Torch intra-op thread count once, using the supplied per-process budget.
!> \param num_threads desired intra-op thread count
! **************************************************************************************************
   SUBROUTINE torch_set_default_threads(num_threads)
      INTEGER, INTENT(IN)                                :: num_threads

      IF (torch_default_threads_set) RETURN
      torch_default_threads_set = .TRUE.
      CALL torch_set_num_threads(num_threads)
   END SUBROUTINE torch_set_default_threads

! **************************************************************************************************
!> \brief Reserve capacity in the shared pinned host staging buffer.
! **************************************************************************************************
   SUBROUTINE torch_pinned_reserve(n_doubles)
      INTEGER(kind=int_8), INTENT(IN)                    :: n_doubles

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_pinned_reserve(n_doubles) BIND(C, name="torch_c_pinned_reserve")
            IMPORT :: C_INT64_T
            INTEGER(kind=C_INT64_T), VALUE               :: n_doubles
         END SUBROUTINE torch_c_pinned_reserve
      END INTERFACE

      CALL torch_c_pinned_reserve(n_doubles=n_doubles)
#else
      MARK_USED(n_doubles)
#endif
   END SUBROUTINE torch_pinned_reserve

! **************************************************************************************************
!> \brief Return this thread's pinned staging buffer to the pool.
! **************************************************************************************************
   SUBROUTINE torch_pinned_finish()

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_pinned_finish() BIND(C, name="torch_c_pinned_finish")
         END SUBROUTINE torch_c_pinned_finish
      END INTERFACE

      CALL torch_c_pinned_finish()
#endif
   END SUBROUTINE torch_pinned_finish

! **************************************************************************************************
!> \brief Release all pinned staging buffers.
! **************************************************************************************************
   SUBROUTINE torch_pinned_release()

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_pinned_release() BIND(C, name="torch_c_pinned_release")
         END SUBROUTINE torch_c_pinned_release
      END INTERFACE

      CALL torch_c_pinned_release()
#endif
   END SUBROUTINE torch_pinned_release

! **************************************************************************************************
!> \brief Wait for all pending asynchronous device copies.
! **************************************************************************************************
   SUBROUTINE torch_stream_synchronize()

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_stream_synchronize() BIND(C, name="torch_c_stream_synchronize")
         END SUBROUTINE torch_c_stream_synchronize
      END INTERFACE

      CALL torch_c_stream_synchronize()
#endif
   END SUBROUTINE torch_stream_synchronize

! **************************************************************************************************
!> \brief Return a size+mtime fingerprint of a model file (empty if unavailable).
! **************************************************************************************************
   SUBROUTINE torch_file_stamp(filename, stamp)
      CHARACTER(len=*), INTENT(IN)                       :: filename
      CHARACTER(len=*), INTENT(OUT)                      :: stamp

#if defined(__LIBTORCH)
      INTEGER(kind=C_INT)                                :: length

      INTERFACE
         SUBROUTINE torch_c_file_stamp(filename, buf, length) BIND(C, name="torch_c_file_stamp")
            IMPORT :: C_CHAR, C_INT
            CHARACTER(kind=C_CHAR), DIMENSION(*)         :: filename
            CHARACTER(kind=C_CHAR), DIMENSION(*)         :: buf
            INTEGER(kind=C_INT)                          :: length
         END SUBROUTINE torch_c_file_stamp
      END INTERFACE

      stamp = ""
      length = INT(LEN(stamp), C_INT)
      CALL torch_c_file_stamp(filename=TRIM(filename)//C_NULL_CHAR, buf=stamp, length=length)
#else
      stamp = ""
      MARK_USED(filename)
#endif
   END SUBROUTINE torch_file_stamp

! **************************************************************************************************
!> \brief Release the cached pinned gradient buffers.
! **************************************************************************************************
   SUBROUTINE torch_grad_pool_clear()

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_grad_pool_clear() BIND(C, name="torch_c_grad_pool_clear")
         END SUBROUTINE torch_c_grad_pool_clear
      END INTERFACE

      CALL torch_c_grad_pool_clear()
#endif
   END SUBROUTINE torch_grad_pool_clear

! **************************************************************************************************
!> \brief Clear the cached TorchScript model metadata.
! **************************************************************************************************
   SUBROUTINE torch_model_metadata_clear()

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_model_metadata_clear() BIND(C, name="torch_c_model_metadata_clear")
         END SUBROUTINE torch_c_model_metadata_clear
      END INTERFACE

      CALL torch_c_model_metadata_clear()
#endif
   END SUBROUTINE torch_model_metadata_clear

! **************************************************************************************************
!> \brief Returns the gradient of a Torch tensor which was computed by autograd.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_tensor_grad(tensor, grad)
      TYPE(torch_tensor_type), INTENT(IN)                :: tensor
      TYPE(torch_tensor_type), INTENT(INOUT)             :: grad

      TYPE(torch_tensor_type), DIMENSION(1)              :: grads, tensors

      CPASSERT(.NOT. C_ASSOCIATED(grad%c_ptr))
      tensors(1) = tensor
      CALL torch_tensor_grad_batch(tensors, grads)
      grad = grads(1)

   END SUBROUTINE torch_tensor_grad

! **************************************************************************************************
!> \brief Copies three autograd gradients to CPU memory.
! **************************************************************************************************
   SUBROUTINE torch_tensor_grad_batch3(tensor1, tensor2, tensor3, grad1, grad2, grad3)
      TYPE(torch_tensor_type), INTENT(IN)                :: tensor1, tensor2, tensor3
      TYPE(torch_tensor_type), INTENT(INOUT)             :: grad1, grad2, grad3

      CALL torch_tensor_grad_batch3_impl(tensor1, tensor2, tensor3, grad1, grad2, grad3, &
                                         deferred=.FALSE.)

   END SUBROUTINE torch_tensor_grad_batch3

! **************************************************************************************************
!> \brief Copies three autograd gradients to CPU memory without synchronizing.
! **************************************************************************************************
   SUBROUTINE torch_tensor_grad_batch3_deferred(tensor1, tensor2, tensor3, grad1, grad2, grad3)
      TYPE(torch_tensor_type), INTENT(IN)                :: tensor1, tensor2, tensor3
      TYPE(torch_tensor_type), INTENT(INOUT)             :: grad1, grad2, grad3

      CALL torch_tensor_grad_batch3_impl(tensor1, tensor2, tensor3, grad1, grad2, grad3, &
                                         deferred=.TRUE.)

   END SUBROUTINE torch_tensor_grad_batch3_deferred

! **************************************************************************************************
!> \brief Shared body of the three-gradient wrappers.
! **************************************************************************************************
   SUBROUTINE torch_tensor_grad_batch3_impl(tensor1, tensor2, tensor3, grad1, grad2, grad3, deferred)
      TYPE(torch_tensor_type), INTENT(IN)                :: tensor1, tensor2, tensor3
      TYPE(torch_tensor_type), INTENT(INOUT)             :: grad1, grad2, grad3
      LOGICAL, INTENT(IN)                                :: deferred

      TYPE(torch_tensor_type), DIMENSION(3)              :: grads, tensors

      tensors(1) = tensor1
      tensors(2) = tensor2
      tensors(3) = tensor3
      CALL torch_tensor_grad_batch_impl(tensors, grads, deferred)
      grad1 = grads(1)
      grad2 = grads(2)
      grad3 = grads(3)

   END SUBROUTINE torch_tensor_grad_batch3_impl

! **************************************************************************************************
!> \brief Copies several autograd gradients to CPU memory with one synchronization.
! **************************************************************************************************
   SUBROUTINE torch_tensor_grad_batch(tensors, grads)
      TYPE(torch_tensor_type), DIMENSION(:), INTENT(IN)  :: tensors
      TYPE(torch_tensor_type), DIMENSION(:), INTENT(INOUT) :: grads

      CALL torch_tensor_grad_batch_impl(tensors, grads, deferred=.FALSE.)

   END SUBROUTINE torch_tensor_grad_batch

! **************************************************************************************************
!> \brief Copies several autograd gradients to the pinned host buffers without synchronizing.
!>        Call torch_stream_synchronize() before reading the gradients.
! **************************************************************************************************
   SUBROUTINE torch_tensor_grad_batch_deferred(tensors, grads)
      TYPE(torch_tensor_type), DIMENSION(:), INTENT(IN)  :: tensors
      TYPE(torch_tensor_type), DIMENSION(:), INTENT(INOUT) :: grads

      CALL torch_tensor_grad_batch_impl(tensors, grads, deferred=.TRUE.)

   END SUBROUTINE torch_tensor_grad_batch_deferred

! **************************************************************************************************
!> \brief Shared implementation of the synchronized and deferred gradient batches.
! **************************************************************************************************
   SUBROUTINE torch_tensor_grad_batch_impl(tensors, grads, deferred)
      TYPE(torch_tensor_type), DIMENSION(:), INTENT(IN)  :: tensors
      TYPE(torch_tensor_type), DIMENSION(:), INTENT(INOUT) :: grads
      LOGICAL, INTENT(IN)                                :: deferred

#if defined(__LIBTORCH)
      INTEGER                                            :: i, n
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:)             :: c_grads, c_tensors

      INTERFACE
         SUBROUTINE torch_c_tensor_grad_batch(tensors, count, grads) &
            BIND(C, name="torch_c_tensor_grad_batch")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), DIMENSION(*)                    :: tensors
            INTEGER(kind=C_INT), VALUE                   :: count
            TYPE(C_PTR), DIMENSION(*)                    :: grads
         END SUBROUTINE torch_c_tensor_grad_batch
         SUBROUTINE torch_c_tensor_grad_batch_deferred(tensors, count, grads) &
            BIND(C, name="torch_c_tensor_grad_batch_deferred")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), DIMENSION(*)                    :: tensors
            INTEGER(kind=C_INT), VALUE                   :: count
            TYPE(C_PTR), DIMENSION(*)                    :: grads
         END SUBROUTINE torch_c_tensor_grad_batch_deferred
      END INTERFACE

      n = SIZE(tensors)
      CPASSERT(SIZE(grads) == n)
      ALLOCATE (c_tensors(n), c_grads(n))
      DO i = 1, n
         CPASSERT(C_ASSOCIATED(tensors(i)%c_ptr))
         c_tensors(i) = tensors(i)%c_ptr
         CPASSERT(.NOT. C_ASSOCIATED(grads(i)%c_ptr))
         c_grads(i) = C_NULL_PTR
      END DO
      IF (deferred) THEN
         CALL torch_c_tensor_grad_batch_deferred(tensors=c_tensors, count=INT(n, C_INT), &
                                                 grads=c_grads)
      ELSE
         CALL torch_c_tensor_grad_batch(tensors=c_tensors, count=INT(n, C_INT), grads=c_grads)
      END IF
      DO i = 1, n
         CPASSERT(C_ASSOCIATED(c_grads(i)))
         grads(i)%c_ptr = c_grads(i)
      END DO
      DEALLOCATE (c_tensors, c_grads)
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(tensors)
      MARK_USED(grads)
      MARK_USED(deferred)
#endif
   END SUBROUTINE torch_tensor_grad_batch_impl

! **************************************************************************************************
!> \brief Returns the weighted sum of two Torch tensors.
! **************************************************************************************************
   SUBROUTINE torch_tensor_weighted_sum(values, weights, result)
      TYPE(torch_tensor_type), INTENT(IN)                :: values, weights
      TYPE(torch_tensor_type), INTENT(INOUT)             :: result

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_tensor_weighted_sum(values, weights, result) &
            BIND(C, name="torch_c_tensor_weighted_sum")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                           :: values
            TYPE(C_PTR), VALUE                           :: weights
            TYPE(C_PTR)                                  :: result
         END SUBROUTINE torch_c_tensor_weighted_sum
      END INTERFACE

      CPASSERT(C_ASSOCIATED(values%c_ptr))
      CPASSERT(C_ASSOCIATED(weights%c_ptr))
      CPASSERT(.NOT. C_ASSOCIATED(result%c_ptr))
      CALL torch_c_tensor_weighted_sum(values=values%c_ptr, weights=weights%c_ptr, result=result%c_ptr)
      CPASSERT(C_ASSOCIATED(result%c_ptr))
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(values)
      MARK_USED(weights)
      MARK_USED(result)
#endif
   END SUBROUTINE torch_tensor_weighted_sum

! **************************************************************************************************
!> \brief Computes all weighted segment sums, returns them to host in one
!>        synchronization, and returns their total as a differentiable scalar.
! **************************************************************************************************
   SUBROUTINE torch_tensor_weighted_sum_segments_total(values, weights, offsets, lengths, total, result)
      TYPE(torch_tensor_type), INTENT(IN)                :: values, weights
      INTEGER(KIND=int_8), DIMENSION(:), INTENT(IN)      :: offsets, lengths
      TYPE(torch_tensor_type), INTENT(INOUT)             :: total
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: result

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_tensor_weighted_sum_segments_total(values, weights, offsets, lengths, nseg, &
                                                               result, total) &
            BIND(C, name="torch_c_tensor_weighted_sum_segments_total")
            IMPORT :: C_INT, C_INT64_T, C_DOUBLE, C_PTR
            TYPE(C_PTR), VALUE                           :: values
            TYPE(C_PTR), VALUE                           :: weights
            INTEGER(KIND=C_INT64_T), DIMENSION(*)        :: offsets
            INTEGER(KIND=C_INT64_T), DIMENSION(*)        :: lengths
            INTEGER(KIND=C_INT), VALUE                   :: nseg
            REAL(KIND=C_DOUBLE), DIMENSION(*)            :: result
            TYPE(C_PTR)                                  :: total
         END SUBROUTINE torch_c_tensor_weighted_sum_segments_total
      END INTERFACE

      CPASSERT(C_ASSOCIATED(values%c_ptr))
      CPASSERT(C_ASSOCIATED(weights%c_ptr))
      CPASSERT(SIZE(offsets) == SIZE(result))
      CPASSERT(SIZE(lengths) == SIZE(result))
      CPASSERT(.NOT. C_ASSOCIATED(total%c_ptr))
      CALL torch_c_tensor_weighted_sum_segments_total(values=values%c_ptr, weights=weights%c_ptr, &
                                                      offsets=offsets, lengths=lengths, &
                                                      nseg=INT(SIZE(result), C_INT), result=result, &
                                                      total=total%c_ptr)
      CPASSERT(C_ASSOCIATED(total%c_ptr))
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(values)
      MARK_USED(weights)
      MARK_USED(offsets)
      MARK_USED(lengths)
      MARK_USED(total)
      MARK_USED(result)
#endif
   END SUBROUTINE torch_tensor_weighted_sum_segments_total

! **************************************************************************************************
!> \brief Adds a device scalar to an accumulator without a host synchronization.
! **************************************************************************************************
   SUBROUTINE torch_tensor_accumulate(accumulator, value)
      TYPE(torch_tensor_type), INTENT(INOUT)             :: accumulator
      TYPE(torch_tensor_type), INTENT(IN)                :: value

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_tensor_accumulate(accumulator, value) &
            BIND(C, name="torch_c_tensor_accumulate")
            IMPORT :: C_PTR
            TYPE(C_PTR)                                  :: accumulator
            TYPE(C_PTR), VALUE                           :: value
         END SUBROUTINE torch_c_tensor_accumulate
      END INTERFACE

      CPASSERT(C_ASSOCIATED(value%c_ptr))
      CALL torch_c_tensor_accumulate(accumulator=accumulator%c_ptr, value=value%c_ptr)
      CPASSERT(C_ASSOCIATED(accumulator%c_ptr))
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(accumulator)
      MARK_USED(value)
#endif
   END SUBROUTINE torch_tensor_accumulate

! **************************************************************************************************
!> \brief Returns a scalar double value from a Torch tensor.
! **************************************************************************************************
   FUNCTION torch_tensor_item_double(tensor) RESULT(value)
      TYPE(torch_tensor_type), INTENT(IN)                :: tensor
      REAL(KIND=dp)                                      :: value

#if defined(__LIBTORCH)
      INTERFACE
         FUNCTION torch_c_tensor_item_double(tensor) RESULT(value) &
            BIND(C, name="torch_c_tensor_item_double")
            IMPORT :: C_DOUBLE, C_PTR
            TYPE(C_PTR), VALUE                           :: tensor
            REAL(KIND=C_DOUBLE)                          :: value
         END FUNCTION torch_c_tensor_item_double
      END INTERFACE

      CPASSERT(C_ASSOCIATED(tensor%c_ptr))
      value = torch_c_tensor_item_double(tensor=tensor%c_ptr)
#else
      value = 0.0_dp
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(tensor)
#endif
   END FUNCTION torch_tensor_item_double

! **************************************************************************************************
!> \brief Releases a Torch tensor and all its ressources.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_tensor_release(tensor)
      TYPE(torch_tensor_type), INTENT(INOUT)               :: tensor

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_tensor_release(tensor) BIND(C, name="torch_c_tensor_release")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: tensor
         END SUBROUTINE torch_c_tensor_release
      END INTERFACE

      CPASSERT(C_ASSOCIATED(tensor%c_ptr))
      CALL torch_c_tensor_release(tensor=tensor%c_ptr)
      tensor%c_ptr = C_NULL_PTR
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(tensor)
#endif
   END SUBROUTINE torch_tensor_release

! **************************************************************************************************
!> \brief Creates an empty Torch dictionary.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_dict_create(dict)
      TYPE(torch_dict_type), INTENT(INOUT)               :: dict

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_dict_create(dict) BIND(C, name="torch_c_dict_create")
            IMPORT :: C_PTR
            TYPE(C_PTR)                               :: dict
         END SUBROUTINE torch_c_dict_create
      END INTERFACE

      CPASSERT(.NOT. C_ASSOCIATED(dict%c_ptr))
      CALL torch_c_dict_create(dict=dict%c_ptr)
      CPASSERT(C_ASSOCIATED(dict%c_ptr))
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(dict)
#endif
   END SUBROUTINE torch_dict_create

! **************************************************************************************************
!> \brief Clones a Torch dictionary.
! **************************************************************************************************
   SUBROUTINE torch_dict_clone(source, target)
      TYPE(torch_dict_type), INTENT(IN)                  :: source
      TYPE(torch_dict_type), INTENT(INOUT)               :: target

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_dict_clone(source, target) BIND(C, name="torch_c_dict_clone")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: source
            TYPE(C_PTR)                               :: target
         END SUBROUTINE torch_c_dict_clone
      END INTERFACE

      CPASSERT(C_ASSOCIATED(source%c_ptr))
      CPASSERT(.NOT. C_ASSOCIATED(target%c_ptr))
      CALL torch_c_dict_clone(source=source%c_ptr, target=target%c_ptr)
      CPASSERT(C_ASSOCIATED(target%c_ptr))
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(source)
      MARK_USED(target)
#endif
   END SUBROUTINE torch_dict_clone

! **************************************************************************************************
!> \brief Inserts a Torch tensor into a Torch dictionary.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_dict_insert(dict, key, tensor)
      TYPE(torch_dict_type), INTENT(INOUT)               :: dict
      CHARACTER(len=*), INTENT(IN)                       :: key
      TYPE(torch_tensor_type), INTENT(IN)                :: tensor

#if defined(__LIBTORCH)

      INTERFACE
         SUBROUTINE torch_c_dict_insert(dict, key, tensor) &
            BIND(C, name="torch_c_dict_insert")
            IMPORT :: C_CHAR, C_PTR
            TYPE(C_PTR), VALUE                           :: dict
            CHARACTER(kind=C_CHAR), DIMENSION(*)         :: key
            TYPE(C_PTR), VALUE                           :: tensor
         END SUBROUTINE torch_c_dict_insert
      END INTERFACE

      CPASSERT(C_ASSOCIATED(dict%c_ptr))
      CPASSERT(C_ASSOCIATED(tensor%c_ptr))
      CALL torch_c_dict_insert(dict=dict%c_ptr, key=TRIM(key)//C_NULL_CHAR, tensor=tensor%c_ptr)
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(dict)
      MARK_USED(key)
      MARK_USED(tensor)
#endif
   END SUBROUTINE torch_dict_insert

! **************************************************************************************************
!> \brief Inserts or replaces a Torch tensor in a Torch dictionary.
! **************************************************************************************************
   SUBROUTINE torch_dict_set(dict, key, tensor)
      TYPE(torch_dict_type), INTENT(INOUT)               :: dict
      CHARACTER(len=*), INTENT(IN)                       :: key
      TYPE(torch_tensor_type), INTENT(IN)                :: tensor

#if defined(__LIBTORCH)

      INTERFACE
         SUBROUTINE torch_c_dict_set(dict, key, tensor) &
            BIND(C, name="torch_c_dict_set")
            IMPORT :: C_CHAR, C_PTR
            TYPE(C_PTR), VALUE                           :: dict
            CHARACTER(kind=C_CHAR), DIMENSION(*)         :: key
            TYPE(C_PTR), VALUE                           :: tensor
         END SUBROUTINE torch_c_dict_set
      END INTERFACE

      CPASSERT(C_ASSOCIATED(dict%c_ptr))
      CPASSERT(C_ASSOCIATED(tensor%c_ptr))
      CALL torch_c_dict_set(dict=dict%c_ptr, key=TRIM(key)//C_NULL_CHAR, tensor=tensor%c_ptr)
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(dict)
      MARK_USED(key)
      MARK_USED(tensor)
#endif
   END SUBROUTINE torch_dict_set

! **************************************************************************************************
!> \brief Retrieves a Torch tensor from a Torch dictionary.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_dict_get(dict, key, tensor)
      TYPE(torch_dict_type), INTENT(IN)                  :: dict
      CHARACTER(len=*), INTENT(IN)                       :: key
      TYPE(torch_tensor_type), INTENT(INOUT)             :: tensor

#if defined(__LIBTORCH)

      INTERFACE
         SUBROUTINE torch_c_dict_get(dict, key, tensor) &
            BIND(C, name="torch_c_dict_get")
            IMPORT :: C_CHAR, C_PTR
            TYPE(C_PTR), VALUE                           :: dict
            CHARACTER(kind=C_CHAR), DIMENSION(*)         :: key
            TYPE(C_PTR)                                  :: tensor
         END SUBROUTINE torch_c_dict_get
      END INTERFACE

      CPASSERT(C_ASSOCIATED(dict%c_ptr))
      CPASSERT(.NOT. C_ASSOCIATED(tensor%c_ptr))
      CALL torch_c_dict_get(dict=dict%c_ptr, key=TRIM(key)//C_NULL_CHAR, tensor=tensor%c_ptr)
      CPASSERT(C_ASSOCIATED(tensor%c_ptr))

#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(dict)
      MARK_USED(key)
      MARK_USED(tensor)
#endif
   END SUBROUTINE torch_dict_get

! **************************************************************************************************
!> \brief Releases a Torch dictionary and all its ressources.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_dict_release(dict)
      TYPE(torch_dict_type), INTENT(INOUT)               :: dict

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_dict_release(dict) BIND(C, name="torch_c_dict_release")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: dict
         END SUBROUTINE torch_c_dict_release
      END INTERFACE

      CPASSERT(C_ASSOCIATED(dict%c_ptr))
      CALL torch_c_dict_release(dict=dict%c_ptr)
      dict%c_ptr = C_NULL_PTR
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(dict)
#endif
   END SUBROUTINE torch_dict_release

! **************************************************************************************************
!> \brief Loads a Torch model from given "*.pth" file. (In Torch lingo models are called modules)
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_model_load(model, filename)
      TYPE(torch_model_type), INTENT(INOUT)              :: model
      CHARACTER(len=*), INTENT(IN)                       :: filename

#if defined(__LIBTORCH)
      CHARACTER(len=*), PARAMETER                        :: routineN = 'torch_model_load'
      INTEGER                                            :: handle

      INTERFACE
         SUBROUTINE torch_c_model_load(model, filename) BIND(C, name="torch_c_model_load")
            IMPORT :: C_PTR, C_CHAR
            TYPE(C_PTR)                               :: model
            CHARACTER(kind=C_CHAR), DIMENSION(*)      :: filename
         END SUBROUTINE torch_c_model_load
      END INTERFACE

      CALL timeset(routineN, handle)
      CPASSERT(.NOT. C_ASSOCIATED(model%c_ptr))
      CALL torch_c_model_load(model=model%c_ptr, filename=TRIM(filename)//C_NULL_CHAR)
      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CALL timestop(handle)
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(model)
      MARK_USED(filename)
#endif
   END SUBROUTINE torch_model_load

! **************************************************************************************************
!> \brief Loads a Torch model and reads two metadata entries in the same archive pass.
! **************************************************************************************************
   SUBROUTINE torch_model_load_with_metadata(model, filename, key1, value1, key2, value2)
      TYPE(torch_model_type), INTENT(INOUT)              :: model
      CHARACTER(len=*), INTENT(IN)                       :: filename, key1, key2
      CHARACTER(:), ALLOCATABLE, INTENT(OUT)             :: value1, value2

#if defined(__LIBTORCH)
      CHARACTER(len=*), PARAMETER                        :: routineN = 'torch_model_load_with_metadata'
      INTEGER                                            :: handle, length1, length2
      TYPE(C_PTR)                                        :: content1_c, content2_c

      INTERFACE
         SUBROUTINE torch_c_model_load_with_metadata(model, filename, key1, key2, &
                                                     content1, length1, content2, length2) &
            BIND(C, name="torch_c_model_load_with_metadata")
            IMPORT :: C_CHAR, C_INT, C_PTR
            TYPE(C_PTR)                                  :: model
            CHARACTER(kind=C_CHAR), DIMENSION(*)         :: filename, key1, key2
            TYPE(C_PTR)                                  :: content1, content2
            INTEGER(kind=C_INT)                          :: length1, length2
         END SUBROUTINE torch_c_model_load_with_metadata
      END INTERFACE

      CALL timeset(routineN, handle)
      CPASSERT(.NOT. C_ASSOCIATED(model%c_ptr))
      content1_c = C_NULL_PTR
      content2_c = C_NULL_PTR
      length1 = -1
      length2 = -1
      CALL torch_c_model_load_with_metadata(model=model%c_ptr, &
                                            filename=TRIM(filename)//C_NULL_CHAR, &
                                            key1=TRIM(key1)//C_NULL_CHAR, &
                                            key2=TRIM(key2)//C_NULL_CHAR, &
                                            content1=content1_c, length1=length1, &
                                            content2=content2_c, length2=length2)
      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CALL c_string_to_allocatable(content1_c, length1, value1)
      CALL c_string_to_allocatable(content2_c, length2, value2)
      CALL timestop(handle)
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(model)
      MARK_USED(filename)
      MARK_USED(key1)
      MARK_USED(value1)
      MARK_USED(key2)
      MARK_USED(value2)
#endif
   END SUBROUTINE torch_model_load_with_metadata

! **************************************************************************************************
!> \brief Maps serialized TorchScript device constants to the active Torch device.
! **************************************************************************************************
   SUBROUTINE torch_model_remap_device_constants(model)
      TYPE(torch_model_type), INTENT(INOUT)              :: model

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_model_remap_device_constants(model) &
            BIND(C, name="torch_c_model_remap_device_constants")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                           :: model
         END SUBROUTINE torch_c_model_remap_device_constants
      END INTERFACE

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CALL torch_c_model_remap_device_constants(model=model%c_ptr)
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(model)
#endif
   END SUBROUTINE torch_model_remap_device_constants

! **************************************************************************************************
!> \brief Promote float32 model state and explicit floating-point casts to double precision.
! **************************************************************************************************
   SUBROUTINE torch_model_promote_float32(model)
      TYPE(torch_model_type), INTENT(INOUT)              :: model

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_model_promote_float32(model) BIND(C, name="torch_c_model_promote_float32")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                           :: model
         END SUBROUTINE torch_c_model_promote_float32
      END INTERFACE

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CALL torch_c_model_promote_float32(model%c_ptr)
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(model)
#endif
   END SUBROUTINE torch_model_promote_float32

! **************************************************************************************************
!> \brief Prepare a model for higher input derivatives.
! **************************************************************************************************
   SUBROUTINE torch_model_prepare_higher_derivatives(model)
      TYPE(torch_model_type), INTENT(INOUT)              :: model

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_model_prepare_higher_derivatives(model) &
            BIND(C, name="torch_c_model_prepare_higher_derivatives")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                           :: model
         END SUBROUTINE torch_c_model_prepare_higher_derivatives
      END INTERFACE

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CALL torch_c_model_prepare_higher_derivatives(model=model%c_ptr)
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(model)
#endif
   END SUBROUTINE torch_model_prepare_higher_derivatives

! **************************************************************************************************
!> \brief Disable gradients for inference-only model parameters.
! **************************************************************************************************
   SUBROUTINE torch_model_disable_parameter_gradients(model)
      TYPE(torch_model_type), INTENT(INOUT)              :: model

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_model_disable_parameter_gradients(model) &
            BIND(C, name="torch_c_model_disable_parameter_gradients")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                           :: model
         END SUBROUTINE torch_c_model_disable_parameter_gradients
      END INTERFACE

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CALL torch_c_model_disable_parameter_gradients(model=model%c_ptr)
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(model)
#endif
   END SUBROUTINE torch_model_disable_parameter_gradients

! **************************************************************************************************
!> \brief Evaluates the given Torch model.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_model_forward(model, inputs, outputs)
      TYPE(torch_model_type), INTENT(INOUT)              :: model
      TYPE(torch_dict_type), INTENT(IN)                  :: inputs
      TYPE(torch_dict_type), INTENT(INOUT)               :: outputs

#if defined(__LIBTORCH)
      CHARACTER(len=*), PARAMETER                        :: routineN = 'torch_model_forward'
      INTEGER                                            :: handle

      INTERFACE
         SUBROUTINE torch_c_model_forward(model, inputs, outputs) BIND(C, name="torch_c_model_forward")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: model
            TYPE(C_PTR), VALUE                        :: inputs
            TYPE(C_PTR), VALUE                        :: outputs
         END SUBROUTINE torch_c_model_forward
      END INTERFACE

      CALL timeset(routineN, handle)
      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CPASSERT(C_ASSOCIATED(inputs%c_ptr))
      CPASSERT(C_ASSOCIATED(outputs%c_ptr))
      CALL torch_c_model_forward(model=model%c_ptr, inputs=inputs%c_ptr, outputs=outputs%c_ptr)
      CALL timestop(handle)
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(model)
      MARK_USED(inputs)
      MARK_USED(outputs)
#endif
   END SUBROUTINE torch_model_forward

! **************************************************************************************************
!> \brief Evaluates a TorchScript model method expecting keyword argument "mol".
! **************************************************************************************************
   SUBROUTINE torch_model_forward_mol_tensor(model, method_name, inputs, output)
      TYPE(torch_model_type), INTENT(INOUT)              :: model
      CHARACTER(len=*), INTENT(IN)                       :: method_name
      TYPE(torch_dict_type), INTENT(IN)                  :: inputs
      TYPE(torch_tensor_type), INTENT(INOUT)             :: output

#if defined(__LIBTORCH)
      CHARACTER(len=*), PARAMETER                        :: routineN = 'torch_model_forward_mol_tensor'
      INTEGER                                            :: handle

      INTERFACE
         SUBROUTINE torch_c_model_forward_mol_tensor(model, method_name, inputs, output) &
            BIND(C, name="torch_c_model_forward_mol_tensor")
            IMPORT :: C_CHAR, C_PTR
            TYPE(C_PTR), VALUE                           :: model
            CHARACTER(kind=C_CHAR), DIMENSION(*)         :: method_name
            TYPE(C_PTR), VALUE                           :: inputs
            TYPE(C_PTR)                                  :: output
         END SUBROUTINE torch_c_model_forward_mol_tensor
      END INTERFACE

      CALL timeset(routineN, handle)
      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CPASSERT(C_ASSOCIATED(inputs%c_ptr))
      CPASSERT(.NOT. C_ASSOCIATED(output%c_ptr))
      CALL torch_c_model_forward_mol_tensor(model=model%c_ptr, &
                                            method_name=TRIM(method_name)//C_NULL_CHAR, &
                                            inputs=inputs%c_ptr, &
                                            output=output%c_ptr)
      CPASSERT(C_ASSOCIATED(output%c_ptr))
      CALL timestop(handle)
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(model)
      MARK_USED(method_name)
      MARK_USED(inputs)
      MARK_USED(output)
#endif
   END SUBROUTINE torch_model_forward_mol_tensor

! **************************************************************************************************
!> \brief Releases a Torch model and all its ressources.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_model_release(model)
      TYPE(torch_model_type), INTENT(INOUT)              :: model

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_model_release(model) BIND(C, name="torch_c_model_release")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: model
         END SUBROUTINE torch_c_model_release
      END INTERFACE

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CALL torch_c_model_release(model=model%c_ptr)
      model%c_ptr = C_NULL_PTR
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(model)
#endif
   END SUBROUTINE torch_model_release

! **************************************************************************************************
!> \brief Reads metadata entry from given "*.pth" file. (In Torch lingo they are called extra files)
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION torch_model_read_metadata(filename, key) RESULT(res)
      CHARACTER(len=*), INTENT(IN)                       :: filename, key
      CHARACTER(:), ALLOCATABLE                           :: res

#if defined(__LIBTORCH)
      CHARACTER(len=*), PARAMETER                        :: routineN = 'torch_model_read_metadata'
      INTEGER                                            :: handle

      INTEGER                                            :: length
      TYPE(C_PTR)                                        :: content_c

      INTERFACE
         SUBROUTINE torch_c_model_read_metadata(filename, key, content, length) &
            BIND(C, name="torch_c_model_read_metadata")
            IMPORT :: C_CHAR, C_PTR, C_INT
            CHARACTER(kind=C_CHAR), DIMENSION(*)      :: filename, key
            TYPE(C_PTR)                               :: content
            INTEGER(kind=C_INT)                       :: length
         END SUBROUTINE torch_c_model_read_metadata
      END INTERFACE

      CALL timeset(routineN, handle)
      content_c = C_NULL_PTR
      length = -1
      CALL torch_c_model_read_metadata(filename=TRIM(filename)//C_NULL_CHAR, &
                                       key=TRIM(key)//C_NULL_CHAR, &
                                       content=content_c, &
                                       length=length)
      CALL c_string_to_allocatable(content_c, length, res)
      CALL timestop(handle)
#else
      res = ""
      MARK_USED(filename)
      MARK_USED(key)
      CPABORT("CP2K was compiled without Torch library.")
#endif
   END FUNCTION torch_model_read_metadata

! **************************************************************************************************
!> \brief Move a C-allocated null-terminated string into an allocatable Fortran string.
! **************************************************************************************************
   SUBROUTINE c_string_to_allocatable(content_c, length, res)
      TYPE(C_PTR), INTENT(INOUT)                         :: content_c
      INTEGER, INTENT(IN)                                :: length
      CHARACTER(:), ALLOCATABLE, INTENT(OUT)             :: res

#if defined(__LIBTORCH)
      CHARACTER(LEN=1, KIND=C_CHAR), DIMENSION(:), &
         POINTER                                         :: content_f
      INTEGER                                            :: i

      INTERFACE
         SUBROUTINE torch_c_free_string(content) BIND(C, name="torch_c_free_string")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: content
         END SUBROUTINE torch_c_free_string
      END INTERFACE

      CPASSERT(C_ASSOCIATED(content_c))
      CPASSERT(length >= 0)

      CALL C_F_POINTER(content_c, content_f, shape=[length + 1])
      CPASSERT(content_f(length + 1) == C_NULL_CHAR)

      ALLOCATE (CHARACTER(LEN=length) :: res)
      DO i = 1, length
         CPASSERT(content_f(i) /= C_NULL_CHAR)
         res(i:i) = content_f(i)
      END DO

      NULLIFY (content_f)
      CALL torch_c_free_string(content_c)
      content_c = C_NULL_PTR

#else
      res = ""
      MARK_USED(content_c)
      MARK_USED(length)
      CPABORT("CP2K was compiled without Torch library.")
#endif
   END SUBROUTINE c_string_to_allocatable

! **************************************************************************************************
!> \brief Returns true iff the Torch CUDA backend is available.
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION torch_cuda_is_available() RESULT(res)
      LOGICAL                                            :: res

#if defined(__LIBTORCH)
      INTERFACE
         FUNCTION torch_c_cuda_is_available() BIND(C, name="torch_c_cuda_is_available")
            IMPORT :: C_BOOL
            LOGICAL(C_BOOL)                           :: torch_c_cuda_is_available
         END FUNCTION torch_c_cuda_is_available
      END INTERFACE

      res = torch_c_cuda_is_available()
#else
      CPABORT("CP2K was compiled without Torch library.")
      res = .FALSE.
#endif
   END FUNCTION torch_cuda_is_available

! **************************************************************************************************
!> \brief Return the number of CUDA devices visible to Torch.
! **************************************************************************************************
   FUNCTION torch_cuda_device_count() RESULT(count)
      INTEGER                                            :: count

#if defined(__LIBTORCH)
      INTERFACE
         FUNCTION torch_c_cuda_device_count() BIND(C, name="torch_c_cuda_device_count")
            IMPORT :: C_INT
            INTEGER(C_INT)                            :: torch_c_cuda_device_count
         END FUNCTION torch_c_cuda_device_count
      END INTERFACE

      count = torch_c_cuda_device_count()
#else
      CPABORT("CP2K was compiled without Torch library.")
      count = 0
#endif
   END FUNCTION torch_cuda_device_count

! **************************************************************************************************
!> \brief Set whether to allow the use of TF32.
!>        Needed due to changes in defaults from pytorch 1.7 to 1.11 to >=1.12
!>        See https://pytorch.org/docs/stable/notes/cuda.html
!> \author Gabriele Tocci
! **************************************************************************************************
   SUBROUTINE torch_allow_tf32(allow_tf32)
      LOGICAL, INTENT(IN)                                  :: allow_tf32

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_allow_tf32(allow_tf32) BIND(C, name="torch_c_allow_tf32")
            IMPORT :: C_BOOL
            LOGICAL(C_BOOL), VALUE                  :: allow_tf32
         END SUBROUTINE torch_c_allow_tf32
      END INTERFACE

      CALL torch_c_allow_tf32(allow_tf32=LOGICAL(allow_tf32, C_BOOL))
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(allow_tf32)
#endif
   END SUBROUTINE torch_allow_tf32

! **************************************************************************************************
!> \brief Freeze the given Torch model: applies generic optimization that speed up model.
!>        See https://pytorch.org/docs/stable/generated/torch.jit.freeze.html
!> \author Gabriele Tocci
! **************************************************************************************************
   SUBROUTINE torch_model_freeze(model)
      TYPE(torch_model_type), INTENT(INOUT)              :: model

#if defined(__LIBTORCH)
      CHARACTER(len=*), PARAMETER                        :: routineN = 'torch_model_freeze'
      INTEGER                                            :: handle

      INTERFACE
         SUBROUTINE torch_c_model_freeze(model) BIND(C, name="torch_c_model_freeze")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: model
         END SUBROUTINE torch_c_model_freeze
      END INTERFACE

      CALL timeset(routineN, handle)
      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CALL torch_c_model_freeze(model=model%c_ptr)
      CALL timestop(handle)
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(model)
#endif
   END SUBROUTINE torch_model_freeze

! **************************************************************************************************
!> \brief Freeze a Torch model while preserving one exported method.
!> \param model ...
!> \param method_name ...
! **************************************************************************************************
   SUBROUTINE torch_model_freeze_preserving_method(model, method_name)
      TYPE(torch_model_type), INTENT(INOUT)              :: model
      CHARACTER(len=*), INTENT(IN)                       :: method_name

#if defined(__LIBTORCH)
      CHARACTER(len=*), PARAMETER                        :: routineN = &
                                                            'torch_model_freeze_preserving_method'
      INTEGER                                            :: handle

      INTERFACE
         SUBROUTINE torch_c_model_freeze_preserving_method(model, method_name) &
            BIND(C, name="torch_c_model_freeze_preserving_method")
            IMPORT :: C_CHAR, C_PTR
            TYPE(C_PTR), VALUE                        :: model
            CHARACTER(kind=C_CHAR), DIMENSION(*)      :: method_name
         END SUBROUTINE torch_c_model_freeze_preserving_method
      END INTERFACE

      CALL timeset(routineN, handle)
      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CALL torch_c_model_freeze_preserving_method( &
         model=model%c_ptr, method_name=TRIM(method_name)//C_NULL_CHAR)
      CALL timestop(handle)
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(method_name)
      MARK_USED(model)
#endif
   END SUBROUTINE torch_model_freeze_preserving_method

   #:set typenames = ['int64', 'double', 'string']
   #:set types_f = ['INTEGER(kind=int_8)', 'REAL(dp)', 'CHARACTER(LEN=default_string_length)']
   #:set types_c = ['INTEGER(kind=C_INT64_T)', 'REAL(kind=C_DOUBLE)', 'CHARACTER(kind=C_CHAR), DIMENSION(*)']
   #:set zeros_f = ['0', '0.0_dp', '""']

   #:for typename, type_f, type_c, zero_f in zip(typenames, types_f, types_c, zeros_f)
! **************************************************************************************************
!> \brief Retrieves an attribute from a Torch model. Must be called before torch_model_freeze.
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE torch_model_get_attr_${typename}$ (model, key, dest)
         TYPE(torch_model_type), INTENT(IN)                 :: model
         CHARACTER(len=*), INTENT(IN)                       :: key
         ${type_f}$, INTENT(OUT)                            :: dest

#if defined(__LIBTORCH)

         INTERFACE
            SUBROUTINE torch_c_model_get_attr_${typename}$ (model, key, dest) &
               BIND(C, name="torch_c_model_get_attr_${typename}$")
               IMPORT :: C_PTR, C_CHAR, C_INT64_T, C_DOUBLE
               TYPE(C_PTR), VALUE                           :: model
               CHARACTER(kind=C_CHAR), DIMENSION(*)         :: key
               ${type_c}$                                   :: dest
            END SUBROUTINE torch_c_model_get_attr_${typename}$
         END INTERFACE

         CALL torch_c_model_get_attr_${typename}$ (model=model%c_ptr, &
                                                   key=TRIM(key)//C_NULL_CHAR, &
                                                   dest=dest)
#else
         dest = ${zero_f}$
         MARK_USED(model)
         MARK_USED(key)
         CPABORT("CP2K compiled without the Torch library.")
#endif
      END SUBROUTINE torch_model_get_attr_${typename}$
   #:endfor

! **************************************************************************************************
!> \brief Retrieves an attribute from a Torch model. Must be called before torch_model_freeze.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_model_get_attr_int32(model, key, dest)
      TYPE(torch_model_type), INTENT(IN)                 :: model
      CHARACTER(len=*), INTENT(IN)                       :: key
      INTEGER, INTENT(OUT)                               :: dest

      INTEGER(kind=int_8)                                :: temp
      LOGICAL                                            :: in_range
      CALL torch_model_get_attr_int64(model, key, temp)
      in_range = temp >= -INT(HUGE(dest), KIND=int_8) - 1_int_8
      in_range = in_range .AND. temp <= INT(HUGE(dest), KIND=int_8)
      CPASSERT(in_range)
      dest = INT(temp)
   END SUBROUTINE torch_model_get_attr_int32

! **************************************************************************************************
!> \brief Retrieves a list attribute from a Torch model. Must be called before torch_model_freeze.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_model_get_attr_strlist(model, key, dest)
      TYPE(torch_model_type), INTENT(IN)                 :: model
      CHARACTER(len=*), INTENT(IN)                       :: key
      CHARACTER(LEN=default_string_length), &
         ALLOCATABLE, DIMENSION(:), INTENT(OUT)          :: dest

#if defined(__LIBTORCH)

      INTEGER :: num_items, i

      INTERFACE
         SUBROUTINE torch_c_model_get_attr_list_size(model, key, size) &
            BIND(C, name="torch_c_model_get_attr_list_size")
            IMPORT :: C_PTR, C_CHAR, C_INT
            TYPE(C_PTR), VALUE                           :: model
            CHARACTER(kind=C_CHAR), DIMENSION(*)         :: key
            INTEGER(kind=C_INT)                          :: size
         END SUBROUTINE torch_c_model_get_attr_list_size
      END INTERFACE

      INTERFACE
         SUBROUTINE torch_c_model_get_attr_strlist(model, key, index, dest) &
            BIND(C, name="torch_c_model_get_attr_strlist")
            IMPORT :: C_PTR, C_CHAR, C_INT
            TYPE(C_PTR), VALUE                           :: model
            CHARACTER(kind=C_CHAR), DIMENSION(*)         :: key
            INTEGER(kind=C_INT), VALUE                   :: index
            CHARACTER(kind=C_CHAR), DIMENSION(*)         :: dest
         END SUBROUTINE torch_c_model_get_attr_strlist
      END INTERFACE

      CALL torch_c_model_get_attr_list_size(model=model%c_ptr, &
                                            key=TRIM(key)//C_NULL_CHAR, &
                                            size=num_items)
      IF (num_items < 0) THEN
         CPABORT("Torch model metadata list returned a negative length.")
      END IF
      ALLOCATE (dest(num_items))
      dest(:) = ""

      DO i = 1, num_items
         CALL torch_c_model_get_attr_strlist(model=model%c_ptr, &
                                             key=TRIM(key)//C_NULL_CHAR, &
                                             index=i - 1, &
                                             dest=dest(i))

      END DO
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(model)
      MARK_USED(key)
      MARK_USED(dest)
#endif

   END SUBROUTINE torch_model_get_attr_strlist

END MODULE torch_api
