!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
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
                            C_INT64_T

   USE kinds, ONLY: sp, int_8, dp, default_string_length

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
         MODULE PROCEDURE torch_tensor_from_array_float_${ndims}$d
         MODULE PROCEDURE torch_tensor_from_array_int64_${ndims}$d
         MODULE PROCEDURE torch_tensor_from_array_double_${ndims}$d
      #:endfor
   END INTERFACE torch_tensor_from_array

   INTERFACE torch_tensor_data_ptr
      #:for ndims  in range(1, max_dim+1)
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

   PUBLIC :: torch_tensor_type, torch_tensor_from_array, torch_tensor_release
   PUBLIC :: torch_tensor_data_ptr, torch_tensor_backward, torch_tensor_grad
   PUBLIC :: torch_dict_type, torch_dict_create, torch_dict_insert, torch_dict_get, torch_dict_release
   PUBLIC :: torch_model_type, torch_model_load, torch_model_forward, torch_model_release
   PUBLIC :: torch_model_get_attr, torch_model_read_metadata
   PUBLIC :: torch_cuda_is_available, torch_allow_tf32, torch_model_freeze

CONTAINS

   #:set typenames = ['float', 'int64', 'double']
   #:set types_f = ['REAL(sp)','INTEGER(kind=int_8)', 'REAL(dp)']
   #:set types_c = ['REAL(kind=C_FLOAT)','INTEGER(kind=C_INT64_T)', 'REAL(kind=C_DOUBLE)']

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
                  IMPORT :: C_PTR, C_INT, C_INT64_T, C_FLOAT, C_DOUBLE, C_BOOL
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
                  IMPORT :: C_CHAR, C_PTR, C_INT, C_INT64_T
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

            CPASSERT(ALL(sizes_c >= 0))
            CPASSERT(C_ASSOCIATED(data_ptr_c))

            #:for axis in range(ndims)
               sizes_f(${axis + 1}$) = sizes_c(${ndims - axis}$) ! C arrays are stored row-major.
            #:endfor
            CALL C_F_POINTER(data_ptr_c, data_ptr, shape=sizes_f)
#else
            CPABORT("CP2K compiled without the Torch library.")
            MARK_USED(tensor)
            MARK_USED(data_ptr)
#endif
         END SUBROUTINE torch_tensor_data_ptr_${typename}$_${ndims}$d

      #:endfor
   #:endfor

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
!> \brief Returns the gradient of a Torch tensor which was computed by autograd.
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_tensor_grad(tensor, grad)
      TYPE(torch_tensor_type), INTENT(IN)                :: tensor
      TYPE(torch_tensor_type), INTENT(INOUT)             :: grad

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_tensor_grad(tensor, grad) &
            BIND(C, name="torch_c_tensor_grad")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                           :: tensor
            TYPE(C_PTR)                                  :: grad
         END SUBROUTINE torch_c_tensor_grad
      END INTERFACE

      CPASSERT(C_ASSOCIATED(tensor%c_ptr))
      CPASSERT(.NOT. C_ASSOCIATED(grad%c_ptr))
      CALL torch_c_tensor_grad(tensor=tensor%c_ptr, grad=grad%c_ptr)
      CPASSERT(C_ASSOCIATED(grad%c_ptr))
#else
      CPABORT("CP2K compiled without the Torch library.")
      MARK_USED(tensor)
      MARK_USED(grad)
#endif
   END SUBROUTINE torch_tensor_grad

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

      CHARACTER(LEN=1, KIND=C_CHAR), DIMENSION(:), &
         POINTER                                         :: content_f
      INTEGER                                            :: i
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
      CPASSERT(C_ASSOCIATED(content_c))
      CPASSERT(length >= 0)

      CALL C_F_POINTER(content_c, content_f, shape=(/length + 1/))
      CPASSERT(content_f(length + 1) == C_NULL_CHAR)

      ALLOCATE (CHARACTER(LEN=length) :: res)
      DO i = 1, length
         CPASSERT(content_f(i) /= C_NULL_CHAR)
         res(i:i) = content_f(i)
      END DO

      DEALLOCATE (content_f) ! Was allocated on the C side.
      CALL timestop(handle)
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(filename)
      MARK_USED(key)
      MARK_USED(res)
#endif
   END FUNCTION torch_model_read_metadata

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

   #:set typenames = ['int64', 'double', 'string']
   #:set types_f = ['INTEGER(kind=int_8)', 'REAL(dp)', 'CHARACTER(LEN=default_string_length)']
   #:set types_c = ['INTEGER(kind=C_INT64_T)', 'REAL(kind=C_DOUBLE)', 'CHARACTER(kind=C_CHAR), DIMENSION(*)']

   #:for typename, type_f, type_c in zip(typenames, types_f, types_c)
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
         CPABORT("CP2K compiled without the Torch library.")
         MARK_USED(model)
         MARK_USED(key)
         MARK_USED(dest)
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
      CALL torch_model_get_attr_int64(model, key, temp)
      CPASSERT(ABS(temp) < HUGE(dest))
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
         ALLOCATABLE, DIMENSION(:)                       :: dest

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
