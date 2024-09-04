!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
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

   USE kinds, ONLY: sp, int_8, dp

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   TYPE torch_dict_type
      PRIVATE
      TYPE(C_PTR)                          :: c_ptr = C_NULL_PTR
   END TYPE torch_dict_type

   TYPE torch_model_type
      PRIVATE
      TYPE(C_PTR)                          :: c_ptr = C_NULL_PTR
   END TYPE torch_model_type

   #:set max_dim = 3
   INTERFACE torch_dict_insert
      #:for ndims  in range(1, max_dim+1)
         MODULE PROCEDURE torch_dict_insert_float_${ndims}$d
         MODULE PROCEDURE torch_dict_insert_int64_${ndims}$d
         MODULE PROCEDURE torch_dict_insert_double_${ndims}$d
      #:endfor
   END INTERFACE torch_dict_insert

   INTERFACE torch_dict_get
      #:for ndims  in range(1, max_dim+1)
         MODULE PROCEDURE torch_dict_get_float_${ndims}$d
         MODULE PROCEDURE torch_dict_get_int64_${ndims}$d
         MODULE PROCEDURE torch_dict_get_double_${ndims}$d
      #:endfor
   END INTERFACE torch_dict_get

   PUBLIC :: torch_dict_type, torch_dict_create, torch_dict_release
   PUBLIC :: torch_dict_insert, torch_dict_get
   PUBLIC :: torch_model_type, torch_model_load, torch_model_eval, torch_model_release
   PUBLIC :: torch_model_read_metadata
   PUBLIC :: torch_cuda_is_available, torch_allow_tf32, torch_model_freeze

CONTAINS

   #:set typenames = ['float', 'int64', 'double']
   #:set types_f = ['REAL(sp)','INTEGER(kind=int_8)', 'REAL(dp)']
   #:set types_c = ['REAL(kind=C_FLOAT)','INTEGER(kind=C_INT64_T)', 'REAL(kind=C_DOUBLE)']

   #:for ndims in range(1, max_dim+1)
      #:for typename, type_f, type_c in zip(typenames, types_f, types_c)

! **************************************************************************************************
!> \brief Inserts array into Torch dictionary. The passed array has to outlive the dictionary!
!> \author Ole Schuett
! **************************************************************************************************
         SUBROUTINE torch_dict_insert_${typename}$_${ndims}$d(dict, key, source)
            TYPE(torch_dict_type), INTENT(INOUT)               :: dict
            CHARACTER(len=*), INTENT(IN)                       :: key
            #:set arraydims = ", ".join(":" for i in range(ndims))
            ${type_f}$, CONTIGUOUS, DIMENSION(${arraydims}$), INTENT(IN)  :: source

#if defined(__LIBTORCH)
            INTEGER(kind=int_8), DIMENSION(${ndims}$)          :: sizes_c

            INTERFACE
               SUBROUTINE torch_c_dict_insert_${typename}$ (dict, key, ndims, sizes, source) &
                  BIND(C, name="torch_c_dict_insert_${typename}$")
                  IMPORT :: C_CHAR, C_PTR, C_INT, C_INT64_T, C_FLOAT, C_DOUBLE
                  TYPE(C_PTR), VALUE                           :: dict
                  CHARACTER(kind=C_CHAR), DIMENSION(*)         :: key
                  INTEGER(kind=C_INT), VALUE                   :: ndims
                  INTEGER(kind=C_INT64_T), DIMENSION(*)        :: sizes
                  ${type_c}$, DIMENSION(*)                     :: source
               END SUBROUTINE torch_c_dict_insert_${typename}$
            END INTERFACE

            #:for axis in range(ndims)
               sizes_c(${axis + 1}$) = SIZE(source, ${ndims - axis}$) ! C arrays are stored row-major.
            #:endfor

            CPASSERT(C_ASSOCIATED(dict%c_ptr))
            CALL torch_c_dict_insert_${typename}$ (dict=dict%c_ptr, &
                                                   key=TRIM(key)//C_NULL_CHAR, &
                                                   ndims=${ndims}$, &
                                                   sizes=sizes_c, &
                                                   source=source)
#else
            CPABORT("CP2K compiled without the Torch library.")
            MARK_USED(dict)
            MARK_USED(key)
            MARK_USED(source)
#endif
         END SUBROUTINE torch_dict_insert_${typename}$_${ndims}$d

! **************************************************************************************************
!> \brief Retrieves array from Torch dictionary. The returned array has to deallocated by caller!
!> \author Ole Schuett
! **************************************************************************************************
         SUBROUTINE torch_dict_get_${typename}$_${ndims}$d(dict, key, dest)
            TYPE(torch_dict_type), INTENT(IN)                  :: dict
            CHARACTER(len=*), INTENT(IN)                       :: key
            #:set arraydims = ", ".join(":" for i in range(ndims))
            ${type_f}$, DIMENSION(${arraydims}$), POINTER      :: dest

#if defined(__LIBTORCH)
            INTEGER(kind=int_8), DIMENSION(${ndims}$)          :: sizes_f, sizes_c
            TYPE(C_PTR)                                        :: dest_c

            INTERFACE
               SUBROUTINE torch_c_dict_get_${typename}$ (dict, key, ndims, sizes, dest) &
                  BIND(C, name="torch_c_dict_get_${typename}$")
                  IMPORT :: C_CHAR, C_PTR, C_INT, C_INT64_T
                  TYPE(C_PTR), VALUE                           :: dict
                  CHARACTER(kind=C_CHAR), DIMENSION(*)         :: key
                  INTEGER(kind=C_INT), VALUE                   :: ndims
                  INTEGER(kind=C_INT64_T), DIMENSION(*)        :: sizes
                  TYPE(C_PTR)                                  :: dest
               END SUBROUTINE torch_c_dict_get_${typename}$
            END INTERFACE

            sizes_c(:) = -1
            dest_c = C_NULL_PTR
            CPASSERT(C_ASSOCIATED(dict%c_ptr))
            CPASSERT(.NOT. ASSOCIATED(dest))
            CALL torch_c_dict_get_${typename}$ (dict=dict%c_ptr, &
                                                key=TRIM(key)//C_NULL_CHAR, &
                                                ndims=${ndims}$, &
                                                sizes=sizes_c, &
                                                dest=dest_c)

            CPASSERT(ALL(sizes_c >= 0))
            CPASSERT(C_ASSOCIATED(dest_c))

            #:for axis in range(ndims)
               sizes_f(${axis + 1}$) = sizes_c(${ndims - axis}$) ! C arrays are stored row-major.
            #:endfor
            CALL C_F_POINTER(dest_c, dest, shape=sizes_f)
#else
            CPABORT("CP2K compiled without the Torch library.")
            MARK_USED(dict)
            MARK_USED(key)
            MARK_USED(dest)
#endif
         END SUBROUTINE torch_dict_get_${typename}$_${ndims}$d

      #:endfor
   #:endfor

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
      INTERFACE
         SUBROUTINE torch_c_model_load(model, filename) BIND(C, name="torch_c_model_load")
            IMPORT :: C_PTR, C_CHAR
            TYPE(C_PTR)                               :: model
            CHARACTER(kind=C_CHAR), DIMENSION(*)      :: filename
         END SUBROUTINE torch_c_model_load
      END INTERFACE

      CPASSERT(.NOT. C_ASSOCIATED(model%c_ptr))
      CALL torch_c_model_load(model=model%c_ptr, filename=TRIM(filename)//C_NULL_CHAR)
      CPASSERT(C_ASSOCIATED(model%c_ptr))
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(model)
      MARK_USED(filename)
#endif
   END SUBROUTINE torch_model_load

! **************************************************************************************************
!> \brief Evaluates the given Torch model. (In Torch lingo this operation is called forward())
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE torch_model_eval(model, inputs, outputs)
      TYPE(torch_model_type), INTENT(INOUT)              :: model
      TYPE(torch_dict_type), INTENT(IN)                  :: inputs
      TYPE(torch_dict_type), INTENT(INOUT)               :: outputs

#if defined(__LIBTORCH)
      INTERFACE
         SUBROUTINE torch_c_model_eval(model, inputs, outputs) BIND(C, name="torch_c_model_eval")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: model
            TYPE(C_PTR), VALUE                        :: inputs
            TYPE(C_PTR), VALUE                        :: outputs
         END SUBROUTINE torch_c_model_eval
      END INTERFACE

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CPASSERT(C_ASSOCIATED(inputs%c_ptr))
      CPASSERT(C_ASSOCIATED(outputs%c_ptr))
      CALL torch_c_model_eval(model=model%c_ptr, &
                              inputs=inputs%c_ptr, &
                              outputs=outputs%c_ptr)
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(model)
      MARK_USED(inputs)
      MARK_USED(outputs)
#endif
   END SUBROUTINE torch_model_eval

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
      INTERFACE
         SUBROUTINE torch_c_model_freeze(model) BIND(C, name="torch_c_model_freeze")
            IMPORT :: C_PTR
            TYPE(C_PTR), VALUE                        :: model
         END SUBROUTINE torch_c_model_freeze
      END INTERFACE

      CPASSERT(C_ASSOCIATED(model%c_ptr))
      CALL torch_c_model_freeze(model=model%c_ptr)
#else
      CPABORT("CP2K was compiled without Torch library.")
      MARK_USED(model)
#endif
   END SUBROUTINE torch_model_freeze

END MODULE torch_api
