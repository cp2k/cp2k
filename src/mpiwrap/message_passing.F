!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the message passing library MPI
!> \par History
!>      JGH (02-Jan-2001): New error handling
!>                         Performance tools
!>      JGH (14-Jan-2001): New routines mp_comm_compare, mp_cart_coords,
!>                                      mp_rank_compare, mp_alltoall
!>      JGH (06-Feb-2001): New routines mp_comm_free
!>      JGH (22-Mar-2001): New routines mp_comm_dup
!>      fawzi (04-NOV-2004): storable performance info (for f77 interface)
!>      Wrapper routine for mpi_gatherv added (22.12.2005,MK)
!>      JGH (13-Feb-2006): Flexible precision
!>      JGH (15-Feb-2006): single precision mp_alltoall
!> \author JGH
! **************************************************************************************************
MODULE message_passing
   USE ISO_C_BINDING, ONLY: C_F_POINTER, &
                            C_PTR
   USE kinds, ONLY: &
      dp, int_4, int_4_size, int_8, int_8_size, real_4, real_4_size, real_8, &
      real_8_size, default_string_length
   USE machine, ONLY: m_abort

#include "../base/base_uses.f90"

! To simplify the transition between the old MPI module and the F08-style module, we introduce these constants to switch between the required handle types
! Unfortunately, Fortran does not offer something like typedef in C++
#if defined(__parallel) && defined(__MPI_F08)
#define MPI_DATA_TYPE TYPE(MPI_Datatype)
#define MPI_COMM_TYPE TYPE(MPI_Comm)
#define MPI_REQUEST_TYPE TYPE(MPI_Request)
#define MPI_WIN_TYPE TYPE(MPI_Win)
#define MPI_FILE_TYPE TYPE(MPI_File)
#define MPI_INFO_TYPE TYPE(MPI_Info)
#define MPI_STATUS_TYPE TYPE(MPI_Status)
#define MPI_GROUP_TYPE TYPE(MPI_Group)
#define MPI_STATUS_EXTRACT(X) %X
#else
#define MPI_DATA_TYPE INTEGER
#define MPI_COMM_TYPE INTEGER
#define MPI_REQUEST_TYPE INTEGER
#define MPI_WIN_TYPE INTEGER
#define MPI_FILE_TYPE INTEGER
#define MPI_INFO_TYPE INTEGER
#define MPI_STATUS_TYPE INTEGER, DIMENSION(MPI_STATUS_SIZE)
#define MPI_GROUP_TYPE INTEGER
#define MPI_STATUS_EXTRACT(X) (X)
#endif

#if defined(__parallel)
#if defined(__MPI_F08)
   USE mpi_f08
#else
   USE mpi
#endif
! subroutines: unfortunately, mpi implementations do not provide interfaces for all subroutines (problems with types and ranks explosion),
!              we do not quite know what is in the module, so we can not include any....
!              to nevertheless get checking for what is included, we use the mpi module without use clause, getting all there is
! USE mpi, ONLY: mpi_allgather, mpi_allgatherv, mpi_alloc_mem, mpi_allreduce, mpi_alltoall, mpi_alltoallv, mpi_bcast,&
!                mpi_cart_coords, mpi_cart_create, mpi_cart_get, mpi_cart_rank, mpi_cart_sub, mpi_dims_create, mpi_file_close,&
!                mpi_file_get_size, mpi_file_open, mpi_file_read_at_all, mpi_file_read_at, mpi_file_write_at_all,&
!                mpi_file_write_at, mpi_free_mem, mpi_gather, mpi_gatherv, mpi_get_address, mpi_group_translate_ranks, mpi_irecv,&
!                mpi_isend, mpi_recv, mpi_reduce, mpi_reduce_scatter, mpi_rget, mpi_scatter, mpi_send,&
!                mpi_sendrecv, mpi_sendrecv_replace, mpi_testany, mpi_waitall, mpi_waitany, mpi_win_create
! functions
! USE mpi, ONLY: mpi_wtime
! constants
! USE mpi, ONLY: MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX, MPI_REAL, MPI_COMPLEX, MPI_ANY_TAG,&
!                MPI_ANY_SOURCE, MPI_COMM_NULL, MPI_REQUEST_NULL, MPI_WIN_NULL, MPI_STATUS_SIZE, MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE, &
!                MPI_ADDRESS_KIND, MPI_OFFSET_KIND, MPI_MODE_CREATE, MPI_MODE_RDONLY, MPI_MODE_WRONLY,&
!                MPI_MODE_RDWR, MPI_MODE_EXCL, MPI_COMM_SELF, MPI_COMM_WORLD, MPI_THREAD_SERIALIZED,&
!                MPI_ERRORS_RETURN, MPI_SUCCESS, MPI_MAX_PROCESSOR_NAME, MPI_MAX_ERROR_STRING, MPI_IDENT,&
!                MPI_UNEQUAL, MPI_MAX, MPI_SUM, MPI_INFO_NULL, MPI_IN_PLACE, MPI_CONGRUENT, MPI_SIMILAR, MPI_MIN, MPI_SOURCE,&
!                MPI_TAG, MPI_INTEGER8, MPI_INTEGER, MPI_MAXLOC, MPI_2INTEGER, MPI_MINLOC, MPI_LOGICAL, MPI_2DOUBLE_PRECISION,&
!                MPI_LOR, MPI_CHARACTER, MPI_BOTTOM, MPI_MODE_NOCHECK, MPI_2REAL
#endif

   IMPLICIT NONE
   PRIVATE

   ! parameters that might be needed
#if defined(__parallel)
   LOGICAL, PARAMETER :: cp2k_is_parallel = .TRUE.
   INTEGER, PARAMETER, PUBLIC :: mp_any_tag = MPI_ANY_TAG
   INTEGER, PARAMETER, PUBLIC :: mp_any_source = MPI_ANY_SOURCE
   MPI_COMM_TYPE, PARAMETER :: mp_comm_null_handle = MPI_COMM_NULL
   MPI_COMM_TYPE, PARAMETER :: mp_comm_self_handle = MPI_COMM_SELF
   MPI_COMM_TYPE, PARAMETER :: mp_comm_world_handle = MPI_COMM_WORLD
   MPI_REQUEST_TYPE, PARAMETER :: mp_request_null_handle = MPI_REQUEST_NULL
   MPI_WIN_TYPE, PARAMETER :: mp_win_null_handle = MPI_WIN_NULL
   MPI_FILE_TYPE, PARAMETER :: mp_file_null_handle = MPI_FILE_NULL
   MPI_INFO_TYPE, PARAMETER :: mp_info_null_handle = MPI_INFO_NULL
   INTEGER, PARAMETER, PUBLIC :: mp_status_size = MPI_STATUS_SIZE
   INTEGER, PARAMETER, PUBLIC :: mp_proc_null = MPI_PROC_NULL
   ! Set max allocatable memory by MPI to 2 GiByte
   INTEGER(KIND=MPI_ADDRESS_KIND), PARAMETER, PRIVATE :: mp_max_memory_size = HUGE(INT(1, KIND=int_4))

   INTEGER, PARAMETER, PUBLIC :: mp_max_library_version_string = MPI_MAX_LIBRARY_VERSION_STRING

   INTEGER, PARAMETER, PUBLIC :: file_offset = MPI_OFFSET_KIND
   INTEGER, PARAMETER, PUBLIC :: address_kind = MPI_ADDRESS_KIND
   INTEGER, PARAMETER, PUBLIC :: file_amode_create = MPI_MODE_CREATE
   INTEGER, PARAMETER, PUBLIC :: file_amode_rdonly = MPI_MODE_RDONLY
   INTEGER, PARAMETER, PUBLIC :: file_amode_wronly = MPI_MODE_WRONLY
   INTEGER, PARAMETER, PUBLIC :: file_amode_rdwr = MPI_MODE_RDWR
   INTEGER, PARAMETER, PUBLIC :: file_amode_excl = MPI_MODE_EXCL
   INTEGER, PARAMETER, PUBLIC :: file_amode_append = MPI_MODE_APPEND
#else
   LOGICAL, PARAMETER :: cp2k_is_parallel = .FALSE.
   INTEGER, PARAMETER, PUBLIC :: mp_any_tag = -1
   INTEGER, PARAMETER, PUBLIC :: mp_any_source = -2
   MPI_COMM_TYPE, PARAMETER :: mp_comm_null_handle = -3
   MPI_COMM_TYPE, PARAMETER :: mp_comm_self_handle = -11
   MPI_COMM_TYPE, PARAMETER :: mp_comm_world_handle = -12
   MPI_REQUEST_TYPE, PARAMETER :: mp_request_null_handle = -4
   MPI_WIN_TYPE, PARAMETER :: mp_win_null_handle = -5
   MPI_FILE_TYPE, PARAMETER :: mp_file_null_handle = -6
   MPI_INFO_TYPE, PARAMETER :: mp_info_null_handle = -7
   INTEGER, PARAMETER, PUBLIC :: mp_status_size = -8
   INTEGER, PARAMETER, PUBLIC :: mp_proc_null = -9
   INTEGER, PARAMETER, PUBLIC :: mp_max_library_version_string = 1

   INTEGER, PARAMETER, PUBLIC :: file_offset = int_8
   INTEGER, PARAMETER, PUBLIC :: address_kind = int_8
   INTEGER, PARAMETER, PUBLIC :: file_amode_create = 1
   INTEGER, PARAMETER, PUBLIC :: file_amode_rdonly = 2
   INTEGER, PARAMETER, PUBLIC :: file_amode_wronly = 4
   INTEGER, PARAMETER, PUBLIC :: file_amode_rdwr = 8
   INTEGER, PARAMETER, PUBLIC :: file_amode_excl = 64
   INTEGER, PARAMETER, PUBLIC :: file_amode_append = 128
#endif

   ! we need to fix this to a given number (crossing fingers)
   ! so that the serial code using Fortran stream IO and the MPI have the same sizes.
   INTEGER, PARAMETER, PUBLIC :: mpi_character_size = 1
   INTEGER, PARAMETER, PUBLIC :: mpi_integer_size = 4

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = 'message_passing'

   ! internal reference counter used to debug communicator leaks
   INTEGER, PRIVATE, SAVE :: debug_comm_count

   PUBLIC :: mp_comm_type
   PUBLIC :: mp_request_type
   PUBLIC :: mp_win_type
   PUBLIC :: mp_file_type
   PUBLIC :: mp_info_type
   PUBLIC :: mp_cart_type

   TYPE mp_comm_type
      PRIVATE
      MPI_COMM_TYPE :: handle = mp_comm_null_handle
   CONTAINS
      ! Setters/Getters
      PROCEDURE, PASS, NON_OVERRIDABLE :: set_handle => mp_comm_type_set_handle
      PROCEDURE, PASS, NON_OVERRIDABLE :: get_handle => mp_comm_type_get_handle
      ! Comparisons
      PROCEDURE, PRIVATE, PASS, NON_OVERRIDABLE :: mp_comm_op_eq
      PROCEDURE, PRIVATE, PASS, NON_OVERRIDABLE :: mp_comm_op_neq
      GENERIC, PUBLIC :: operator(.EQ.) => mp_comm_op_eq
      GENERIC, PUBLIC :: operator(.NE.) => mp_comm_op_neq
      ! Communication routines

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: &
         mp_sendrecv_iv, mp_sendrecv_im2, mp_sendrecv_im3, mp_sendrecv_im4, &
         mp_sendrecv_lv, mp_sendrecv_lm2, mp_sendrecv_lm3, mp_sendrecv_lm4, &
         mp_sendrecv_rv, mp_sendrecv_rm2, mp_sendrecv_rm3, mp_sendrecv_rm4, &
         mp_sendrecv_dv, mp_sendrecv_dm2, mp_sendrecv_dm3, mp_sendrecv_dm4, &
         mp_sendrecv_cv, mp_sendrecv_cm2, mp_sendrecv_cm3, mp_sendrecv_cm4, &
         mp_sendrecv_zv, mp_sendrecv_zm2, mp_sendrecv_zm3, mp_sendrecv_zm4
      GENERIC, PUBLIC :: sendrecv => mp_sendrecv_iv, mp_sendrecv_im2, mp_sendrecv_im3, mp_sendrecv_im4, &
         mp_sendrecv_lv, mp_sendrecv_lm2, mp_sendrecv_lm3, mp_sendrecv_lm4, &
         mp_sendrecv_rv, mp_sendrecv_rm2, mp_sendrecv_rm3, mp_sendrecv_rm4, &
         mp_sendrecv_dv, mp_sendrecv_dm2, mp_sendrecv_dm3, mp_sendrecv_dm4, &
         mp_sendrecv_cv, mp_sendrecv_cm2, mp_sendrecv_cm3, mp_sendrecv_cm4, &
         mp_sendrecv_zv, mp_sendrecv_zm2, mp_sendrecv_zm3, mp_sendrecv_zm4

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_minloc_iv, &
         mp_minloc_lv, mp_minloc_rv, mp_minloc_dv
      GENERIC, PUBLIC :: minloc => mp_minloc_iv, &
         mp_minloc_lv, mp_minloc_rv, mp_minloc_dv

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_maxloc_iv, &
         mp_maxloc_lv, mp_maxloc_rv, mp_maxloc_dv
      GENERIC, PUBLIC :: maxloc => mp_maxloc_iv, &
         mp_maxloc_lv, mp_maxloc_rv, mp_maxloc_dv

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_shift_im, mp_shift_i, &
         mp_shift_lm, mp_shift_l, mp_shift_rm, mp_shift_r, &
         mp_shift_dm, mp_shift_d, mp_shift_cm, mp_shift_c, &
         mp_shift_zm, mp_shift_z
      GENERIC, PUBLIC :: shift => mp_shift_im, mp_shift_i, &
         mp_shift_lm, mp_shift_l, mp_shift_rm, mp_shift_r, &
         mp_shift_dm, mp_shift_d, mp_shift_cm, mp_shift_c, &
         mp_shift_zm, mp_shift_z

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_bcast_i, mp_bcast_iv, mp_bcast_im, mp_bcast_i3, &
         mp_bcast_l, mp_bcast_lv, mp_bcast_lm, mp_bcast_l3, &
         mp_bcast_r, mp_bcast_rv, mp_bcast_rm, mp_bcast_r3, &
         mp_bcast_d, mp_bcast_dv, mp_bcast_dm, mp_bcast_d3, &
         mp_bcast_c, mp_bcast_cv, mp_bcast_cm, mp_bcast_c3, &
         mp_bcast_z, mp_bcast_zv, mp_bcast_zm, mp_bcast_z3, &
         mp_bcast_b, mp_bcast_bv, mp_bcast_av, mp_bcast_am
      GENERIC, PUBLIC :: bcast => mp_bcast_i, mp_bcast_iv, mp_bcast_im, mp_bcast_i3, &
         mp_bcast_l, mp_bcast_lv, mp_bcast_lm, mp_bcast_l3, &
         mp_bcast_r, mp_bcast_rv, mp_bcast_rm, mp_bcast_r3, &
         mp_bcast_d, mp_bcast_dv, mp_bcast_dm, mp_bcast_d3, &
         mp_bcast_c, mp_bcast_cv, mp_bcast_cm, mp_bcast_c3, &
         mp_bcast_z, mp_bcast_zv, mp_bcast_zm, mp_bcast_z3, &
         mp_bcast_b, mp_bcast_bv, mp_bcast_av, mp_bcast_am

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_ibcast_i, mp_ibcast_iv, &
         mp_ibcast_l, mp_ibcast_lv, mp_ibcast_r, mp_ibcast_rv, &
         mp_ibcast_d, mp_ibcast_dv, mp_ibcast_c, mp_ibcast_cv, &
         mp_ibcast_z, mp_ibcast_zv
      GENERIC, PUBLIC :: ibcast => mp_ibcast_i, mp_ibcast_iv, &
         mp_ibcast_l, mp_ibcast_lv, mp_ibcast_r, mp_ibcast_rv, &
         mp_ibcast_d, mp_ibcast_dv, mp_ibcast_c, mp_ibcast_cv, &
         mp_ibcast_z, mp_ibcast_zv

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: &
         mp_sum_i, mp_sum_iv, mp_sum_im, mp_sum_im3, mp_sum_im4, &
         mp_sum_l, mp_sum_lv, mp_sum_lm, mp_sum_lm3, mp_sum_lm4, &
         mp_sum_r, mp_sum_rv, mp_sum_rm, mp_sum_rm3, mp_sum_rm4, &
         mp_sum_d, mp_sum_dv, mp_sum_dm, mp_sum_dm3, mp_sum_dm4, &
         mp_sum_c, mp_sum_cv, mp_sum_cm, mp_sum_cm3, mp_sum_cm4, &
         mp_sum_z, mp_sum_zv, mp_sum_zm, mp_sum_zm3, mp_sum_zm4, &
         mp_sum_root_iv, mp_sum_root_im, mp_sum_root_lv, mp_sum_root_lm, &
         mp_sum_root_rv, mp_sum_root_rm, mp_sum_root_dv, mp_sum_root_dm, &
         mp_sum_root_cv, mp_sum_root_cm, mp_sum_root_zv, mp_sum_root_zm, &
         mp_sum_b, mp_sum_bv
      GENERIC, PUBLIC :: sum => mp_sum_i, mp_sum_iv, mp_sum_im, mp_sum_im3, mp_sum_im4, &
         mp_sum_l, mp_sum_lv, mp_sum_lm, mp_sum_lm3, mp_sum_lm4, &
         mp_sum_r, mp_sum_rv, mp_sum_rm, mp_sum_rm3, mp_sum_rm4, &
         mp_sum_d, mp_sum_dv, mp_sum_dm, mp_sum_dm3, mp_sum_dm4, &
         mp_sum_c, mp_sum_cv, mp_sum_cm, mp_sum_cm3, mp_sum_cm4, &
         mp_sum_z, mp_sum_zv, mp_sum_zm, mp_sum_zm3, mp_sum_zm4, &
         mp_sum_root_iv, mp_sum_root_im, mp_sum_root_lv, mp_sum_root_lm, &
         mp_sum_root_rv, mp_sum_root_rm, mp_sum_root_dv, mp_sum_root_dm, &
         mp_sum_root_cv, mp_sum_root_cm, mp_sum_root_zv, mp_sum_root_zm, &
         mp_sum_b, mp_sum_bv

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_isum_iv, &
         mp_isum_lv, mp_isum_rv, mp_isum_dv, mp_isum_cv, &
         mp_isum_zv, mp_isum_bv
      GENERIC, PUBLIC :: isum => mp_isum_iv, &
         mp_isum_lv, mp_isum_rv, mp_isum_dv, mp_isum_cv, &
         mp_isum_zv, mp_isum_bv

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_sum_partial_im, &
         mp_sum_partial_lm, mp_sum_partial_rm, mp_sum_partial_dm, &
         mp_sum_partial_cm, mp_sum_partial_zm
      GENERIC, PUBLIC :: sum_partial => mp_sum_partial_im, &
         mp_sum_partial_lm, mp_sum_partial_rm, mp_sum_partial_dm, &
         mp_sum_partial_cm, mp_sum_partial_zm

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_max_i, mp_max_iv, &
         mp_max_l, mp_max_lv, mp_max_r, mp_max_rv, &
         mp_max_d, mp_max_dv, mp_max_c, mp_max_cv, &
         mp_max_z, mp_max_zv
      GENERIC, PUBLIC :: max => mp_max_i, mp_max_iv, &
         mp_max_l, mp_max_lv, mp_max_r, mp_max_rv, &
         mp_max_d, mp_max_dv, mp_max_c, mp_max_cv, &
         mp_max_z, mp_max_zv

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_min_i, mp_min_iv, &
         mp_min_l, mp_min_lv, mp_min_r, mp_min_rv, &
         mp_min_d, mp_min_dv, mp_min_c, mp_min_cv, &
         mp_min_z, mp_min_zv
      GENERIC, PUBLIC :: min => mp_min_i, mp_min_iv, &
         mp_min_l, mp_min_lv, mp_min_r, mp_min_rv, &
         mp_min_d, mp_min_dv, mp_min_c, mp_min_cv, &
         mp_min_z, mp_min_zv

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_prod_r, mp_prod_d, mp_prod_c, mp_prod_z
      GENERIC, PUBLIC :: prod => mp_prod_r, mp_prod_d, mp_prod_c, mp_prod_z

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_gather_i, mp_gather_iv, mp_gather_im, &
         mp_gather_l, mp_gather_lv, mp_gather_lm, &
         mp_gather_r, mp_gather_rv, mp_gather_rm, &
         mp_gather_d, mp_gather_dv, mp_gather_dm, &
         mp_gather_c, mp_gather_cv, mp_gather_cm, &
         mp_gather_z, mp_gather_zv, mp_gather_zm
      GENERIC, PUBLIC :: gather => mp_gather_i, mp_gather_iv, mp_gather_im, &
         mp_gather_l, mp_gather_lv, mp_gather_lm, &
         mp_gather_r, mp_gather_rv, mp_gather_rm, &
         mp_gather_d, mp_gather_dv, mp_gather_dm, &
         mp_gather_c, mp_gather_cv, mp_gather_cm, &
         mp_gather_z, mp_gather_zv, mp_gather_zm

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_gatherv_iv, &
         mp_gatherv_lv, mp_gatherv_rv, mp_gatherv_dv, &
         mp_gatherv_cv, mp_gatherv_zv, mp_gatherv_lm2, mp_gatherv_rm2, &
         mp_gatherv_dm2, mp_gatherv_cm2, mp_gatherv_zm2
      GENERIC, PUBLIC :: gatherv => mp_gatherv_iv, &
         mp_gatherv_lv, mp_gatherv_rv, mp_gatherv_dv, &
         mp_gatherv_cv, mp_gatherv_zv, mp_gatherv_lm2, mp_gatherv_rm2, &
         mp_gatherv_dm2, mp_gatherv_cm2, mp_gatherv_zm2

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_igatherv_iv, &
         mp_igatherv_lv, mp_igatherv_rv, mp_igatherv_dv, &
         mp_igatherv_cv, mp_igatherv_zv
      GENERIC, PUBLIC :: igatherv => mp_igatherv_iv, &
         mp_igatherv_lv, mp_igatherv_rv, mp_igatherv_dv, &
         mp_igatherv_cv, mp_igatherv_zv

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_allgather_i, mp_allgather_i2, &
         mp_allgather_i12, mp_allgather_i23, mp_allgather_i34, &
         mp_allgather_i22, mp_allgather_l, mp_allgather_l2, &
         mp_allgather_l12, mp_allgather_l23, mp_allgather_l34, &
         mp_allgather_l22, mp_allgather_r, mp_allgather_r2, &
         mp_allgather_r12, mp_allgather_r23, mp_allgather_r34, &
         mp_allgather_r22, mp_allgather_d, mp_allgather_d2, &
         mp_allgather_d12, mp_allgather_d23, mp_allgather_d34, &
         mp_allgather_d22, mp_allgather_c, mp_allgather_c2, &
         mp_allgather_c12, mp_allgather_c23, mp_allgather_c34, &
         mp_allgather_c22, mp_allgather_z, mp_allgather_z2, &
         mp_allgather_z12, mp_allgather_z23, mp_allgather_z34, &
         mp_allgather_z22
      GENERIC, PUBLIC :: allgather => mp_allgather_i, mp_allgather_i2, &
         mp_allgather_i12, mp_allgather_i23, mp_allgather_i34, &
         mp_allgather_i22, mp_allgather_l, mp_allgather_l2, &
         mp_allgather_l12, mp_allgather_l23, mp_allgather_l34, &
         mp_allgather_l22, mp_allgather_r, mp_allgather_r2, &
         mp_allgather_r12, mp_allgather_r23, mp_allgather_r34, &
         mp_allgather_r22, mp_allgather_d, mp_allgather_d2, &
         mp_allgather_d12, mp_allgather_d23, mp_allgather_d34, &
         mp_allgather_d22, mp_allgather_c, mp_allgather_c2, &
         mp_allgather_c12, mp_allgather_c23, mp_allgather_c34, &
         mp_allgather_c22, mp_allgather_z, mp_allgather_z2, &
         mp_allgather_z12, mp_allgather_z23, mp_allgather_z34, &
         mp_allgather_z22

      PROCEDURE, PRIVATE, PASS(comm) ::  mp_allgatherv_iv, mp_allgatherv_lv, &
         mp_allgatherv_rv, mp_allgatherv_dv, mp_allgatherv_cv, mp_allgatherv_zv, &
         mp_allgatherv_im2, mp_allgatherv_lm2, mp_allgatherv_rm2, &
         mp_allgatherv_dm2, mp_allgatherv_cm2, mp_allgatherv_zm2
      GENERIC, PUBLIC :: allgatherv => mp_allgatherv_iv, mp_allgatherv_lv, &
         mp_allgatherv_rv, mp_allgatherv_dv, mp_allgatherv_cv, mp_allgatherv_zv, &
         mp_allgatherv_im2, mp_allgatherv_lm2, mp_allgatherv_rm2, &
         mp_allgatherv_dm2, mp_allgatherv_cm2, mp_allgatherv_zm2

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_iallgather_i, mp_iallgather_l, &
         mp_iallgather_r, mp_iallgather_d, mp_iallgather_c, mp_iallgather_z, &
         mp_iallgather_i11, mp_iallgather_l11, mp_iallgather_r11, mp_iallgather_d11, &
         mp_iallgather_c11, mp_iallgather_z11, mp_iallgather_i13, mp_iallgather_l13, &
         mp_iallgather_r13, mp_iallgather_d13, mp_iallgather_c13, mp_iallgather_z13, &
         mp_iallgather_i22, mp_iallgather_l22, mp_iallgather_r22, mp_iallgather_d22, &
         mp_iallgather_c22, mp_iallgather_z22, mp_iallgather_i24, mp_iallgather_l24, &
         mp_iallgather_r24, mp_iallgather_d24, mp_iallgather_c24, mp_iallgather_z24, &
         mp_iallgather_i33, mp_iallgather_l33, mp_iallgather_r33, mp_iallgather_d33, &
         mp_iallgather_c33, mp_iallgather_z33
      GENERIC, PUBLIC :: iallgather => mp_iallgather_i, mp_iallgather_l, &
         mp_iallgather_r, mp_iallgather_d, mp_iallgather_c, mp_iallgather_z, &
         mp_iallgather_i11, mp_iallgather_l11, mp_iallgather_r11, mp_iallgather_d11, &
         mp_iallgather_c11, mp_iallgather_z11, mp_iallgather_i13, mp_iallgather_l13, &
         mp_iallgather_r13, mp_iallgather_d13, mp_iallgather_c13, mp_iallgather_z13, &
         mp_iallgather_i22, mp_iallgather_l22, mp_iallgather_r22, mp_iallgather_d22, &
         mp_iallgather_c22, mp_iallgather_z22, mp_iallgather_i24, mp_iallgather_l24, &
         mp_iallgather_r24, mp_iallgather_d24, mp_iallgather_c24, mp_iallgather_z24, &
         mp_iallgather_i33, mp_iallgather_l33, mp_iallgather_r33, mp_iallgather_d33, &
         mp_iallgather_c33, mp_iallgather_z33

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_iallgatherv_iv, mp_iallgatherv_iv2, &
         mp_iallgatherv_lv, mp_iallgatherv_lv2, mp_iallgatherv_rv, mp_iallgatherv_rv2, &
         mp_iallgatherv_dv, mp_iallgatherv_dv2, mp_iallgatherv_cv, mp_iallgatherv_cv2, &
         mp_iallgatherv_zv, mp_iallgatherv_zv2
      GENERIC, PUBLIC :: iallgatherv => mp_iallgatherv_iv, mp_iallgatherv_iv2, &
         mp_iallgatherv_lv, mp_iallgatherv_lv2, mp_iallgatherv_rv, mp_iallgatherv_rv2, &
         mp_iallgatherv_dv, mp_iallgatherv_dv2, mp_iallgatherv_cv, mp_iallgatherv_cv2, &
         mp_iallgatherv_zv, mp_iallgatherv_zv2

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_scatter_iv, mp_scatter_lv, &
         mp_scatter_rv, mp_scatter_dv, mp_scatter_cv, mp_scatter_zv
      GENERIC, PUBLIC :: scatter => mp_scatter_iv, mp_scatter_lv, &
         mp_scatter_rv, mp_scatter_dv, mp_scatter_cv, mp_scatter_zv

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_iscatter_i, mp_iscatter_l, &
         mp_iscatter_r, mp_iscatter_d, mp_iscatter_c, mp_iscatter_z, &
         mp_iscatter_iv2, mp_iscatter_lv2, mp_iscatter_rv2, mp_iscatter_dv2, &
         mp_iscatter_cv2, mp_iscatter_zv2
      GENERIC, PUBLIC :: iscatter => mp_iscatter_i, mp_iscatter_l, &
         mp_iscatter_r, mp_iscatter_d, mp_iscatter_c, mp_iscatter_z, &
         mp_iscatter_iv2, mp_iscatter_lv2, mp_iscatter_rv2, mp_iscatter_dv2, &
         mp_iscatter_cv2, mp_iscatter_zv2

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_iscatterv_iv, mp_iscatterv_lv, &
         mp_iscatterv_rv, mp_iscatterv_dv, mp_iscatterv_cv, mp_iscatterv_zv
      GENERIC, PUBLIC :: iscatterv => mp_iscatterv_iv, mp_iscatterv_lv, &
         mp_iscatterv_rv, mp_iscatterv_dv, mp_iscatterv_cv, mp_iscatterv_zv

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_alltoall_i, mp_alltoall_i22, mp_alltoall_i33, &
         mp_alltoall_i44, mp_alltoall_i55, mp_alltoall_i45, mp_alltoall_i34, &
         mp_alltoall_i11v, mp_alltoall_i22v, mp_alltoall_i54, &
         mp_alltoall_l, mp_alltoall_l22, mp_alltoall_l33, &
         mp_alltoall_l44, mp_alltoall_l55, mp_alltoall_l45, mp_alltoall_l34, &
         mp_alltoall_l11v, mp_alltoall_l22v, mp_alltoall_l54, &
         mp_alltoall_r, mp_alltoall_r22, mp_alltoall_r33, &
         mp_alltoall_r44, mp_alltoall_r55, mp_alltoall_r45, mp_alltoall_r34, &
         mp_alltoall_r11v, mp_alltoall_r22v, mp_alltoall_r54, &
         mp_alltoall_d, mp_alltoall_d22, mp_alltoall_d33, &
         mp_alltoall_d44, mp_alltoall_d55, mp_alltoall_d45, mp_alltoall_d34, &
         mp_alltoall_d11v, mp_alltoall_d22v, mp_alltoall_d54, &
         mp_alltoall_c, mp_alltoall_c22, mp_alltoall_c33, &
         mp_alltoall_c44, mp_alltoall_c55, mp_alltoall_c45, mp_alltoall_c34, &
         mp_alltoall_c11v, mp_alltoall_c22v, mp_alltoall_c54, &
         mp_alltoall_z, mp_alltoall_z22, mp_alltoall_z33, &
         mp_alltoall_z44, mp_alltoall_z55, mp_alltoall_z45, mp_alltoall_z34, &
         mp_alltoall_z11v, mp_alltoall_z22v, mp_alltoall_z54
      GENERIC, PUBLIC :: alltoall => mp_alltoall_i, mp_alltoall_i22, mp_alltoall_i33, &
         mp_alltoall_i44, mp_alltoall_i55, mp_alltoall_i45, mp_alltoall_i34, &
         mp_alltoall_i11v, mp_alltoall_i22v, mp_alltoall_i54, &
         mp_alltoall_l, mp_alltoall_l22, mp_alltoall_l33, &
         mp_alltoall_l44, mp_alltoall_l55, mp_alltoall_l45, mp_alltoall_l34, &
         mp_alltoall_l11v, mp_alltoall_l22v, mp_alltoall_l54, &
         mp_alltoall_r, mp_alltoall_r22, mp_alltoall_r33, &
         mp_alltoall_r44, mp_alltoall_r55, mp_alltoall_r45, mp_alltoall_r34, &
         mp_alltoall_r11v, mp_alltoall_r22v, mp_alltoall_r54, &
         mp_alltoall_d, mp_alltoall_d22, mp_alltoall_d33, &
         mp_alltoall_d44, mp_alltoall_d55, mp_alltoall_d45, mp_alltoall_d34, &
         mp_alltoall_d11v, mp_alltoall_d22v, mp_alltoall_d54, &
         mp_alltoall_c, mp_alltoall_c22, mp_alltoall_c33, &
         mp_alltoall_c44, mp_alltoall_c55, mp_alltoall_c45, mp_alltoall_c34, &
         mp_alltoall_c11v, mp_alltoall_c22v, mp_alltoall_c54, &
         mp_alltoall_z, mp_alltoall_z22, mp_alltoall_z33, &
         mp_alltoall_z44, mp_alltoall_z55, mp_alltoall_z45, mp_alltoall_z34, &
         mp_alltoall_z11v, mp_alltoall_z22v, mp_alltoall_z54

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_send_i, mp_send_iv, mp_send_im2, mp_send_im3, &
         mp_send_l, mp_send_lv, mp_send_lm2, mp_send_lm3, &
         mp_send_r, mp_send_rv, mp_send_rm2, mp_send_rm3, &
         mp_send_d, mp_send_dv, mp_send_dm2, mp_send_dm3, &
         mp_send_c, mp_send_cv, mp_send_cm2, mp_send_cm3, &
         mp_send_z, mp_send_zv, mp_send_zm2, mp_send_zm3
      GENERIC, PUBLIC :: send => mp_send_i, mp_send_iv, mp_send_im2, mp_send_im3, &
         mp_send_l, mp_send_lv, mp_send_lm2, mp_send_lm3, &
         mp_send_r, mp_send_rv, mp_send_rm2, mp_send_rm3, &
         mp_send_d, mp_send_dv, mp_send_dm2, mp_send_dm3, &
         mp_send_c, mp_send_cv, mp_send_cm2, mp_send_cm3, &
         mp_send_z, mp_send_zv, mp_send_zm2, mp_send_zm3

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_recv_i, mp_recv_iv, mp_recv_im2, mp_recv_im3, &
         mp_recv_l, mp_recv_lv, mp_recv_lm2, mp_recv_lm3, &
         mp_recv_r, mp_recv_rv, mp_recv_rm2, mp_recv_rm3, &
         mp_recv_d, mp_recv_dv, mp_recv_dm2, mp_recv_dm3, &
         mp_recv_c, mp_recv_cv, mp_recv_cm2, mp_recv_cm3, &
         mp_recv_z, mp_recv_zv, mp_recv_zm2, mp_recv_zm3
      GENERIC, PUBLIC :: recv => mp_recv_i, mp_recv_iv, mp_recv_im2, mp_recv_im3, &
         mp_recv_l, mp_recv_lv, mp_recv_lm2, mp_recv_lm3, &
         mp_recv_r, mp_recv_rv, mp_recv_rm2, mp_recv_rm3, &
         mp_recv_d, mp_recv_dv, mp_recv_dm2, mp_recv_dm3, &
         mp_recv_c, mp_recv_cv, mp_recv_cm2, mp_recv_cm3, &
         mp_recv_z, mp_recv_zv, mp_recv_zm2, mp_recv_zm3

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_isendrecv_i, mp_isendrecv_iv, &
         mp_isendrecv_l, mp_isendrecv_lv, mp_isendrecv_r, mp_isendrecv_rv, &
         mp_isendrecv_d, mp_isendrecv_dv, mp_isendrecv_c, mp_isendrecv_cv, &
         mp_isendrecv_z, mp_isendrecv_zv
      GENERIC, PUBLIC :: isendrecv => mp_isendrecv_i, mp_isendrecv_iv, &
         mp_isendrecv_l, mp_isendrecv_lv, mp_isendrecv_r, mp_isendrecv_rv, &
         mp_isendrecv_d, mp_isendrecv_dv, mp_isendrecv_c, mp_isendrecv_cv, &
         mp_isendrecv_z, mp_isendrecv_zv

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_isend_iv, mp_isend_im2, mp_isend_im3, mp_isend_im4, &
         mp_isend_lv, mp_isend_lm2, mp_isend_lm3, mp_isend_lm4, &
         mp_isend_rv, mp_isend_rm2, mp_isend_rm3, mp_isend_rm4, &
         mp_isend_dv, mp_isend_dm2, mp_isend_dm3, mp_isend_dm4, &
         mp_isend_cv, mp_isend_cm2, mp_isend_cm3, mp_isend_cm4, &
         mp_isend_zv, mp_isend_zm2, mp_isend_zm3, mp_isend_zm4, &
         mp_isend_bv, mp_isend_bm3, mp_isend_custom
      GENERIC, PUBLIC :: isend => mp_isend_iv, mp_isend_im2, mp_isend_im3, mp_isend_im4, &
         mp_isend_lv, mp_isend_lm2, mp_isend_lm3, mp_isend_lm4, &
         mp_isend_rv, mp_isend_rm2, mp_isend_rm3, mp_isend_rm4, &
         mp_isend_dv, mp_isend_dm2, mp_isend_dm3, mp_isend_dm4, &
         mp_isend_cv, mp_isend_cm2, mp_isend_cm3, mp_isend_cm4, &
         mp_isend_zv, mp_isend_zm2, mp_isend_zm3, mp_isend_zm4, &
         mp_isend_bv, mp_isend_bm3, mp_isend_custom

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_irecv_iv, mp_irecv_im2, mp_irecv_im3, mp_irecv_im4, &
         mp_irecv_lv, mp_irecv_lm2, mp_irecv_lm3, mp_irecv_lm4, &
         mp_irecv_rv, mp_irecv_rm2, mp_irecv_rm3, mp_irecv_rm4, &
         mp_irecv_dv, mp_irecv_dm2, mp_irecv_dm3, mp_irecv_dm4, &
         mp_irecv_cv, mp_irecv_cm2, mp_irecv_cm3, mp_irecv_cm4, &
         mp_irecv_zv, mp_irecv_zm2, mp_irecv_zm3, mp_irecv_zm4, &
         mp_irecv_bv, mp_irecv_bm3, mp_irecv_custom
      GENERIC, PUBLIC :: irecv => mp_irecv_iv, mp_irecv_im2, mp_irecv_im3, mp_irecv_im4, &
         mp_irecv_lv, mp_irecv_lm2, mp_irecv_lm3, mp_irecv_lm4, &
         mp_irecv_rv, mp_irecv_rm2, mp_irecv_rm3, mp_irecv_rm4, &
         mp_irecv_dv, mp_irecv_dm2, mp_irecv_dm3, mp_irecv_dm4, &
         mp_irecv_cv, mp_irecv_cm2, mp_irecv_cm3, mp_irecv_cm4, &
         mp_irecv_zv, mp_irecv_zm2, mp_irecv_zm3, mp_irecv_zm4, &
         mp_irecv_bv, mp_irecv_bm3, mp_irecv_custom

      PROCEDURE, PUBLIC, PASS(comm), NON_OVERRIDABLE :: probe => mp_probe

      PROCEDURE, PUBLIC, PASS(comm), NON_OVERRIDABLE :: sync => mp_sync
      PROCEDURE, PUBLIC, PASS(comm), NON_OVERRIDABLE :: isync => mp_isync

      PROCEDURE, PUBLIC, PASS(comm1), NON_OVERRIDABLE :: compare => mp_comm_compare
      PROCEDURE, PUBLIC, PASS(comm1), NON_OVERRIDABLE :: rank_compare => mp_rank_compare

      PROCEDURE, PUBLIC, PASS(comm2), NON_OVERRIDABLE :: dup => mp_comm_dup
      PROCEDURE, PUBLIC, PASS(comm), NON_OVERRIDABLE :: free => mp_comm_free

      PROCEDURE, PUBLIC, PASS(comm), NON_OVERRIDABLE :: environ => mp_environ_l

      PROCEDURE, PRIVATE, PASS(sub_comm), NON_OVERRIDABLE :: mp_comm_split, mp_comm_split_direct
      GENERIC, PUBLIC :: from_split => mp_comm_split, mp_comm_split_direct
      PROCEDURE, PUBLIC, PASS(mp_new_comm), NON_OVERRIDABLE :: from_reordering => mp_reordering
      PROCEDURE, PRIVATE, PASS(comm_new), NON_OVERRIDABLE :: mp_comm_copy
      GENERIC, PUBLIC :: ASSIGNMENT(=) => mp_comm_copy
   END TYPE

   TYPE mp_request_type
      PRIVATE
      MPI_REQUEST_TYPE :: handle = mp_request_null_handle
   CONTAINS
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: set_handle => mp_request_type_set_handle
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: get_handle => mp_request_type_get_handle
      PROCEDURE, PRIVATE, NON_OVERRIDABLE :: mp_request_op_eq
      PROCEDURE, PRIVATE, NON_OVERRIDABLE :: mp_request_op_neq
      GENERIC, PUBLIC :: OPERATOR(.EQ.) => mp_request_op_eq
      GENERIC, PUBLIC :: OPERATOR(.NE.) => mp_request_op_neq

      PROCEDURE, PUBLIC, PASS(request), NON_OVERRIDABLE :: test => mp_test_1

      PROCEDURE, PUBLIC, PASS(request), NON_OVERRIDABLE :: wait => mp_wait
   END TYPE

   TYPE mp_win_type
      PRIVATE
      MPI_WIN_TYPE :: handle = mp_win_null_handle
   CONTAINS
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: set_handle => mp_win_type_set_handle
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: get_handle => mp_win_type_get_handle
      PROCEDURE, PRIVATE, NON_OVERRIDABLE :: mp_win_op_eq
      PROCEDURE, PRIVATE, NON_OVERRIDABLE :: mp_win_op_neq
      GENERIC, PUBLIC :: OPERATOR(.EQ.) => mp_win_op_eq
      GENERIC, PUBLIC :: OPERATOR(.NE.) => mp_win_op_neq

      PROCEDURE, PRIVATE, PASS(win), NON_OVERRIDABLE :: mp_win_create_iv, mp_win_create_lv, &
         mp_win_create_rv, mp_win_create_dv, mp_win_create_cv, mp_win_create_zv
      GENERIC, PUBLIC :: create => mp_win_create_iv, mp_win_create_lv, &
         mp_win_create_rv, mp_win_create_dv, mp_win_create_cv, mp_win_create_zv

      PROCEDURE, PRIVATE, PASS(win), NON_OVERRIDABLE :: mp_rget_iv, mp_rget_lv, &
         mp_rget_rv, mp_rget_dv, mp_rget_cv, mp_rget_zv
      GENERIC, PUBLIC :: rget => mp_rget_iv, mp_rget_lv, &
         mp_rget_rv, mp_rget_dv, mp_rget_cv, mp_rget_zv

      PROCEDURE, PUBLIC, PASS(win), NON_OVERRIDABLE :: free => mp_win_free
      PROCEDURE, PUBLIC, PASS(win_new), NON_OVERRIDABLE :: mp_win_assign
      GENERIC, PUBLIC :: ASSIGNMENT(=) => mp_win_assign

      PROCEDURE, PUBLIC, PASS(win), NON_OVERRIDABLE :: lock_all => mp_win_lock_all
      PROCEDURE, PUBLIC, PASS(win), NON_OVERRIDABLE :: unlock_all => mp_win_unlock_all
      PROCEDURE, PUBLIC, PASS(win), NON_OVERRIDABLE :: flush_all => mp_win_flush_all
   END TYPE

   TYPE mp_file_type
      PRIVATE
      MPI_FILE_TYPE :: handle = mp_file_null_handle
   CONTAINS
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: set_handle => mp_file_type_set_handle
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: get_handle => mp_file_type_get_handle
      PROCEDURE, PRIVATE, NON_OVERRIDABLE :: mp_file_op_eq
      PROCEDURE, PRIVATE, NON_OVERRIDABLE :: mp_file_op_neq
      GENERIC, PUBLIC :: OPERATOR(.EQ.) => mp_file_op_eq
      GENERIC, PUBLIC :: OPERATOR(.NE.) => mp_file_op_neq

      PROCEDURE, PRIVATE, PASS(fh), NON_OVERRIDABLE :: mp_file_write_at_ch, mp_file_write_at_chv, &
         mp_file_write_at_i, mp_file_write_at_iv, mp_file_write_at_r, mp_file_write_at_rv, &
         mp_file_write_at_d, mp_file_write_at_dv, mp_file_write_at_c, mp_file_write_at_cv, &
         mp_file_write_at_z, mp_file_write_at_zv, mp_file_write_at_l, mp_file_write_at_lv
      GENERIC, PUBLIC :: write_at => mp_file_write_at_ch, mp_file_write_at_chv, &
         mp_file_write_at_i, mp_file_write_at_iv, mp_file_write_at_r, mp_file_write_at_rv, &
         mp_file_write_at_d, mp_file_write_at_dv, mp_file_write_at_c, mp_file_write_at_cv, &
         mp_file_write_at_z, mp_file_write_at_zv, mp_file_write_at_l, mp_file_write_at_lv

      PROCEDURE, PRIVATE, PASS(fh), NON_OVERRIDABLE :: mp_file_write_at_all_ch, mp_file_write_at_all_chv, &
         mp_file_write_at_all_i, mp_file_write_at_all_iv, mp_file_write_at_all_l, mp_file_write_at_all_lv, &
         mp_file_write_at_all_r, mp_file_write_at_all_rv, mp_file_write_at_all_d, mp_file_write_at_all_dv, &
         mp_file_write_at_all_c, mp_file_write_at_all_cv, mp_file_write_at_all_z, mp_file_write_at_all_zv
      GENERIC, PUBLIC :: write_at_all => mp_file_write_at_all_ch, mp_file_write_at_all_chv, &
         mp_file_write_at_all_i, mp_file_write_at_all_iv, mp_file_write_at_all_l, mp_file_write_at_all_lv, &
         mp_file_write_at_all_r, mp_file_write_at_all_rv, mp_file_write_at_all_d, mp_file_write_at_all_dv, &
         mp_file_write_at_all_c, mp_file_write_at_all_cv, mp_file_write_at_all_z, mp_file_write_at_all_zv

      PROCEDURE, PRIVATE, PASS(fh), NON_OVERRIDABLE :: mp_file_read_at_ch, mp_file_read_at_chv, &
         mp_file_read_at_i, mp_file_read_at_iv, mp_file_read_at_r, mp_file_read_at_rv, &
         mp_file_read_at_d, mp_file_read_at_dv, mp_file_read_at_c, mp_file_read_at_cv, &
         mp_file_read_at_z, mp_file_read_at_zv, mp_file_read_at_l, mp_file_read_at_lv
      GENERIC, PUBLIC :: read_at => mp_file_read_at_ch, mp_file_read_at_chv, &
         mp_file_read_at_i, mp_file_read_at_iv, mp_file_read_at_r, mp_file_read_at_rv, &
         mp_file_read_at_d, mp_file_read_at_dv, mp_file_read_at_c, mp_file_read_at_cv, &
         mp_file_read_at_z, mp_file_read_at_zv, mp_file_read_at_l, mp_file_read_at_lv

      PROCEDURE, PRIVATE, PASS(fh), NON_OVERRIDABLE :: mp_file_read_at_all_ch, mp_file_read_at_all_chv, &
         mp_file_read_at_all_i, mp_file_read_at_all_iv, mp_file_read_at_all_l, mp_file_read_at_all_lv, &
         mp_file_read_at_all_r, mp_file_read_at_all_rv, mp_file_read_at_all_d, mp_file_read_at_all_dv, &
         mp_file_read_at_all_c, mp_file_read_at_all_cv, mp_file_read_at_all_z, mp_file_read_at_all_zv
      GENERIC, PUBLIC :: read_at_all => mp_file_read_at_all_ch, mp_file_read_at_all_chv, &
         mp_file_read_at_all_i, mp_file_read_at_all_iv, mp_file_read_at_all_l, mp_file_read_at_all_lv, &
         mp_file_read_at_all_r, mp_file_read_at_all_rv, mp_file_read_at_all_d, mp_file_read_at_all_dv, &
         mp_file_read_at_all_c, mp_file_read_at_all_cv, mp_file_read_at_all_z, mp_file_read_at_all_zv

      PROCEDURE, PUBLIC, PASS(fh), NON_OVERRIDABLE :: open => mp_file_open
      PROCEDURE, PUBLIC, PASS(fh), NON_OVERRIDABLE :: close => mp_file_close
      PROCEDURE, PRIVATE, PASS(fh_new), NON_OVERRIDABLE :: mp_file_assign
      GENERIC, PUBLIC :: ASSIGNMENT(=) => mp_file_assign

      PROCEDURE, PUBLIC, PASS(fh), NON_OVERRIDABLE :: get_size => mp_file_get_size
      PROCEDURE, PUBLIC, PASS(fh), NON_OVERRIDABLE :: get_position => mp_file_get_position

      PROCEDURE, PUBLIC, PASS(fh), NON_OVERRIDABLE :: read_all => mp_file_read_all_chv
      PROCEDURE, PUBLIC, PASS(fh), NON_OVERRIDABLE :: write_all => mp_file_write_all_chv
   END TYPE

   TYPE mp_info_type
      PRIVATE
      MPI_INFO_TYPE :: handle = mp_info_null_handle
   CONTAINS
      PROCEDURE, NON_OVERRIDABLE :: set_handle => mp_info_type_set_handle
      PROCEDURE, NON_OVERRIDABLE :: get_handle => mp_info_type_get_handle
      PROCEDURE, PRIVATE, NON_OVERRIDABLE :: mp_info_op_eq
      PROCEDURE, PRIVATE, NON_OVERRIDABLE :: mp_info_op_neq
      GENERIC, PUBLIC :: OPERATOR(.EQ.) => mp_info_op_eq
      GENERIC, PUBLIC :: OPERATOR(.NE.) => mp_info_op_neq
   END TYPE

   TYPE, EXTENDS(mp_comm_type) :: mp_cart_type
   CONTAINS
      PROCEDURE, PUBLIC, PASS(comm_cart), NON_OVERRIDABLE :: create => mp_cart_create
      PROCEDURE, PUBLIC, PASS(sub_comm), NON_OVERRIDABLE :: from_sub => mp_cart_sub

      PROCEDURE, PRIVATE, PASS(comm), NON_OVERRIDABLE :: mp_environ_c, mp_environ_c2
      GENERIC, PUBLIC :: environ_cart => mp_environ_c, mp_environ_c2

      PROCEDURE, PUBLIC, PASS(comm), NON_OVERRIDABLE :: coords => mp_cart_coords
      PROCEDURE, PUBLIC, PASS(comm), NON_OVERRIDABLE :: rank_cart => mp_cart_rank
   END TYPE

   ! Create the constants from the corresponding handles
   TYPE(mp_comm_type), PARAMETER, PUBLIC :: mp_comm_null = mp_comm_type(mp_comm_null_handle)
   TYPE(mp_comm_type), PARAMETER, PUBLIC :: mp_comm_self = mp_comm_type(mp_comm_self_handle)
   TYPE(mp_comm_type), PARAMETER, PUBLIC :: mp_comm_world = mp_comm_type(mp_comm_world_handle)
   TYPE(mp_request_type), PARAMETER, PUBLIC :: mp_request_null = mp_request_type(mp_request_null_handle)
   TYPE(mp_win_type), PARAMETER, PUBLIC :: mp_win_null = mp_win_type(mp_win_null_handle)
   TYPE(mp_file_type), PARAMETER, PUBLIC :: mp_file_null = mp_file_type(mp_file_null_handle)
   TYPE(mp_info_type), PARAMETER, PUBLIC :: mp_info_null = mp_info_type(mp_info_null_handle)

#if !defined(__parallel)
   ! This communicator is to be used in serial mode to emulate a valid communicator which is not a compiler constant
   INTEGER, PARAMETER, PRIVATE :: mp_comm_default_handle = 1
   TYPE(mp_comm_type), PARAMETER, PRIVATE :: mp_comm_default = mp_comm_type(mp_comm_default_handle)
#endif

   ! Constants to compare communicators
   INTEGER, PARAMETER, PUBLIC :: mp_comm_ident = 0
   INTEGER, PARAMETER, PUBLIC :: mp_comm_congruent = 1
   INTEGER, PARAMETER, PUBLIC :: mp_comm_similar = 2
   INTEGER, PARAMETER, PUBLIC :: mp_comm_unequal = 3
   INTEGER, PARAMETER, PUBLIC :: mp_comm_compare_default = -1

   ! init and error
   PUBLIC :: mp_world_init, mp_world_finalize
   PUBLIC :: mp_abort

   ! performance gathering
   PUBLIC :: mp_perf_env_type
   PUBLIC :: mp_perf_env_retain, mp_perf_env_release
   PUBLIC :: add_mp_perf_env, rm_mp_perf_env, get_mp_perf_env, describe_mp_perf_env

   ! informational / generation of sub comms
   PUBLIC :: mp_dims_create
   PUBLIC :: cp2k_is_parallel
   PUBLIC :: mp_get_node_global_rank

   ! message passing
   PUBLIC :: mp_waitall, mp_waitany
   PUBLIC :: mp_testall, mp_testany

   ! Memory management
   PUBLIC :: mp_allocate, mp_deallocate

   ! I/O
   PUBLIC :: mp_file_delete
   PUBLIC :: mp_file_get_amode

   ! some 'advanced types' currently only used for dbcsr
   PUBLIC :: mp_type_descriptor_type
   PUBLIC :: mp_type_make
   PUBLIC :: mp_type_size

   ! some benchmarking code
   PUBLIC :: mpi_perf_test

   ! vector types
   PUBLIC :: mp_type_indexed_make_r, mp_type_indexed_make_d, &
             mp_type_indexed_make_c, mp_type_indexed_make_z

   ! More I/O types and routines: variable spaced data using bytes for spacings
   PUBLIC :: mp_file_descriptor_type
   PUBLIC :: mp_file_type_free
   PUBLIC :: mp_file_type_hindexed_make_chv
   PUBLIC :: mp_file_type_set_view_chv

   PUBLIC :: mp_get_library_version

   ! assumed to be private

   INTERFACE mp_waitall
      MODULE PROCEDURE mp_waitall_1, mp_waitall_2
   END INTERFACE

   INTERFACE mp_testall
      MODULE PROCEDURE mp_testall_tv
   END INTERFACE

   INTERFACE mp_testany
      MODULE PROCEDURE mp_testany_1, mp_testany_2
   END INTERFACE

   INTERFACE mp_type_free
      MODULE PROCEDURE mp_type_free_m, mp_type_free_v
   END INTERFACE

   !
   ! interfaces to deal easily with scalars / vectors / matrices / ...
   ! of the different types (integers, doubles, logicals, characters)
   !
   INTERFACE mp_allocate
      MODULE PROCEDURE mp_allocate_i, &
         mp_allocate_l, &
         mp_allocate_r, &
         mp_allocate_d, &
         mp_allocate_c, &
         mp_allocate_z
   END INTERFACE

   INTERFACE mp_deallocate
      MODULE PROCEDURE mp_deallocate_i, &
         mp_deallocate_l, &
         mp_deallocate_r, &
         mp_deallocate_d, &
         mp_deallocate_c, &
         mp_deallocate_z
   END INTERFACE

   INTERFACE mp_type_make
      MODULE PROCEDURE mp_type_make_struct
      MODULE PROCEDURE mp_type_make_i, mp_type_make_l, &
         mp_type_make_r, mp_type_make_d, &
         mp_type_make_c, mp_type_make_z
   END INTERFACE

   INTERFACE mp_alloc_mem
      MODULE PROCEDURE mp_alloc_mem_i, mp_alloc_mem_l, &
         mp_alloc_mem_d, mp_alloc_mem_z, &
         mp_alloc_mem_r, mp_alloc_mem_c
   END INTERFACE

   INTERFACE mp_free_mem
      MODULE PROCEDURE mp_free_mem_i, mp_free_mem_l, &
         mp_free_mem_d, mp_free_mem_z, &
         mp_free_mem_r, mp_free_mem_c
   END INTERFACE

! Type declarations
   TYPE mp_indexing_meta_type
      INTEGER, DIMENSION(:), POINTER :: index, chunks
   END TYPE mp_indexing_meta_type

   TYPE mp_type_descriptor_type
      MPI_DATA_TYPE :: type_handle
      INTEGER :: length
#if defined(__parallel)
      INTEGER(kind=mpi_address_kind) :: base
#endif
      INTEGER(kind=int_4), DIMENSION(:), POINTER :: data_i
      INTEGER(kind=int_8), DIMENSION(:), POINTER :: data_l
      REAL(kind=real_4), DIMENSION(:), POINTER :: data_r
      REAL(kind=real_8), DIMENSION(:), POINTER :: data_d
      COMPLEX(kind=real_4), DIMENSION(:), POINTER :: data_c
      COMPLEX(kind=real_8), DIMENSION(:), POINTER :: data_z
      TYPE(mp_type_descriptor_type), DIMENSION(:), POINTER :: subtype
      INTEGER :: vector_descriptor(2)
      LOGICAL :: has_indexing
      TYPE(mp_indexing_meta_type) :: index_descriptor
   END TYPE mp_type_descriptor_type

   TYPE mp_file_indexing_meta_type
      INTEGER, DIMENSION(:), POINTER   :: index
      INTEGER(kind=file_offset), &
         DIMENSION(:), POINTER         :: chunks
   END TYPE mp_file_indexing_meta_type

   TYPE mp_file_descriptor_type
      MPI_DATA_TYPE :: type_handle
      INTEGER                          :: length
      LOGICAL                          :: has_indexing = .FALSE.
      TYPE(mp_file_indexing_meta_type) :: index_descriptor
   END TYPE

   ! type internally used to store message passing performance indicators
! **************************************************************************************************
   TYPE mp_perf_type
      CHARACTER(LEN=20) :: name
      INTEGER :: count
      REAL(KIND=dp) :: msg_size
   END TYPE mp_perf_type

   INTEGER, PARAMETER :: MAX_PERF = 28

! **************************************************************************************************
   TYPE mp_perf_env_type
      !private
      INTEGER :: ref_count
      TYPE(mp_perf_type), DIMENSION(MAX_PERF) :: mp_perfs
   END TYPE mp_perf_env_type

! **************************************************************************************************
   TYPE mp_perf_env_p_type
      TYPE(mp_perf_env_type), POINTER         :: mp_perf_env => Null()
   END TYPE mp_perf_env_p_type

   ! introduce a stack of mp_perfs, first index is the stack pointer, for convenience is replacing
   INTEGER, PARAMETER :: max_stack_size = 10
   INTEGER            :: stack_pointer = 0
   ! target attribute needed as a hack around ifc 7.1 bug
   TYPE(mp_perf_env_p_type), DIMENSION(max_stack_size), TARGET, SAVE :: mp_perf_stack

   CHARACTER(LEN=20), PARAMETER :: sname(MAX_PERF) = &
                                   (/"MP_Group            ", "MP_Bcast            ", "MP_Allreduce        ", &
                                     "MP_Gather           ", "MP_Sync             ", "MP_Alltoall         ", &
                                     "MP_SendRecv         ", "MP_ISendRecv        ", "MP_Wait             ", &
                                     "MP_comm_split       ", "MP_ISend            ", "MP_IRecv            ", &
                                     "MP_Send             ", "MP_Recv             ", "MP_Memory           ", &
                                     "MP_Put              ", "MP_Get              ", "MP_Fence            ", &
                                     "MP_Win_Lock         ", "MP_Win_Create       ", "MP_Win_Free         ", &
                                     "MP_IBcast           ", "MP_IAllreduce       ", "MP_IScatter         ", &
                                     "MP_RGet             ", "MP_Isync            ", "MP_Read_All         ", &
                                     "MP_Write_All        "/)

   ! we make some assumptions on the length of INTEGERS, REALS and LOGICALS
   INTEGER, PARAMETER :: intlen = BIT_SIZE(0)/8
   INTEGER, PARAMETER :: reallen = 8
   INTEGER, PARAMETER :: loglen = BIT_SIZE(0)/8
   INTEGER, PARAMETER :: charlen = 1

   LOGICAL, PUBLIC, SAVE :: mp_collect_timings = .FALSE.

CONTAINS

   #:mute
      #:set types = ["comm", "request", "win", "file", "info"]
   #:endmute
   #:for type in types
      LOGICAL FUNCTION mp_${type}$_op_eq(${type}$1, ${type}$2)
         CLASS(mp_${type}$_type), INTENT(IN) :: ${type}$1, ${type}$2
         mp_${type}$_op_eq = (${type}$1%handle .EQ. ${type}$2%handle)
      END FUNCTION mp_${type}$_op_eq

      LOGICAL FUNCTION mp_${type}$_op_neq(${type}$1, ${type}$2)
         CLASS(mp_${type}$_type), INTENT(IN) :: ${type}$1, ${type}$2
         mp_${type}$_op_neq = (${type}$1%handle .NE. ${type}$2%handle)
      END FUNCTION mp_${type}$_op_neq

      ELEMENTAL SUBROUTINE mp_${type}$_type_set_handle(this, handle)
         CLASS(mp_${type}$_type), INTENT(INOUT) :: this
         INTEGER, INTENT(IN) :: handle

#if defined(__parallel) && defined(__MPI_F08)
         this%handle%mpi_val = handle
#else
         this%handle = handle
#endif
      END SUBROUTINE mp_${type}$_type_set_handle

      ELEMENTAL FUNCTION mp_${type}$_type_get_handle(this) RESULT(handle)
         CLASS(mp_${type}$_type), INTENT(IN) :: this
         INTEGER :: handle

#if defined(__parallel) && defined(__MPI_F08)
         handle = this%handle%mpi_val
#else
         handle = this%handle
#endif
      END FUNCTION mp_${type}$_type_get_handle
   #:endfor

! **************************************************************************************************
!> \brief initializes the system default communicator
!> \param mp_comm [output] : handle of the default communicator
!> \par History
!>      2.2004 created [Joost VandeVondele ]
!> \note
!>      should only be called once
! **************************************************************************************************
   SUBROUTINE mp_world_init(mp_comm)
      TYPE(mp_comm_type), INTENT(OUT)                     :: mp_comm
#if defined(__parallel)
      INTEGER                                  :: ierr
!$    INTEGER                                  :: provided_tsl
!$    LOGICAL                                  :: no_threading_support

#if defined(__NO_MPI_THREAD_SUPPORT_CHECK)
      ! Hack that does not request or check MPI thread support level.
      ! User asserts that the MPI library will work correctly with
      ! threads.
!
!$    no_threading_support = .TRUE.
#else
      ! Does the right thing when using OpenMP: requests that the MPI
      ! library supports serialized mode and verifies that the MPI library
      ! provides that support.
      !
      ! Developers: Only the master thread will ever make calls to the
      ! MPI library.
!
!$    no_threading_support = .FALSE.
#endif
!$    IF (no_threading_support) THEN
         CALL mpi_init(ierr)
         IF (ierr /= 0) CALL mp_stop(ierr, "mpi_init @ mp_world_init")
!$    ELSE
!$OMP MASTER
!$       CALL mpi_init_thread(MPI_THREAD_SERIALIZED, provided_tsl, ierr)
!$       IF (ierr /= 0) CALL mp_stop(ierr, "mpi_init_thread @ mp_world_init")
!$       IF (provided_tsl .LT. MPI_THREAD_SERIALIZED) THEN
!$          CALL mp_stop(0, "MPI library does not support the requested level of threading (MPI_THREAD_SERIALIZED).")
!$       END IF
!$OMP END MASTER
!$    END IF
      CALL mpi_comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_comm_set_errhandler @ mp_world_init")
#endif
      debug_comm_count = 1
      mp_comm = mp_comm_world
      CALL add_mp_perf_env()
   END SUBROUTINE mp_world_init

! **************************************************************************************************
!> \brief re-create the system default communicator with a different MPI
!>        rank order
!> \param mp_comm [output] : handle of the default communicator
!> \param mp_new_comm ...
!> \param ranks_order ...
!> \par History
!>      1.2012 created [ Christiane Pousa ]
!> \note
!>      should only be called once, at very beginning of CP2K run
! **************************************************************************************************
   SUBROUTINE mp_reordering(mp_comm, mp_new_comm, ranks_order)
      CLASS(mp_comm_type), INTENT(IN)                      :: mp_comm
      CLASS(mp_comm_type), INTENT(out)                     :: mp_new_comm
      INTEGER, DIMENSION(:), CONTIGUOUS, INTENT(IN)                    :: ranks_order

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_reordering'

      INTEGER                                  :: handle, ierr
#if defined(__parallel)
      MPI_GROUP_TYPE                                  :: newgroup, oldgroup
#endif

      CALL mp_timeset(routineN, handle)
      ierr = 0
#if defined(__parallel)

      CALL mpi_comm_group(mp_comm%handle, oldgroup, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_comm_group @ mp_reordering")
      CALL mpi_group_incl(oldgroup, SIZE(ranks_order), ranks_order, newgroup, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_group_incl @ mp_reordering")

      CALL mpi_comm_create(mp_comm%handle, newgroup, mp_new_comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_comm_create @ mp_reordering")

      CALL mpi_group_free(oldgroup, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_group_free @ mp_reordering")
      CALL mpi_group_free(newgroup, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_group_free @ mp_reordering")

      CALL add_perf(perf_id=1, count=1)
#else
      MARK_USED(mp_comm)
      MARK_USED(ranks_order)
      mp_new_comm%handle = mp_comm_default_handle
#endif
      debug_comm_count = debug_comm_count + 1
      CALL mp_timestop(handle)
   END SUBROUTINE mp_reordering

! **************************************************************************************************
!> \brief finalizes the system default communicator
!> \par History
!>      2.2004 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE mp_world_finalize()

      CHARACTER(LEN=default_string_length) :: debug_comm_count_char
#if defined(__parallel)
      INTEGER                                  :: ierr
      CALL mpi_barrier(MPI_COMM_WORLD, ierr) ! call mpi directly to avoid 0 stack pointer
#endif
      CALL rm_mp_perf_env()

      debug_comm_count = debug_comm_count - 1
#if defined(__parallel)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_barrier @ mp_world_finalize")
#endif
      IF (debug_comm_count .NE. 0) THEN
         ! A bug, we're leaking or double-freeing communicators. Needs to be fixed where the leak happens.
         ! Memory leak checking might be helpful to locate the culprit
         WRITE (unit=debug_comm_count_char, FMT='(I0)')
         CALL cp_abort(__LOCATION__, "mp_world_finalize: assert failed:"// &
                       "  leaking communicators "//TRIM(debug_comm_count_char))
      END IF
#if defined(__parallel)
      CALL mpi_finalize(ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_finalize @ mp_world_finalize")
#endif

   END SUBROUTINE mp_world_finalize

! all the following routines should work for a given communicator, not MPI_WORLD

! **************************************************************************************************
!> \brief start and stop the performance indicators
!>      for every call to start there has to be (exactly) one call to stop
!> \param perf_env ...
!> \par History
!>      2.2004 created [Joost VandeVondele]
!> \note
!>      can be used to measure performance of a sub-part of a program.
!>      timings measured here will not show up in the outer start/stops
!>      Doesn't need a fresh communicator
! **************************************************************************************************
   SUBROUTINE add_mp_perf_env(perf_env)
      TYPE(mp_perf_env_type), OPTIONAL, POINTER          :: perf_env

      stack_pointer = stack_pointer + 1
      IF (stack_pointer > max_stack_size) THEN
         CPABORT("stack_pointer too large : message_passing @ add_mp_perf_env")
      END IF
      NULLIFY (mp_perf_stack(stack_pointer)%mp_perf_env)
      IF (PRESENT(perf_env)) THEN
         mp_perf_stack(stack_pointer)%mp_perf_env => perf_env
         IF (ASSOCIATED(perf_env)) CALL mp_perf_env_retain(perf_env)
      END IF
      IF (.NOT. ASSOCIATED(mp_perf_stack(stack_pointer)%mp_perf_env)) THEN
         CALL mp_perf_env_create(mp_perf_stack(stack_pointer)%mp_perf_env)
      END IF
   END SUBROUTINE add_mp_perf_env

! **************************************************************************************************
!> \brief ...
!> \param perf_env ...
! **************************************************************************************************
   SUBROUTINE mp_perf_env_create(perf_env)
      TYPE(mp_perf_env_type), OPTIONAL, POINTER          :: perf_env

      INTEGER                                            :: i, stat

      NULLIFY (perf_env)
      ALLOCATE (perf_env, stat=stat)
      IF (stat /= 0) THEN
         CPABORT("allocation failed in mp_perf_env_create")
      END IF
      perf_env%ref_count = 1
      DO i = 1, MAX_PERF
         perf_env%mp_perfs(i)%name = sname(i)
         perf_env%mp_perfs(i)%count = 0
         perf_env%mp_perfs(i)%msg_size = 0.0_dp
      END DO

   END SUBROUTINE mp_perf_env_create

! **************************************************************************************************
!> \brief ...
!> \param perf_env ...
! **************************************************************************************************
   SUBROUTINE mp_perf_env_release(perf_env)
      TYPE(mp_perf_env_type), POINTER                    :: perf_env

      IF (ASSOCIATED(perf_env)) THEN
         IF (perf_env%ref_count < 1) THEN
            CPABORT("invalid ref_count: message_passing @ mp_perf_env_release")
         END IF
         perf_env%ref_count = perf_env%ref_count - 1
         IF (perf_env%ref_count == 0) THEN
            DEALLOCATE (perf_env)
         END IF
      END IF
      NULLIFY (perf_env)
   END SUBROUTINE mp_perf_env_release

! **************************************************************************************************
!> \brief ...
!> \param perf_env ...
! **************************************************************************************************
   SUBROUTINE mp_perf_env_retain(perf_env)
      TYPE(mp_perf_env_type), POINTER                    :: perf_env

      IF (.NOT. ASSOCIATED(perf_env)) THEN
         CPABORT("unassociated perf_env: message_passing @ mp_perf_env_retain")
      END IF
      IF (perf_env%ref_count < 1) THEN
         CPABORT("invalid ref_count: message_passing @ mp_perf_env_retain")
      END IF
      perf_env%ref_count = perf_env%ref_count + 1
   END SUBROUTINE mp_perf_env_retain

!.. reports the performance counters for the MPI run
! **************************************************************************************************
!> \brief ...
!> \param perf_env ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE mp_perf_env_describe(perf_env, iw)
      TYPE(mp_perf_env_type), POINTER          :: perf_env
      INTEGER, INTENT(IN)                      :: iw

#if defined(__parallel)
      INTEGER                                  :: i
      REAL(KIND=dp)                            :: vol
#endif

      IF (.NOT. ASSOCIATED(perf_env)) THEN
         CPABORT("unassociated perf_env : message_passing @ mp_perf_env_describe")
      END IF
      IF (perf_env%ref_count < 1) THEN
         CPABORT("invalid perf_env%ref_count : message_passing @ mp_perf_env_describe")
      END IF
#if defined(__parallel)
      IF (iw > 0) THEN
         WRITE (iw, '( /, 1X, 79("-") )')
         WRITE (iw, '( " -", 77X, "-" )')
         WRITE (iw, '( " -", 24X, A, 24X, "-" )') ' MESSAGE PASSING PERFORMANCE '
         WRITE (iw, '( " -", 77X, "-" )')
         WRITE (iw, '( 1X, 79("-"), / )')
         WRITE (iw, '( A, A, A )') ' ROUTINE', '             CALLS ', &
            '     AVE VOLUME [Bytes]'
         DO i = 1, MAX_PERF

            IF (perf_env%mp_perfs(i)%count > 0) THEN
               vol = perf_env%mp_perfs(i)%msg_size/REAL(perf_env%mp_perfs(i)%count, KIND=dp)
               IF (vol < 1.0_dp) THEN
                  WRITE (iw, '(1X,A15,T17,I10)') &
                     ADJUSTL(perf_env%mp_perfs(i)%name), perf_env%mp_perfs(i)%count
               ELSE
                  WRITE (iw, '(1X,A15,T17,I10,T40,F11.0)') &
                     ADJUSTL(perf_env%mp_perfs(i)%name), perf_env%mp_perfs(i)%count, &
                     vol
               END IF
            END IF

         END DO
         WRITE (iw, '( 1X, 79("-"), / )')
      END IF
#else
      MARK_USED(iw)
#endif
   END SUBROUTINE mp_perf_env_describe

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE rm_mp_perf_env()
      IF (stack_pointer < 1) THEN
         CPABORT("no perf_env in the stack : message_passing @ rm_mp_perf_env")
      END IF
      CALL mp_perf_env_release(mp_perf_stack(stack_pointer)%mp_perf_env)
      stack_pointer = stack_pointer - 1
   END SUBROUTINE rm_mp_perf_env

! **************************************************************************************************
!> \brief ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_mp_perf_env() RESULT(res)
      TYPE(mp_perf_env_type), POINTER                    :: res

      IF (stack_pointer < 1) THEN
         CPABORT("no perf_env in the stack : message_passing @ get_mp_perf_env")
      END IF
      res => mp_perf_stack(stack_pointer)%mp_perf_env
   END FUNCTION get_mp_perf_env

! **************************************************************************************************
!> \brief ...
!> \param scr ...
! **************************************************************************************************
   SUBROUTINE describe_mp_perf_env(scr)
      INTEGER, INTENT(in)                                :: scr

      TYPE(mp_perf_env_type), POINTER                    :: perf_env

      perf_env => get_mp_perf_env()
      CALL mp_perf_env_describe(perf_env, scr)
   END SUBROUTINE describe_mp_perf_env

! **************************************************************************************************
!> \brief adds the performance informations of one call
!> \param perf_id ...
!> \param count ...
!> \param msg_size ...
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE add_perf(perf_id, count, msg_size)
      INTEGER, INTENT(in)                      :: perf_id
      INTEGER, INTENT(in), OPTIONAL            :: count
      INTEGER, INTENT(in), OPTIONAL            :: msg_size

#if defined(__parallel)
      TYPE(mp_perf_type), POINTER              :: mp_perf

      IF (.NOT. ASSOCIATED(mp_perf_stack(stack_pointer)%mp_perf_env)) return

      mp_perf => mp_perf_stack(stack_pointer)%mp_perf_env%mp_perfs(perf_id)
      IF (PRESENT(count)) THEN
         mp_perf%count = mp_perf%count + count
      END IF
      IF (PRESENT(msg_size)) THEN
         mp_perf%msg_size = mp_perf%msg_size + REAL(msg_size, dp)
      END IF
#else
      MARK_USED(perf_id)
      MARK_USED(count)
      MARK_USED(msg_size)
#endif

   END SUBROUTINE add_perf

! **************************************************************************************************
!> \brief globally stops all tasks
!>       this is intended to be low level, most of CP2K should call cp_abort()
! **************************************************************************************************
   SUBROUTINE mp_abort()
      INTEGER                                            :: ierr

      ierr = 0

#if !defined(__NO_ABORT)
#if defined(__parallel)
      CALL mpi_abort(MPI_COMM_WORLD, 1, ierr)
#else
      CALL m_abort()
#endif
#endif
      ! this routine never returns and levels with non-zero exit code
      STOP 1
   END SUBROUTINE mp_abort

! **************************************************************************************************
!> \brief stops *after an mpi error* translating the error code
!> \param ierr an error code * returned by an mpi call *
!> \param prg_code ...
!> \note
!>       this function is private to message_passing.F
! **************************************************************************************************
   SUBROUTINE mp_stop(ierr, prg_code)
      INTEGER, INTENT(IN)                      :: ierr
      CHARACTER(LEN=*), INTENT(IN)                         :: prg_code

#if defined(__parallel)
      INTEGER                                  :: istat, len
      CHARACTER(LEN=MPI_MAX_ERROR_STRING)     :: error_string
      CHARACTER(LEN=MPI_MAX_ERROR_STRING + 512)  :: full_error
#else
      CHARACTER(LEN=512)                       :: full_error
#endif

#if defined(__parallel)
      CALL mpi_error_string(ierr, error_string, len, istat)
      WRITE (full_error, '(A,I0,A)') ' MPI error ', ierr, ' in '//TRIM(prg_code)//' : '//error_string(1:len)
#else
      WRITE (full_error, '(A,I0,A)') ' MPI error (!?) ', ierr, ' in '//TRIM(prg_code)
#endif

      CPABORT(full_error)

   END SUBROUTINE mp_stop

! **************************************************************************************************
!> \brief synchronizes with a barrier a given group of mpi tasks
!> \param group mpi communicator
! **************************************************************************************************
   SUBROUTINE mp_sync(comm)
      CLASS(mp_comm_type), INTENT(IN)                                :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_sync'

      INTEGER                                            :: handle, ierr

      ierr = 0
      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      CALL mpi_barrier(comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_barrier @ mp_sync")
      CALL add_perf(perf_id=5, count=1)
#else
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)

   END SUBROUTINE mp_sync

! **************************************************************************************************
!> \brief synchronizes with a barrier a given group of mpi tasks
!> \param comm mpi communicator
!> \param request ...
! **************************************************************************************************
   SUBROUTINE mp_isync(comm, request)
      CLASS(mp_comm_type), INTENT(IN)                                :: comm
      TYPE(mp_request_type), INTENT(OUT)                               :: request

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_isync'

      INTEGER                                            :: handle, ierr

      ierr = 0
      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      CALL mpi_ibarrier(comm%handle, request%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_ibarrier @ mp_isync")
      CALL add_perf(perf_id=26, count=1)
#else
      MARK_USED(comm)
      request = mp_request_null
#endif
      CALL mp_timestop(handle)

   END SUBROUTINE mp_isync

! **************************************************************************************************
!> \brief returns number of tasks and task id for a given mpi group
!>       simple and cartesian version.. recursive needed in case of failing mpi_comm_rank.
!> \param numtask ...
!> \param taskid ...
!> \param comm mpi communicator
!> \note
!>         ..mp_world_setup is gone, use mp_environ instead (i.e. give a groupid explicitly)
! **************************************************************************************************
   SUBROUTINE mp_environ_l(numtask, taskid, comm)

      INTEGER, INTENT(OUT)                               :: numtask, taskid
      CLASS(mp_comm_type), INTENT(IN)                                :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_environ_l'

      INTEGER                                            :: handle, ierr

      ierr = 0
      CALL mp_timeset(routineN, handle)

      numtask = 1
      taskid = 0
#if defined(__parallel)
      CALL mpi_comm_rank(comm%handle, taskid, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_comm_rank @ mp_environ_l")

      CALL mpi_comm_size(comm%handle, numtask, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_comm_size @ mp_environ_l")
#else
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)

   END SUBROUTINE mp_environ_l

! **************************************************************************************************
!> \brief ...
!> \param numtask ...
!> \param dims ...
!> \param task_coor ...
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE mp_environ_c(numtask, dims, task_coor, comm)

      INTEGER, INTENT(OUT)                     :: numtask, dims(2), task_coor(2)
      CLASS(mp_cart_type), INTENT(IN)                      :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_environ_c'

      INTEGER                                  :: handle, ierr
#if defined(__parallel)
      LOGICAL, DIMENSION(2)                    :: periods
#endif

      ierr = 0
      CALL mp_timeset(routineN, handle)
      numtask = 1
      task_coor = 0
      dims = 1
#if defined(__parallel)
      CALL mpi_comm_size(comm%handle, numtask, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_comm_size @ mp_environ_c")

      CALL mpi_cart_get(comm%handle, 2, dims, periods, task_coor, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_cart_get @ mp_environ_c")
#else
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)

   END SUBROUTINE mp_environ_c

! **************************************************************************************************
!> \brief ...
!> \param comm ...
!> \param ndims ...
!> \param dims ...
!> \param task_coor ...
!> \param periods ...
! **************************************************************************************************
   SUBROUTINE mp_environ_c2(comm, ndims, dims, task_coor, periods)

      CLASS(mp_cart_type), INTENT(IN) :: comm
      INTEGER, INTENT(IN)                                :: ndims
      INTEGER, INTENT(OUT)                               :: dims(ndims), task_coor(ndims)
      LOGICAL, INTENT(out)                               :: periods(ndims)

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_environ_c2'

      INTEGER                                            :: handle, ierr

      ierr = 0
      CALL mp_timeset(routineN, handle)

      task_coor = 0
      dims = 1
      periods = .FALSE.
#if defined(__parallel)
      CALL mpi_cart_get(comm%handle, ndims, dims, periods, task_coor, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_cart_get @ mp_environ_c")
#else
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)

   END SUBROUTINE mp_environ_c2

! **************************************************************************************************
!> \brief ...
!> \param comm_old ...
!> \param ndims ...
!> \param dims ...
!> \param pos ...
!> \param comm_cart ...
! **************************************************************************************************
   SUBROUTINE mp_cart_create(comm_old, ndims, dims, pos, comm_cart)

      CLASS(mp_comm_type), INTENT(IN) :: comm_old
      INTEGER, INTENT(IN)                      :: ndims
      INTEGER, INTENT(INOUT), CONTIGUOUS                   :: dims(:)
      INTEGER, INTENT(OUT), CONTIGUOUS                     :: pos(:)
      CLASS(mp_cart_type), INTENT(OUT) :: comm_cart

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_cart_create'

      INTEGER                                  :: handle, ierr, nodes
#if defined(__parallel)
      LOGICAL, DIMENSION(1:ndims)              :: period
      LOGICAL                                  :: reorder
#endif

      ierr = 0
      CALL mp_timeset(routineN, handle)

      nodes = 0
      pos(1:ndims) = -1
      comm_cart%handle = comm_old%handle
#if defined(__parallel)

      CALL mpi_comm_size(comm_old%handle, nodes, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_comm_size @ mp_cart_create")

      IF (ANY(dims == 0)) CALL mpi_dims_create(nodes, ndims, dims, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_dims_create @ mp_cart_create")

      ! FIX ME.  Quick hack to avoid problems with realspace grids for compilers
      ! like IBM that actually reorder the processors when creating the new
      ! communicator
      reorder = .FALSE.
      period = .TRUE.
      CALL mpi_cart_create(comm_old%handle, ndims, dims, period, reorder, comm_cart%handle, &
                           ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_cart_create @ mp_cart_create")

      IF (comm_cart%handle /= MPI_COMM_NULL) THEN
         CALL mpi_cart_get(comm_cart%handle, ndims, dims, period, pos, ierr)
         IF (ierr /= 0) CALL mp_stop(ierr, "mpi_cart_get @ mp_cart_create")
      END IF
      CALL add_perf(perf_id=1, count=1)
#else
      pos(1:ndims) = 0
      dims = 1
      comm_cart%handle = mp_comm_default_handle
#endif
      debug_comm_count = debug_comm_count + 1
      CALL mp_timestop(handle)

   END SUBROUTINE mp_cart_create

! **************************************************************************************************
!> \brief ...
!> \param comm ...
!> \param rank ...
!> \param coords ...
! **************************************************************************************************
   SUBROUTINE mp_cart_coords(comm, rank, coords)

      CLASS(mp_cart_type), INTENT(IN) :: comm
      INTEGER, INTENT(IN)                                :: rank
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: coords

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_cart_coords'

      INTEGER                                            :: handle, ierr, m

      ierr = 0
      CALL mp_timeset(routineN, handle)

      m = SIZE(coords)
#if defined(__parallel)
      CALL mpi_cart_coords(comm%handle, rank, m, coords, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_cart_coords @ mp_cart_coords")
#else
      coords = 0
      MARK_USED(rank)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)

   END SUBROUTINE mp_cart_coords

! **************************************************************************************************
!> \brief ...
!> \param comm1 ...
!> \param comm2 ...
!> \param res ...
! **************************************************************************************************
   FUNCTION mp_comm_compare(comm1, comm2) RESULT(res)

      CLASS(mp_comm_type), INTENT(IN)                                :: comm1, comm2
      INTEGER                               :: res

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_comm_compare'

      INTEGER                                            :: handle
#if defined(__parallel)
      INTEGER :: ierr, iout
#endif

      CALL mp_timeset(routineN, handle)

      res = 0
#if defined(__parallel)
      CALL mpi_comm_compare(comm1%handle, comm2%handle, iout, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_comm_compare @ mp_comm_compare")
      SELECT CASE (iout)
      CASE (MPI_IDENT)
         res = mp_comm_ident
      CASE (MPI_CONGRUENT)
         res = mp_comm_congruent
      CASE (MPI_SIMILAR)
         res = mp_comm_similar
      CASE (MPI_UNEQUAL)
         res = mp_comm_unequal
      CASE default
         CPABORT("Unknown comparison state of the communicators!")
      END SELECT
#else
      MARK_USED(comm1)
      MARK_USED(comm2)
#endif
      CALL mp_timestop(handle)

   END FUNCTION mp_comm_compare

! **************************************************************************************************
!> \brief ...
!> \param comm ...
!> \param rdim ...
!> \param sub_comm ...
! **************************************************************************************************
   SUBROUTINE mp_cart_sub(comm, rdim, sub_comm)

      CLASS(mp_cart_type), INTENT(IN)                                :: comm
      LOGICAL, DIMENSION(:), CONTIGUOUS, INTENT(IN)                  :: rdim
      CLASS(mp_cart_type), INTENT(OUT)                               :: sub_comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_cart_sub'

      INTEGER                                            :: handle
#if defined(__parallel)
      INTEGER :: ierr
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      CALL mpi_cart_sub(comm%handle, rdim, sub_comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_cart_sub @ mp_cart_sub")
#else
      MARK_USED(comm)
      MARK_USED(rdim)
      sub_comm%handle = mp_comm_default_handle
#endif
      debug_comm_count = debug_comm_count + 1
      CALL mp_timestop(handle)

   END SUBROUTINE mp_cart_sub

! **************************************************************************************************
!> \brief ...
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE mp_comm_free(comm)

      CLASS(mp_comm_type), INTENT(INOUT)                             :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_comm_free'

      INTEGER                                            :: handle
#if defined(__parallel)
      INTEGER :: ierr
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      CALL mpi_comm_free(comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_comm_free @ mp_comm_free")
#else
      comm%handle = mp_comm_null_handle
#endif
      debug_comm_count = debug_comm_count - 1
      CALL mp_timestop(handle)

   END SUBROUTINE mp_comm_free

! **************************************************************************************************
!> \brief ...
!> \param comm1 ...
!> \param comm2 ...
! **************************************************************************************************
   SUBROUTINE mp_comm_dup(comm1, comm2)

      CLASS(mp_comm_type), INTENT(IN)                                :: comm1
      CLASS(mp_comm_type), INTENT(OUT)                               :: comm2

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_comm_dup'

      INTEGER                                            :: handle
#if defined(__parallel)
      INTEGER :: ierr
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      CALL mpi_comm_dup(comm1%handle, comm2%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_comm_dup @ mp_comm_dup")
#else
      MARK_USED(comm1)
      comm2%handle = mp_comm_default_handle
#endif
      debug_comm_count = debug_comm_count + 1
      CALL mp_timestop(handle)

   END SUBROUTINE mp_comm_dup

   SUBROUTINE mp_comm_copy(comm_new, comm_old)
      CLASS(mp_comm_type), INTENT(IN) :: comm_old
      CLASS(mp_comm_type), INTENT(OUT) :: comm_new

      comm_new%handle = comm_old%handle
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param comm1 ...
!> \param comm2 ...
!> \param rank ...
! **************************************************************************************************
   SUBROUTINE mp_rank_compare(comm1, comm2, rank)

      CLASS(mp_comm_type), INTENT(IN)                      :: comm1, comm2
      INTEGER, DIMENSION(:), CONTIGUOUS, INTENT(OUT)       :: rank

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_rank_compare'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: i, ierr, n, n1, n2
      INTEGER, ALLOCATABLE, DIMENSION(:)       :: rin
      MPI_GROUP_TYPE :: g1, g2
#endif

      CALL mp_timeset(routineN, handle)

      rank = 0
#if defined(__parallel)
      CALL mpi_comm_size(comm1%handle, n1, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_comm_size @ mp_rank_compare")
      CALL mpi_comm_size(comm2%handle, n2, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_comm_size @ mp_rank_compare")
      n = MAX(n1, n2)
      CALL mpi_comm_group(comm1%handle, g1, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_comm_group @ mp_rank_compare")
      CALL mpi_comm_group(comm2%handle, g2, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_comm_group @ mp_rank_compare")
      ALLOCATE (rin(0:n - 1), STAT=ierr)
      IF (ierr /= 0) &
         CPABORT("allocate @ mp_rank_compare")
      DO i = 0, n - 1
         rin(i) = i
      END DO
      CALL mpi_group_translate_ranks(g1, n, rin, g2, rank, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, &
                                  "mpi_group_translate_rank @ mp_rank_compare")
      CALL mpi_group_free(g1, ierr)
      IF (ierr /= 0) &
         CPABORT("group_free @ mp_rank_compare")
      CALL mpi_group_free(g2, ierr)
      IF (ierr /= 0) &
         CPABORT("group_free @ mp_rank_compare")
      DEALLOCATE (rin)
#else
      MARK_USED(comm1)
      MARK_USED(comm2)
#endif
      CALL mp_timestop(handle)

   END SUBROUTINE mp_rank_compare

! **************************************************************************************************
!> \brief ...
!> \param nodes ...
!> \param dims ...
! **************************************************************************************************
   SUBROUTINE mp_dims_create(nodes, dims)

      INTEGER, INTENT(IN)                                :: nodes
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: dims

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_dims_create'

      INTEGER                                            :: handle, ndim
#if defined(__parallel)
      INTEGER :: ierr
#endif

      CALL mp_timeset(routineN, handle)

      ndim = SIZE(dims)
#if defined(__parallel)
      IF (ANY(dims == 0)) CALL mpi_dims_create(nodes, ndim, dims, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_dims_create @ mp_dims_create")
#else
      dims = 1
      MARK_USED(nodes)
#endif
      CALL mp_timestop(handle)

   END SUBROUTINE mp_dims_create

! **************************************************************************************************
!> \brief ...
!> \param comm ...
!> \param pos ...
!> \param rank ...
! **************************************************************************************************
   SUBROUTINE mp_cart_rank(comm, pos, rank)
      CLASS(mp_cart_type), INTENT(IN)                                :: comm
      INTEGER, DIMENSION(:), INTENT(IN)                  :: pos
      INTEGER, INTENT(OUT)                               :: rank

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_cart_rank'

      INTEGER                                            :: handle
#if defined(__parallel)
      INTEGER :: ierr
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      CALL mpi_cart_rank(comm%handle, pos, rank, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_cart_rank @ mp_cart_rank")
#else
      rank = 0
      MARK_USED(comm)
      MARK_USED(pos)
#endif
      CALL mp_timestop(handle)

   END SUBROUTINE mp_cart_rank

! **************************************************************************************************
!> \brief waits for completion of the given request
!> \param request ...
!> \par History
!>      08.2003 created [f&j]
!> \author joost & fawzi
!> \note
!>      see isendrecv
! **************************************************************************************************
   SUBROUTINE mp_wait(request)
      CLASS(mp_request_type), INTENT(inout)                             :: request

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_wait'

      INTEGER                                            :: handle
#if defined(__parallel)
      INTEGER :: ierr
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)

      CALL mpi_wait(request%handle, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_wait @ mp_wait")

      CALL add_perf(perf_id=9, count=1)
#else
      request%handle = mp_request_null_handle
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_wait

! **************************************************************************************************
!> \brief waits for completion of the given requests
!> \param requests ...
!> \par History
!>      08.2003 created [f&j]
!> \author joost & fawzi
!> \note
!>      see isendrecv
! **************************************************************************************************
   SUBROUTINE mp_waitall_1(requests)
      TYPE(mp_request_type), DIMENSION(:), INTENT(inout)     :: requests

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_waitall_1'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: count, ierr
#if !defined(__MPI_F08)
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: status
#else
      TYPE(MPI_Status), ALLOCATABLE, DIMENSION(:)    :: status
#endif
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      count = SIZE(requests)
#if !defined(__MPI_F08)
      ALLOCATE (status(MPI_STATUS_SIZE, count))
#else
      ALLOCATE (status(count))
#endif
      CALL mpi_waitall_internal(count, requests, status, ierr) ! MPI_STATUSES_IGNORE openmpi workaround
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_waitall @ mp_waitall_1")
      DEALLOCATE (status)
      CALL add_perf(perf_id=9, count=1)
#else
      requests = mp_request_null
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_waitall_1

! **************************************************************************************************
!> \brief waits for completion of the given requests
!> \param requests ...
!> \par History
!>      08.2003 created [f&j]
!> \author joost & fawzi
! **************************************************************************************************
   SUBROUTINE mp_waitall_2(requests)
      TYPE(mp_request_type), DIMENSION(:, :), INTENT(inout)  :: requests

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_waitall_2'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: count, ierr
#if !defined(__MPI_F08)
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: status
#else
      TYPE(MPI_Status), ALLOCATABLE, DIMENSION(:)    :: status
#endif
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      count = SIZE(requests)
#if !defined(__MPI_F08)
      ALLOCATE (status(MPI_STATUS_SIZE, count))
#else
      ALLOCATE (status(count))
#endif

      CALL mpi_waitall_internal(count, requests, status, ierr) ! MPI_STATUSES_IGNORE openmpi workaround
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_waitall @ mp_waitall_2")
      DEALLOCATE (status)

      CALL add_perf(perf_id=9, count=1)
#else
      requests = mp_request_null
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_waitall_2

! **************************************************************************************************
!> \brief wrapper needed to deal with interfaces as present in openmpi 1.8.1
!>        the issue is with the rank or requests
!> \param count ...
!> \param array_of_requests ...
!> \param array_of_statuses ...
!> \param ierr ...
!> \author Joost VandeVondele
! **************************************************************************************************
#if defined(__parallel)
   SUBROUTINE mpi_waitall_internal(count, array_of_requests, array_of_statuses, ierr)
      INTEGER, INTENT(in)                                :: count
      TYPE(mp_request_type), DIMENSION(count), INTENT(inout)           :: array_of_requests
#if !defined(__MPI_F08)
      INTEGER, DIMENSION(MPI_STATUS_SIZE, count), &
         INTENT(out)                                     :: array_of_statuses
#else
      TYPE(MPI_Status), DIMENSION(count), &
         INTENT(out)                                     :: array_of_statuses
#endif
      INTEGER, INTENT(out)                               :: ierr

      INTEGER :: i
      MPI_REQUEST_TYPE, ALLOCATABLE, DIMENSION(:) :: request_handles

      ALLOCATE (request_handles(count))
      DO i = 1, count
         request_handles(i) = array_of_requests(i)%handle
      END DO

      CALL mpi_waitall(count, request_handles, array_of_statuses, ierr)

      DO i = 1, count
         array_of_requests(i)%handle = request_handles(i)
      END DO

   END SUBROUTINE mpi_waitall_internal
#endif

! **************************************************************************************************
!> \brief waits for completion of any of the given requests
!> \param requests ...
!> \param completed ...
!> \par History
!>      09.2008 created
!> \author Iain Bethune (c) The Numerical Algorithms Group (NAG) Ltd, 2008 on behalf of the HECToR project
! **************************************************************************************************
   SUBROUTINE mp_waitany(requests, completed)
      TYPE(mp_request_type), DIMENSION(:), INTENT(inout)     :: requests
      INTEGER, INTENT(out)                     :: completed

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_waitany'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: count, i, ierr
      MPI_REQUEST_TYPE, ALLOCATABLE, DIMENSION(:) :: request_handles
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      count = SIZE(requests)
! Convert CP2K's request_handles to the plane handle for the library
! (Maybe, the compiler optimizes it away)
      ALLOCATE (request_handles(count))
      DO i = 1, count
         request_handles(i) = requests(i)%handle
      END DO
      CALL mpi_waitany(count, request_handles, completed, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_waitany @ mp_waitany")
! Convert the plane handles to CP2K handles
      DO i = 1, count
         requests(i)%handle = request_handles(i)
      END DO
      CALL add_perf(perf_id=9, count=1)
#else
      requests = mp_request_null
      completed = 1
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_waitany

! **************************************************************************************************
!> \brief Tests for completion of the given requests.
!> \brief We use mpi_test so that we can use a single status.
!> \param requests the list of requests to test
!> \return logical which determines if requests are complete
!> \par History
!>      3.2016 adapted to any shape [Nico Holmberg]
!> \author Alfio Lazzaro
! **************************************************************************************************
   FUNCTION mp_testall_tv(requests) RESULT(flag)
      TYPE(mp_request_type), DIMENSION(:), INTENT(INOUT) :: requests
      LOGICAL                               :: flag

#if defined(__parallel)
      INTEGER                               :: i, ierr
      LOGICAL, DIMENSION(:), POINTER        :: flags
#endif

      flag = .TRUE.

#if defined(__parallel)
      ALLOCATE (flags(SIZE(requests)))
      DO i = 1, SIZE(requests)
         CALL mpi_test(requests(i)%handle, flags(i), MPI_STATUS_IGNORE, ierr)
         IF (ierr /= 0) CALL mp_stop(ierr, "mpi_testall @ mp_testall_tv")
         flag = flag .AND. flags(i)
      END DO
      DEALLOCATE (flags)
#else
      requests = mp_request_null
#endif
   END FUNCTION mp_testall_tv

! **************************************************************************************************
!> \brief Tests for completion of the given request.
!> \param request the request
!> \param flag logical which determines if the request is completed
!> \par History
!>      3.2016 created
!> \author Nico Holmberg
! **************************************************************************************************
   FUNCTION mp_test_1(request) RESULT(flag)
      CLASS(mp_request_type), INTENT(inout)                             :: request
      LOGICAL                                            :: flag

#if defined(__parallel)
      INTEGER                                            :: ierr

      CALL mpi_test(request%handle, flag, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_test @ mp_test_1")
#else
      MARK_USED(request)
      flag = .TRUE.
#endif
   END FUNCTION mp_test_1

! **************************************************************************************************
!> \brief tests for completion of the given requests
!> \param requests ...
!> \param completed ...
!> \param flag ...
!> \par History
!>      08.2011 created
!> \author Iain Bethune
! **************************************************************************************************
   SUBROUTINE mp_testany_1(requests, completed, flag)
      TYPE(mp_request_type), DIMENSION(:), INTENT(inout)  :: requests
      INTEGER, INTENT(out), OPTIONAL           :: completed
      LOGICAL, INTENT(out), OPTIONAL           :: flag

#if defined(__parallel)
      INTEGER                                  :: completed_l, count, ierr
      LOGICAL                                  :: flag_l

      count = SIZE(requests)

      CALL mpi_testany_internal(count, requests, completed_l, flag_l, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_testany_1 @ mp_testany")

      IF (PRESENT(completed)) completed = completed_l
      IF (PRESENT(flag)) flag = flag_l
#else
      MARK_USED(requests)
      IF (PRESENT(completed)) completed = 1
      IF (PRESENT(flag)) flag = .TRUE.
#endif
   END SUBROUTINE mp_testany_1

! **************************************************************************************************
!> \brief tests for completion of the given requests
!> \param requests ...
!> \param completed ...
!> \param flag ...
!> \par History
!>      08.2011 created
!> \author Iain Bethune
! **************************************************************************************************
   SUBROUTINE mp_testany_2(requests, completed, flag)
      TYPE(mp_request_type), DIMENSION(:, :), INTENT(inout)   :: requests
      INTEGER, INTENT(out), OPTIONAL           :: completed
      LOGICAL, INTENT(out), OPTIONAL           :: flag

#if defined(__parallel)
      INTEGER                                  :: completed_l, count, ierr
      LOGICAL                                  :: flag_l

      count = SIZE(requests)

      CALL mpi_testany_internal(count, requests, completed_l, flag_l, MPI_STATUS_IGNORE, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_testany_2 @ mp_testany")

      IF (PRESENT(completed)) completed = completed_l
      IF (PRESENT(flag)) flag = flag_l
#else
      MARK_USED(requests)
      IF (PRESENT(completed)) completed = 1
      IF (PRESENT(flag)) flag = .TRUE.
#endif
   END SUBROUTINE mp_testany_2

! **************************************************************************************************
!> \brief wrapper needed to deal with interfaces as present in openmpi 1.8.1
!>        the issue is with the rank or requests
!> \param count ...
!> \param array_of_requests ...
!> \param index ...
!> \param flag ...
!> \param status ...
!> \param ierr ...
!> \author Joost VandeVondele
! **************************************************************************************************
#if defined(__parallel)
   SUBROUTINE mpi_testany_internal(count, array_of_requests, index, flag, status, ierr)
      INTEGER, INTENT(in)                                :: count
      TYPE(mp_request_type), DIMENSION(count), INTENT(inout)           :: array_of_requests
      INTEGER, INTENT(out)                               :: index
      LOGICAL, INTENT(out)                               :: flag
      MPI_STATUS_TYPE, INTENT(out)   :: status
      INTEGER, INTENT(out)                               :: ierr

      INTEGER :: i
      MPI_REQUEST_TYPE, ALLOCATABLE, DIMENSION(:) :: request_handles

      ALLOCATE (request_handles(count))
      DO i = 1, count
         request_handles(i) = array_of_requests(i)%handle
      END DO

      CALL mpi_testany(count, request_handles, index, flag, status, ierr)

      DO i = 1, count
         array_of_requests(i)%handle = request_handles(i)
      END DO

   END SUBROUTINE mpi_testany_internal
#endif

! **************************************************************************************************
!> \brief the direct way to split a communicator each color is a sub_comm,
!>        the rank order is according to the order in the orig comm
!> \param comm ...
!> \param sub_comm ...
!> \param color ...
!> \param key ...
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE mp_comm_split_direct(comm, sub_comm, color, key)
      CLASS(mp_comm_type), INTENT(in)                                :: comm
      CLASS(mp_comm_type), INTENT(OUT)                               :: sub_comm
      INTEGER, INTENT(in)                                :: color
      INTEGER, INTENT(in), OPTIONAL                      :: key

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_comm_split_direct'

      INTEGER                                            :: handle
#if defined(__parallel)
      INTEGER :: ierr, my_key
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      my_key = 0
      IF (PRESENT(key)) my_key = key
      CALL mpi_comm_split(comm%handle, color, my_key, sub_comm%handle, ierr)
      IF (ierr /= mpi_success) CALL mp_stop(ierr, routineN)
      CALL add_perf(perf_id=10, count=1)
#else
      sub_comm%handle = mp_comm_default_handle
      MARK_USED(comm)
      MARK_USED(color)
      MARK_USED(key)
#endif
      debug_comm_count = debug_comm_count + 1
      CALL mp_timestop(handle)

   END SUBROUTINE mp_comm_split_direct
! **************************************************************************************************
!> \brief splits the given communicator in group in subgroups trying to organize
!>      them in a way that the communication within each subgroup is
!>      efficient (but not necessarily the communication between subgroups)
!> \param comm the mpi communicator that you want to split
!> \param sub_comm the communicator for the subgroup (created, needs to be freed later)
!> \param ngroups actual number of groups
!> \param group_distribution input  : allocated with array with the nprocs entries (0 .. nprocs-1)
!> \param subgroup_min_size the minimum size of the subgroup
!> \param n_subgroups the number of subgroups wanted
!> \param group_partition n_subgroups sized array containing the number of cpus wanted per group.
!>                         should match the total number of cpus (only used if present and associated) (0..ngroups-1)
!> \param stride create groups using a stride (default=1) through the ranks of the comm to be split.
!> \par History
!>      10.2003 created [fawzi]
!>      02.2004 modified [Joost VandeVondele]
!> \author Fawzi Mohamed
!> \note
!>      at least one of subgroup_min_size and n_subgroups is needed,
!>      the other default to the value needed to use most processors.
!>      if less cpus are present than needed for subgroup min size, n_subgroups,
!>      just one comm is created that contains all cpus
! **************************************************************************************************
   SUBROUTINE mp_comm_split(comm, sub_comm, ngroups, group_distribution, &
                            subgroup_min_size, n_subgroups, group_partition, stride)
      CLASS(mp_comm_type), INTENT(in)                      :: comm
      CLASS(mp_comm_type), INTENT(out) :: sub_comm
      INTEGER, INTENT(out)                     :: ngroups
      INTEGER, DIMENSION(0:), INTENT(INOUT)                    :: group_distribution
      INTEGER, INTENT(in), OPTIONAL            :: subgroup_min_size, n_subgroups
      INTEGER, DIMENSION(0:), INTENT(IN), OPTIONAL          :: group_partition
      INTEGER, OPTIONAL, INTENT(IN)                        :: stride

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp_comm_split', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: handle, mepos, nnodes
#if defined(__parallel)
      INTEGER                                  :: color, i, ierr, j, k, &
                                                  my_subgroup_min_size, &
                                                  istride, local_stride, irank
      INTEGER, DIMENSION(:), ALLOCATABLE       :: rank_permutation
#endif

      CALL mp_timeset(routineN, handle)

      ! actual number of groups

      IF (.NOT. PRESENT(subgroup_min_size) .AND. .NOT. PRESENT(n_subgroups)) THEN
         CPABORT(routineP//" missing arguments")
      END IF
      IF (PRESENT(subgroup_min_size) .AND. PRESENT(n_subgroups)) THEN
         CPABORT(routineP//" too many arguments")
      END IF

      CALL comm%environ(nnodes, mepos)

      IF (UBOUND(group_distribution, 1) .NE. nnodes - 1) THEN
         CPABORT(routineP//" group_distribution wrong bounds")
      END IF

#if defined(__parallel)
      IF (PRESENT(subgroup_min_size)) THEN
         IF (subgroup_min_size < 0 .OR. subgroup_min_size > nnodes) THEN
            CPABORT(routineP//" subgroup_min_size too small or too large")
         END IF
         ngroups = nnodes/subgroup_min_size
         my_subgroup_min_size = subgroup_min_size
      ELSE ! n_subgroups
         IF (n_subgroups <= 0) THEN
            CPABORT(routineP//" n_subgroups too small")
         END IF
         IF (nnodes/n_subgroups > 0) THEN ! we have a least one cpu per group
            ngroups = n_subgroups
         ELSE ! well, only one group then
            ngroups = 1
         END IF
         my_subgroup_min_size = nnodes/ngroups
      END IF

      ! rank_permutation: is a permutation of ranks, so that groups are not necessarily continuous in rank of the master group
      ! while the order is not critical (we only color ranks), it can e.g. be used to make groups that have just 1 rank per node
      ! (by setting stride equal to the number of mpi ranks per node), or by sharing  a node between two groups (stride 2).
      ALLOCATE (rank_permutation(0:nnodes - 1))
      local_stride = 1
      IF (PRESENT(stride)) local_stride = stride
      k = 0
      DO istride = 1, local_stride
         DO irank = istride - 1, nnodes - 1, local_stride
            rank_permutation(k) = irank
            k = k + 1
         END DO
      END DO

      DO i = 0, nnodes - 1
         group_distribution(rank_permutation(i)) = MIN(i/my_subgroup_min_size, ngroups - 1)
      END DO
      ! even the user gave a partition, see if we can use it to overwrite this choice
      IF (PRESENT(group_partition)) THEN
         IF (ALL(group_partition > 0) .AND. (SUM(group_partition) .EQ. nnodes) .AND. (ngroups == SIZE(group_partition))) THEN
            k = 0
            DO i = 0, SIZE(group_partition) - 1
               DO j = 1, group_partition(i)
                  group_distribution(rank_permutation(k)) = i
                  k = k + 1
               END DO
            END DO
         ELSE
            ! just ignore silently as we have reasonable defaults. Probably a warning would not be to bad
         END IF
      END IF
      color = group_distribution(mepos)
      CALL mpi_comm_split(comm%handle, color, 0, sub_comm%handle, ierr)
      IF (ierr /= mpi_success) CALL mp_stop(ierr, "in "//routineP//" split")

      CALL add_perf(perf_id=10, count=1)
#else
      sub_comm%handle = mp_comm_default_handle
      group_distribution(0) = 0
      ngroups = 1
      MARK_USED(comm)
      MARK_USED(stride)
      MARK_USED(group_partition)
#endif
      debug_comm_count = debug_comm_count + 1
      CALL mp_timestop(handle)

   END SUBROUTINE mp_comm_split

! **************************************************************************************************
!> \brief Get the local rank on the node according to the global communicator
!> \return Node Rank id
!> \author Alfio Lazzaro
! **************************************************************************************************
   FUNCTION mp_get_node_global_rank() &
      RESULT(node_rank)

      INTEGER                                            :: node_rank

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_get_node_global_rank'
      INTEGER                                            :: handle
#if defined(__parallel)
      INTEGER :: ierr, rank
      TYPE(mp_comm_type) :: comm
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      CALL mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
      IF (ierr /= mpi_success) CALL mp_stop(ierr, routineN)
      CALL mpi_comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, comm%handle, ierr)
      IF (ierr /= mpi_success) CALL mp_stop(ierr, routineN)
      CALL mpi_comm_rank(comm%handle, node_rank, ierr)
      IF (ierr /= mpi_success) CALL mp_stop(ierr, routineN)
      CALL mpi_comm_free(comm%handle, ierr)
      IF (ierr /= mpi_success) CALL mp_stop(ierr, routineN)
#else
      node_rank = 0
#endif
      CALL mp_timestop(handle)

   END FUNCTION mp_get_node_global_rank

! **************************************************************************************************
!> \brief probes for an incoming message with any tag
!> \param[inout] source the source of the possible incoming message,
!>        if MP_ANY_SOURCE it is a blocking one and return value is the source
!>        of the next incoming message
!>        if source is a different value it is a non-blocking probe returning
!>        MP_ANY_SOURCE if there is no incoming message
!> \param[in] comm the communicator
!> \param[out] tag the tag of the incoming message
!> \author Mandes
! **************************************************************************************************
   SUBROUTINE mp_probe(source, comm, tag)
      INTEGER, INTENT(INOUT)                      :: source
      CLASS(mp_comm_type), INTENT(IN)           :: comm
      INTEGER, INTENT(OUT)                     :: tag

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_probe'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER :: ierr
      MPI_STATUS_TYPE     :: status_single
      LOGICAL                                  :: flag
#endif

!   ---------------------------------------------------------------------------

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      IF (source .EQ. mp_any_source) THEN
         CALL mpi_probe(mp_any_source, mp_any_tag, comm%handle, status_single, ierr)
         IF (ierr /= 0) CALL mp_stop(ierr, "mpi_probe @ mp_probe")
         source = status_single MPI_STATUS_EXTRACT(MPI_SOURCE)
         tag = status_single MPI_STATUS_EXTRACT(MPI_TAG)
      ELSE
         flag = .FALSE.
         CALL mpi_iprobe(source, mp_any_tag, comm%handle, flag, status_single, ierr)
         IF (ierr /= 0) CALL mp_stop(ierr, "mpi_iprobe @ mp_probe")
         IF (flag .EQV. .FALSE.) THEN
            source = mp_any_source
            tag = -1 !status_single(MPI_TAG) ! in case of flag==false status is undefined
         ELSE
            tag = status_single MPI_STATUS_EXTRACT(MPI_TAG)
         END IF
      END IF
#else
      tag = -1
      MARK_USED(comm)
      MARK_USED(source)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_probe

! **************************************************************************************************
! Here come the data routines with none of the standard data types.
! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param msg ...
!> \param source ...
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE mp_bcast_b(msg, source, comm)
      LOGICAL, INTENT(INOUT)                                            :: msg
      INTEGER, INTENT(IN)                                            :: source
      CLASS(mp_comm_type), INTENT(IN) :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_b'

      INTEGER                                            :: handle
#if defined(__parallel)
      INTEGER :: ierr, msglen
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      msglen = 1
      CALL mpi_bcast(msg, msglen, MPI_LOGICAL, source, comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_bcast @ "//routineN)
      CALL add_perf(perf_id=2, count=1, msg_size=msglen*loglen)
#else
      MARK_USED(msg)
      MARK_USED(source)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_bcast_b

! **************************************************************************************************
!> \brief ...
!> \param msg ...
!> \param source ...
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE mp_bcast_bv(msg, source, comm)
      LOGICAL, CONTIGUOUS, INTENT(INOUT)                 :: msg(:)
      INTEGER, INTENT(IN)                                            :: source
      CLASS(mp_comm_type), INTENT(IN) :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_bv'

      INTEGER                                            :: handle
#if defined(__parallel)
      INTEGER :: ierr, msglen
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      msglen = SIZE(msg)
      CALL mpi_bcast(msg, msglen, MPI_LOGICAL, source, comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_bcast @ "//routineN)
      CALL add_perf(perf_id=2, count=1, msg_size=msglen*loglen)
#else
      MARK_USED(msg)
      MARK_USED(source)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_bcast_bv

! **************************************************************************************************
!> \brief Non-blocking send of logical vector data
!> \param msgin the input message
!> \param dest the destination processor
!> \param comm  the communicator object
!> \param request communication request index
!> \param tag message tag
!> \par History
!>      3.2016 added _bv subroutine [Nico Holmberg]
!> \author fawzi
!> \note see mp_irecv_iv
!> \note
!>      arrays can be pointers or assumed shape, but they must be contiguous!
! **************************************************************************************************
   SUBROUTINE mp_isend_bv(msgin, dest, comm, request, tag)
      LOGICAL, DIMENSION(:), INTENT(IN)        :: msgin
      INTEGER, INTENT(IN)                      :: dest
      CLASS(mp_comm_type), INTENT(IN) :: comm
      TYPE(mp_request_type), INTENT(out)                     :: request
      INTEGER, INTENT(in), OPTIONAL            :: tag

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_isend_bv'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: ierr, msglen, my_tag
      LOGICAL                                  :: foo(1)
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
#if !defined(__GNUC__) || __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(msgin))
#endif

      my_tag = 0
      IF (PRESENT(tag)) my_tag = tag

      msglen = SIZE(msgin, 1)
      IF (msglen > 0) THEN
         CALL mpi_isend(msgin(1), msglen, MPI_LOGICAL, dest, my_tag, &
                        comm%handle, request%handle, ierr)
      ELSE
         CALL mpi_isend(foo, msglen, MPI_LOGICAL, dest, my_tag, &
                        comm%handle, request%handle, ierr)
      END IF
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_isend @ "//routineN)

      CALL add_perf(perf_id=11, count=1, msg_size=msglen*loglen)
#else
      CPABORT("mp_isend called in non parallel case")
      MARK_USED(msgin)
      MARK_USED(dest)
      MARK_USED(comm)
      MARK_USED(tag)
      request = mp_request_null
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_isend_bv

! **************************************************************************************************
!> \brief Non-blocking receive of logical vector data
!> \param msgout the received message
!> \param source the source processor
!> \param comm  the communicator object
!> \param request communication request index
!> \param tag message tag
!> \par History
!>      3.2016 added _bv subroutine [Nico Holmberg]
!> \author fawzi
!> \note see mp_irecv_iv
!> \note
!>      arrays can be pointers or assumed shape, but they must be contiguous!
! **************************************************************************************************
   SUBROUTINE mp_irecv_bv(msgout, source, comm, request, tag)
      LOGICAL, DIMENSION(:), INTENT(INOUT)                    :: msgout
      INTEGER, INTENT(IN)                      :: source
      CLASS(mp_comm_type), INTENT(IN) :: comm
      TYPE(mp_request_type), INTENT(out)                     :: request
      INTEGER, INTENT(in), OPTIONAL            :: tag

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_irecv_bv'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: ierr, msglen, my_tag
      LOGICAL                                  :: foo(1)
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
#if !defined(__GNUC__) || __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(msgout))
#endif

      my_tag = 0
      IF (PRESENT(tag)) my_tag = tag

      msglen = SIZE(msgout, 1)
      IF (msglen > 0) THEN
         CALL mpi_irecv(msgout(1), msglen, MPI_LOGICAL, source, my_tag, &
                        comm%handle, request%handle, ierr)
      ELSE
         CALL mpi_irecv(foo, msglen, MPI_LOGICAL, source, my_tag, &
                        comm%handle, request%handle, ierr)
      END IF
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_ircv @ "//routineN)

      CALL add_perf(perf_id=12, count=1, msg_size=msglen*loglen)
#else
      CPABORT("mp_irecv called in non parallel case")
      MARK_USED(msgout)
      MARK_USED(source)
      MARK_USED(comm)
      MARK_USED(tag)
      request = mp_request_null
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_irecv_bv

! **************************************************************************************************
!> \brief Non-blocking send of rank-3 logical data
!> \param msgin the input message
!> \param dest the destination processor
!> \param comm  the communicator object
!> \param request communication request index
!> \param tag message tag
!> \par History
!>      2.2016 added _bm3 subroutine [Nico Holmberg]
!> \author fawzi
!> \note see mp_irecv_iv
!> \note
!>      arrays can be pointers or assumed shape, but they must be contiguous!
! **************************************************************************************************
   SUBROUTINE mp_isend_bm3(msgin, dest, comm, request, tag)
      LOGICAL, DIMENSION(:, :, :), INTENT(INOUT)              :: msgin
      INTEGER, INTENT(IN)                      :: dest
      CLASS(mp_comm_type), INTENT(IN) :: comm
      TYPE(mp_request_type), INTENT(out)                     :: request
      INTEGER, INTENT(in), OPTIONAL            :: tag

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_isend_bm3'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: ierr, msglen, my_tag
      LOGICAL                                  :: foo(1)
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
#if !defined(__GNUC__) || __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(msgin))
#endif

      my_tag = 0
      IF (PRESENT(tag)) my_tag = tag

      msglen = SIZE(msgin, 1)*SIZE(msgin, 2)*SIZE(msgin, 3)
      IF (msglen > 0) THEN
         CALL mpi_isend(msgin(1, 1, 1), msglen, MPI_LOGICAL, dest, my_tag, &
                        comm%handle, request%handle, ierr)
      ELSE
         CALL mpi_isend(foo, msglen, MPI_LOGICAL, dest, my_tag, &
                        comm%handle, request%handle, ierr)
      END IF
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_isend @ "//routineN)

      CALL add_perf(perf_id=11, count=1, msg_size=msglen*loglen)
#else
      CPABORT("mp_isend called in non parallel case")
      MARK_USED(msgin)
      MARK_USED(dest)
      MARK_USED(comm)
      MARK_USED(tag)
      request = mp_request_null
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_isend_bm3

! **************************************************************************************************
!> \brief Non-blocking receive of rank-3 logical data
!> \param msgout the received message
!> \param source the source processor
!> \param comm  the communicator object
!> \param request communication request index
!> \param tag message tag
!> \par History
!>      2.2016 added _bm3 subroutine [Nico Holmberg]
!> \author fawzi
!> \note see mp_irecv_iv
!> \note
!>      arrays can be pointers or assumed shape, but they must be contiguous!
! **************************************************************************************************
   SUBROUTINE mp_irecv_bm3(msgout, source, comm, request, tag)
      LOGICAL, DIMENSION(:, :, :), INTENT(INOUT)  :: msgout
      INTEGER, INTENT(IN)                      :: source
      CLASS(mp_comm_type), INTENT(IN) :: comm
      TYPE(mp_request_type), INTENT(out)                     :: request
      INTEGER, INTENT(in), OPTIONAL            :: tag

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_irecv_bm3'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: ierr, msglen, my_tag
      LOGICAL                                  :: foo(1)
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
#if !defined(__GNUC__) || __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(msgout))
#endif

      my_tag = 0
      IF (PRESENT(tag)) my_tag = tag

      msglen = SIZE(msgout, 1)*SIZE(msgout, 2)*SIZE(msgout, 3)
      IF (msglen > 0) THEN
         CALL mpi_irecv(msgout(1, 1, 1), msglen, MPI_LOGICAL, source, my_tag, &
                        comm%handle, request%handle, ierr)
      ELSE
         CALL mpi_irecv(foo, msglen, MPI_LOGICAL, source, my_tag, &
                        comm%handle, request%handle, ierr)
      END IF
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_ircv @ "//routineN)

      CALL add_perf(perf_id=12, count=1, msg_size=msglen*loglen)
#else
      CPABORT("mp_irecv called in non parallel case")
      MARK_USED(msgout)
      MARK_USED(source)
      MARK_USED(comm)
      MARK_USED(request)
      MARK_USED(tag)
      request = mp_request_null
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_irecv_bm3

! **************************************************************************************************
!> \brief ...
!> \param msg ...
!> \param source ...
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE mp_bcast_av(msg, source, comm)
      CHARACTER(LEN=*), INTENT(INOUT)          :: msg
      INTEGER, INTENT(IN)                                  :: source
      CLASS(mp_comm_type), INTENT(IN) :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_av'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: i, ierr, msglen, numtask, taskid
      INTEGER, DIMENSION(:), ALLOCATABLE       :: imsg
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)

      CALL comm%environ(numtask, taskid)
      IF (taskid == source) msglen = LEN_TRIM(msg)

      CALL comm%bcast(msglen, source)
      ! this is a workaround to avoid problems on the T3E
      ! at the moment we have a data alignment error when trying to
      ! broadcast characters on the T3E (not always!)
      ! JH 19/3/99 on galileo
      ! CALL mpi_bcast(msg,msglen,MPI_CHARACTER,source,gid%handle,ierr)
      ALLOCATE (imsg(1:msglen))
      DO i = 1, msglen
         imsg(i) = ICHAR(msg(i:i))
      END DO
      CALL mpi_bcast(imsg, msglen, MPI_INTEGER, source, comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_bcast @ "//routineN)
      msg = ""
      DO i = 1, msglen
         msg(i:i) = CHAR(imsg(i))
      END DO
      DEALLOCATE (imsg)
      CALL add_perf(perf_id=2, count=1, msg_size=msglen*charlen)
#else
      MARK_USED(msg)
      MARK_USED(source)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_bcast_av

! **************************************************************************************************
!> \brief ...
!> \param msg ...
!> \param source ...
!> \param comm ...
! **************************************************************************************************
   SUBROUTINE mp_bcast_am(msg, source, comm)
      CHARACTER(LEN=*), CONTIGUOUS, INTENT(INOUT)      :: msg(:)
      INTEGER, INTENT(IN)                             :: source
      CLASS(mp_comm_type), INTENT(IN) :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_bcast_am'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: i, ierr, j, k, msglen, msgsiz, &
                                                  numtask, taskid
      INTEGER, ALLOCATABLE                     :: imsg(:), imsglen(:)
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      CALL comm%environ(numtask, taskid)
      msgsiz = SIZE(msg)
      ! Determine size of the minimum array of integers to broadcast the string
      ALLOCATE (imsglen(1:msgsiz))
      DO j = 1, msgsiz
         IF (taskid == source) imsglen(j) = LEN_TRIM(msg(j))
      END DO
      CALL comm%bcast(imsglen, source)
      msglen = SUM(imsglen)
      ! this is a workaround to avoid problems on the T3E
      ! at the moment we have a data alignment error when trying to
      ! broadcast characters on the T3E (not always!)
      ! JH 19/3/99 on galileo
      ALLOCATE (imsg(1:msglen))
      k = 0
      DO j = 1, msgsiz
         DO i = 1, imsglen(j)
            k = k + 1
            imsg(k) = ICHAR(msg(j) (i:i))
         END DO
      END DO
      CALL mpi_bcast(imsg, msglen, MPI_INTEGER, source, comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_bcast @ "//routineN)
      msg = ""
      k = 0
      DO j = 1, msgsiz
         DO i = 1, imsglen(j)
            k = k + 1
            msg(j) (i:i) = CHAR(imsg(k))
         END DO
      END DO
      DEALLOCATE (imsg)
      DEALLOCATE (imsglen)
      CALL add_perf(perf_id=2, count=1, msg_size=msglen*charlen*msgsiz)
#else
      MARK_USED(msg)
      MARK_USED(source)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_bcast_am

! **************************************************************************************************
!> \brief Finds the location of the minimal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] comm             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MINLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! **************************************************************************************************
   SUBROUTINE mp_minloc_dv(msg, comm)
      REAL(kind=real_8), CONTIGUOUS, INTENT(INOUT)         :: msg(:)
      CLASS(mp_comm_type), INTENT(IN)                      :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_minloc_dv'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: ierr, msglen
      REAL(kind=real_8), ALLOCATABLE           :: res(:)
#endif

      IF ("d" .EQ. "l" .AND. real_8 .EQ. int_8) THEN
         CPABORT("Minimal location not available with long integers @ "//routineN)
      END IF
      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      msglen = SIZE(msg)
      ALLOCATE (res(1:msglen), STAT=ierr)
      IF (ierr /= 0) &
         CPABORT("allocate @ "//routineN)
      CALL mpi_allreduce(msg, res, msglen/2, MPI_2DOUBLE_PRECISION, MPI_MINLOC, comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_allreduce @ "//routineN)
      msg = res
      DEALLOCATE (res)
      CALL add_perf(perf_id=3, count=1, msg_size=msglen*real_8_size)
#else
      MARK_USED(msg)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_minloc_dv

! **************************************************************************************************
!> \brief Finds the location of the minimal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] comm             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MINLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! **************************************************************************************************
   SUBROUTINE mp_minloc_iv(msg, comm)
      INTEGER(KIND=int_4), CONTIGUOUS, INTENT(INOUT)       :: msg(:)
      CLASS(mp_comm_type), INTENT(IN)                      :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_minloc_iv'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: ierr, msglen
      INTEGER(KIND=int_4), ALLOCATABLE         :: res(:)
#endif

      IF ("i" .EQ. "l" .AND. int_4 .EQ. int_8) THEN
         CPABORT("Minimal location not available with long integers @ "//routineN)
      END IF
      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      msglen = SIZE(msg)
      ALLOCATE (res(1:msglen))
      CALL mpi_allreduce(msg, res, msglen/2, MPI_2INTEGER, MPI_MINLOC, comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_allreduce @ "//routineN)
      msg = res
      DEALLOCATE (res)
      CALL add_perf(perf_id=3, count=1, msg_size=msglen*int_4_size)
#else
      MARK_USED(msg)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_minloc_iv

! **************************************************************************************************
!> \brief Finds the location of the minimal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] comm             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MINLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! **************************************************************************************************
   SUBROUTINE mp_minloc_lv(msg, comm)
      INTEGER(KIND=int_8), CONTIGUOUS, INTENT(INOUT)       :: msg(:)
      CLASS(mp_comm_type), INTENT(IN)                      :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_minloc_lv'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: ierr, msglen
      INTEGER(KIND=int_8), ALLOCATABLE         :: res(:)
#endif

      IF ("l" .EQ. "l" .AND. int_8 .EQ. int_8) THEN
         CPABORT("Minimal location not available with long integers @ "//routineN)
      END IF
      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      msglen = SIZE(msg)
      ALLOCATE (res(1:msglen))
      CALL mpi_allreduce(msg, res, msglen/2, MPI_INTEGER8, MPI_MINLOC, comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_allreduce @ "//routineN)
      msg = res
      DEALLOCATE (res)
      CALL add_perf(perf_id=3, count=1, msg_size=msglen*int_8_size)
#else
      MARK_USED(msg)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_minloc_lv

! **************************************************************************************************
!> \brief Finds the location of the minimal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] comm             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MINLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! **************************************************************************************************
   SUBROUTINE mp_minloc_rv(msg, comm)
      REAL(kind=real_4), CONTIGUOUS, INTENT(INOUT)         :: msg(:)
      CLASS(mp_comm_type), INTENT(IN)                      :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_minloc_rv'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: ierr, msglen
      REAL(kind=real_4), ALLOCATABLE           :: res(:)
#endif

      IF ("r" .EQ. "l" .AND. real_4 .EQ. int_8) THEN
         CPABORT("Minimal location not available with long integers @ "//routineN)
      END IF
      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      msglen = SIZE(msg)
      ALLOCATE (res(1:msglen))
      CALL mpi_allreduce(msg, res, msglen/2, MPI_2REAL, MPI_MINLOC, comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_allreduce @ "//routineN)
      msg = res
      DEALLOCATE (res)
      CALL add_perf(perf_id=3, count=1, msg_size=msglen*real_4_size)
#else
      MARK_USED(msg)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_minloc_rv

! **************************************************************************************************
!> \brief Finds the location of the maximal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] comm             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MAXLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! **************************************************************************************************
   SUBROUTINE mp_maxloc_dv(msg, comm)
      REAL(kind=real_8), CONTIGUOUS, INTENT(INOUT)         :: msg(:)
      CLASS(mp_comm_type), INTENT(IN)                      :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_maxloc_dv'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: ierr, msglen
      REAL(kind=real_8), ALLOCATABLE           :: res(:)
#endif

      IF ("d" .EQ. "l" .AND. real_8 .EQ. int_8) THEN
         CPABORT("Maximal location not available with long integers @ "//routineN)
      END IF
      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      msglen = SIZE(msg)
      ALLOCATE (res(1:msglen))
      CALL mpi_allreduce(msg, res, msglen/2, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_allreduce @ "//routineN)
      msg = res
      DEALLOCATE (res)
      CALL add_perf(perf_id=3, count=1, msg_size=msglen*real_8_size)
#else
      MARK_USED(msg)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_maxloc_dv

! **************************************************************************************************
!> \brief Finds the location of the maximal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] comm             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MAXLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! **************************************************************************************************
   SUBROUTINE mp_maxloc_iv(msg, comm)
      INTEGER(KIND=int_4), CONTIGUOUS, INTENT(INOUT)       :: msg(:)
      CLASS(mp_comm_type), INTENT(IN)                      :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_maxloc_iv'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: ierr, msglen
      INTEGER(KIND=int_4), ALLOCATABLE         :: res(:)
#endif

      IF ("i" .EQ. "l" .AND. int_4 .EQ. int_8) THEN
         CPABORT("Maximal location not available with long integers @ "//routineN)
      END IF
      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      msglen = SIZE(msg)
      ALLOCATE (res(1:msglen))
      CALL mpi_allreduce(msg, res, msglen/2, MPI_2INTEGER, MPI_MAXLOC, comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_allreduce @ "//routineN)
      msg = res
      DEALLOCATE (res)
      CALL add_perf(perf_id=3, count=1, msg_size=msglen*int_4_size)
#else
      MARK_USED(msg)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_maxloc_iv

! **************************************************************************************************
!> \brief Finds the location of the maximal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] comm             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MAXLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! **************************************************************************************************
   SUBROUTINE mp_maxloc_lv(msg, comm)
      INTEGER(KIND=int_8), CONTIGUOUS, INTENT(INOUT)       :: msg(:)
      CLASS(mp_comm_type), INTENT(IN)                      :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_maxloc_lv'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: ierr, msglen
      INTEGER(KIND=int_8), ALLOCATABLE         :: res(:)
#endif

      IF ("l" .EQ. "l" .AND. int_8 .EQ. int_8) THEN
         CPABORT("Maximal location not available with long integers @ "//routineN)
      END IF
      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      msglen = SIZE(msg)
      ALLOCATE (res(1:msglen))
      CALL mpi_allreduce(msg, res, msglen/2, MPI_INTEGER8, MPI_MAXLOC, comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_allreduce @ "//routineN)
      msg = res
      DEALLOCATE (res)
      CALL add_perf(perf_id=3, count=1, msg_size=msglen*int_8_size)
#else
      MARK_USED(msg)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_maxloc_lv

! **************************************************************************************************
!> \brief Finds the location of the maximal element in a vector.
!> \param[in,out] msg         Find location of maximum element among these
!>                            data (input).
!> \param[in] comm             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce with the MPI_MAXLOC reduction function identifier
!> \par Invalid data types
!>      This routine is invalid for (int_8) data!
! **************************************************************************************************
   SUBROUTINE mp_maxloc_rv(msg, comm)
      REAL(kind=real_4), CONTIGUOUS, INTENT(INOUT)         :: msg(:)
      CLASS(mp_comm_type), INTENT(IN)                      :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_maxloc_rv'

      INTEGER                                  :: handle
#if defined(__parallel)
      INTEGER                                  :: ierr, msglen
      REAL(kind=real_4), ALLOCATABLE           :: res(:)
#endif

      IF ("r" .EQ. "l" .AND. real_4 .EQ. int_8) THEN
         CPABORT("Maximal location not available with long integers @ "//routineN)
      END IF
      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      msglen = SIZE(msg)
      ALLOCATE (res(1:msglen))
      CALL mpi_allreduce(msg, res, msglen/2, MPI_2REAL, MPI_MAXLOC, comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_allreduce @ "//routineN)
      msg = res
      DEALLOCATE (res)
      CALL add_perf(perf_id=3, count=1, msg_size=msglen*real_4_size)
#else
      MARK_USED(msg)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_maxloc_rv

! **************************************************************************************************
!> \brief Logical OR reduction
!> \param[in,out] msg         Datum to perform inclusive disjunction (input)
!>                            and resultant inclusive disjunction (output)
!> \param[in] comm             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce
! **************************************************************************************************
   SUBROUTINE mp_sum_b(msg, comm)
      LOGICAL, INTENT(INOUT)                             :: msg
      CLASS(mp_comm_type), INTENT(IN)                                :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_b'

      INTEGER                                            :: handle
#if defined(__parallel)
      INTEGER :: ierr, msglen
#endif

      CALL mp_timeset(routineN, handle)
#if defined(__parallel)
      msglen = 1
      CALL mpi_allreduce(MPI_IN_PLACE, msg, msglen, MPI_LOGICAL, MPI_LOR, comm%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_allreduce @ "//routineN)
#else
      MARK_USED(msg)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_sum_b

! **************************************************************************************************
!> \brief Logical OR reduction
!> \param[in,out] msg         Datum to perform inclusive disjunction (input)
!>                            and resultant inclusive disjunction (output)
!> \param[in] comm             Message passing environment identifier
!> \par MPI mapping
!>      mpi_allreduce
! **************************************************************************************************
   SUBROUTINE mp_sum_bv(msg, comm)
      LOGICAL, DIMENSION(:), CONTIGUOUS, INTENT(INOUT)               :: msg
      CLASS(mp_comm_type), INTENT(IN)                                :: comm

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_sum_bv'

      INTEGER                                            :: handle
#if defined(__parallel)
      INTEGER :: ierr, msglen
#endif

      CALL mp_timeset(routineN, handle)
#if defined(__parallel)
      msglen = SIZE(msg)
      IF (msglen .GT. 0) THEN
         CALL mpi_allreduce(MPI_IN_PLACE, msg, msglen, MPI_LOGICAL, MPI_LOR, comm%handle, ierr)
         IF (ierr /= 0) CALL mp_stop(ierr, "mpi_allreduce @ "//routineN)
      END IF
#else
      MARK_USED(msg)
      MARK_USED(comm)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_sum_bv

! **************************************************************************************************
!> \brief Logical OR reduction
!> \param[in,out] msg         Datum to perform inclusive disjunction (input)
!>                            and resultant inclusive disjunction (output)
!> \param[in] comm             Message passing environment identifier
!> \param request ...
!> \par MPI mapping
!>      mpi_allreduce
! **************************************************************************************************
   SUBROUTINE mp_isum_bv(msg, comm, request)
      LOGICAL, DIMENSION(:), INTENT(INOUT)               :: msg
      CLASS(mp_comm_type), INTENT(IN)                                :: comm
      TYPE(mp_request_type), INTENT(INOUT)                             :: request

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_isum_bv'

      INTEGER                                            :: handle
#if defined(__parallel)
      INTEGER :: ierr, msglen
#endif

      CALL mp_timeset(routineN, handle)
#if defined(__parallel)
      msglen = SIZE(msg)
#if !defined(__GNUC__) || __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(msg))
#endif

      IF (msglen .GT. 0) THEN
         CALL mpi_iallreduce(MPI_IN_PLACE, msg, msglen, MPI_LOGICAL, MPI_LOR, comm%handle, request%handle, ierr)
         IF (ierr /= 0) CALL mp_stop(ierr, "mpi_allreduce @ "//routineN)
      ELSE
         request = mp_request_null
      END IF
#else
      MARK_USED(msg)
      MARK_USED(comm)
      request = mp_request_null
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_isum_bv

! **************************************************************************************************
!> \brief Get Version of the MPI Library (MPI 3)
!> \param[out] version        Version of the library,
!>                            declared as CHARACTER(LEN=mp_max_library_version_string)
!> \param[out] resultlen      Length (in printable characters) of
!>                            the result returned in version (integer)
! **************************************************************************************************
   SUBROUTINE mp_get_library_version(version, resultlen)
      CHARACTER(len=*), INTENT(OUT)                      :: version
      INTEGER, INTENT(OUT)                               :: resultlen

#if defined(__parallel)
      INTEGER                                            :: ierr
#endif

      version = ''

#if defined(__parallel)
      ierr = 0
      CALL mpi_get_library_version(version, resultlen, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_get_library_version @ mp_get_library_version")
#else
      resultlen = 0
#endif
   END SUBROUTINE mp_get_library_version

! **************************************************************************************************
!> \brief Opens a file
!> \param[in] groupid    message passing environment identifier
!> \param[out] fh        file handle (file storage unit)
!> \param[in] filepath   path to the file
!> \param amode_status   access mode
!> \param info ...
!> \par MPI-I/O mapping  mpi_file_open
!> \par STREAM-I/O mapping  OPEN
!>
!> \param[in](optional) info   info object
!> \par History
!>      11.2012 created [Hossein Bani-Hashemian]
! **************************************************************************************************
   SUBROUTINE mp_file_open(groupid, fh, filepath, amode_status, info)
      TYPE(mp_comm_type), INTENT(IN)                      :: groupid
      CLASS(mp_file_type), INTENT(OUT)                     :: fh
      CHARACTER(len=*), INTENT(IN)             :: filepath
      INTEGER, INTENT(IN)                      :: amode_status
      TYPE(mp_info_type), INTENT(IN), OPTIONAL            :: info

#if defined(__parallel)
      INTEGER                                  :: ierr
      MPI_INFO_TYPE                           :: my_info
#else
      CHARACTER(LEN=10)                        :: fstatus, fposition
      INTEGER                                  :: amode, handle, istat
      LOGICAL                                  :: exists, is_open
#endif

#if defined(__parallel)
      ierr = 0
      my_info = mpi_info_null
      IF (PRESENT(info)) my_info = info%handle
      CALL mpi_file_open(groupid%handle, filepath, amode_status, my_info, fh%handle, ierr)
      CALL mpi_file_set_errhandler(fh%handle, MPI_ERRORS_RETURN, ierr)
      IF (ierr .NE. 0) CALL mp_stop(ierr, "mpi_file_set_errhandler @ mp_file_open")
#else
      MARK_USED(groupid)
      MARK_USED(info)
      amode = amode_status
      IF (amode .GT. file_amode_append) THEN
         fposition = "APPEND"
         amode = amode - file_amode_append
      ELSE
         fposition = "REWIND"
      END IF
      IF ((amode .EQ. file_amode_create) .OR. &
          (amode .EQ. file_amode_create + file_amode_wronly) .OR. &
          (amode .EQ. file_amode_create + file_amode_wronly + file_amode_excl)) THEN
         fstatus = "UNKNOWN"
      ELSE
         fstatus = "OLD"
      END IF
      ! Get a new unit number
      DO handle = 1, 999
         INQUIRE (UNIT=handle, EXIST=exists, OPENED=is_open, IOSTAT=istat)
         IF (exists .AND. (.NOT. is_open) .AND. (istat == 0)) EXIT
      END DO
      OPEN (UNIT=handle, FILE=filepath, STATUS=fstatus, ACCESS="STREAM", POSITION=fposition)
      fh%handle = handle
#endif
   END SUBROUTINE mp_file_open

! **************************************************************************************************
!> \brief Deletes a file. Auxiliary routine to emulate 'replace' action for mp_file_open.
!>        Only the master processor should call this routine.
!> \param[in] filepath   path to the file
!> \param[in](optional) info   info object
!> \par History
!>      11.2017 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mp_file_delete(filepath, info)
      CHARACTER(len=*), INTENT(IN)             :: filepath
      TYPE(mp_info_type), INTENT(IN), OPTIONAL            :: info

#if defined(__parallel)
      INTEGER                                  :: ierr
      MPI_INFO_TYPE                           :: my_info
      LOGICAL                                  :: exists

      ierr = 0
      my_info = mpi_info_null
      IF (PRESENT(info)) my_info = info%handle
      INQUIRE (FILE=filepath, EXIST=exists)
      IF (exists) CALL mpi_file_delete(filepath, my_info, ierr)
      IF (ierr .NE. 0) CALL mp_stop(ierr, "mpi_file_set_errhandler @ mp_file_delete")
#else
      MARK_USED(filepath)
      MARK_USED(info)
      ! Explicit file delete not necessary, handled by subsequent call to open_file with action 'replace'
#endif

   END SUBROUTINE mp_file_delete

! **************************************************************************************************
!> \brief Closes a file
!> \param[in] fh   file handle (file storage unit)
!> \par MPI-I/O mapping   mpi_file_close
!> \par STREAM-I/O mapping   CLOSE
!>
!> \par History
!>      11.2012 created [Hossein Bani-Hashemian]
! **************************************************************************************************
   SUBROUTINE mp_file_close(fh)
      CLASS(mp_file_type), INTENT(INOUT)                             :: fh

#if defined(__parallel)
      INTEGER                                            :: ierr

      ierr = 0
#endif
      IF (fh%handle == mp_file_null_handle) &
         CPABORT("Invalid file handle")
#if defined(__parallel)
      CALL mpi_file_set_errhandler(fh%handle, MPI_ERRORS_RETURN, ierr)
      CALL mpi_file_close(fh%handle, ierr)
      IF (ierr .NE. 0) CALL mp_stop(ierr, "mpi_file_set_errhandler @ mp_file_close")
#else
      CLOSE (fh%handle)
      fh%handle = mp_file_null_handle
#endif
   END SUBROUTINE mp_file_close

   SUBROUTINE mp_file_assign(fh_new, fh_old)
      CLASS(mp_file_type), INTENT(OUT) :: fh_new
      CLASS(mp_file_type), INTENT(IN) :: fh_old

      fh_new%handle = fh_old%handle

   END SUBROUTINE

! **************************************************************************************************
!> \brief Returns the file size
!> \param[in] fh file handle (file storage unit)
!> \param[out] file_size  the file size
!> \par MPI-I/O mapping   mpi_file_get_size
!> \par STREAM-I/O mapping   INQUIRE
!>
!> \par History
!>      12.2012 created [Hossein Bani-Hashemian]
! **************************************************************************************************
   SUBROUTINE mp_file_get_size(fh, file_size)
      CLASS(mp_file_type), INTENT(IN)                                :: fh
      INTEGER(kind=file_offset), INTENT(OUT)             :: file_size

#if defined(__parallel)
      INTEGER                                            :: ierr
#endif

      IF (fh%handle == mp_file_null_handle) &
         CPABORT("Invalid file handle")

#if defined(__parallel)
      ierr = 0
      CALL mpi_file_set_errhandler(fh%handle, MPI_ERRORS_RETURN, ierr)
      CALL mpi_file_get_size(fh%handle, file_size, ierr)
      IF (ierr .NE. 0) CALL mp_stop(ierr, "mpi_file_set_errhandler @ mp_file_get_size")
#else
      INQUIRE (UNIT=fh%handle, SIZE=file_size)
#endif
   END SUBROUTINE mp_file_get_size

! **************************************************************************************************
!> \brief Returns the file position
!> \param[in] fh file handle (file storage unit)
!> \param[out] file_size  the file position
!> \par MPI-I/O mapping   mpi_file_get_position
!> \par STREAM-I/O mapping   INQUIRE
!>
!> \par History
!>      11.2017 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE mp_file_get_position(fh, pos)
      CLASS(mp_file_type), INTENT(IN)                                :: fh
      INTEGER(kind=file_offset), INTENT(OUT)             :: pos

#if defined(__parallel)
      INTEGER                                            :: ierr
#endif

      IF (fh%handle == mp_file_null_handle) &
         CPABORT("Invalid file handle")

#if defined(__parallel)
      ierr = 0
      CALL mpi_file_set_errhandler(fh%handle, MPI_ERRORS_RETURN, ierr)
      CALL mpi_file_get_position(fh%handle, pos, ierr)
      IF (ierr .NE. 0) CALL mp_stop(ierr, "mpi_file_set_errhandler @ mp_file_get_position")
#else
      INQUIRE (UNIT=fh%handle, POS=pos)
#endif
   END SUBROUTINE mp_file_get_position

! **************************************************************************************************
!> \brief (parallel) Blocking individual file write using explicit offsets
!>        (serial) Unformatted stream write
!> \param[in] fh     file handle (file storage unit)
!> \param[in] offset file offset (position)
!> \param[in] msg    data to be written to the file
!> \param msglen ...
!> \par MPI-I/O mapping   mpi_file_write_at
!> \par STREAM-I/O mapping   WRITE
!> \param[in](optional) msglen number of the elements of data
! **************************************************************************************************
   SUBROUTINE mp_file_write_at_chv(fh, offset, msg, msglen)
      CHARACTER, CONTIGUOUS, INTENT(IN)                      :: msg(:)
      CLASS(mp_file_type), INTENT(IN)                        :: fh
      INTEGER, INTENT(IN), OPTIONAL              :: msglen
      INTEGER(kind=file_offset), INTENT(IN)      :: offset

#if defined(__parallel)
      INTEGER                                    :: ierr, msg_len
#endif

      IF (fh%handle == mp_file_null_handle) &
         CPABORT("Invalid file handle")

#if defined(__parallel)
      msg_len = SIZE(msg)
      IF (PRESENT(msglen)) msg_len = msglen
      CALL MPI_FILE_WRITE_AT(fh%handle, offset, msg, msg_len, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      IF (ierr .NE. 0) &
         CPABORT("mpi_file_write_at_chv @ mp_file_write_at_chv")
#else
      MARK_USED(msglen)
      WRITE (UNIT=fh%handle, POS=offset + 1) msg
#endif
   END SUBROUTINE mp_file_write_at_chv

! **************************************************************************************************
!> \brief ...
!> \param fh ...
!> \param offset ...
!> \param msg ...
! **************************************************************************************************
   SUBROUTINE mp_file_write_at_ch(fh, offset, msg)
      CHARACTER(LEN=*), INTENT(IN)               :: msg
      CLASS(mp_file_type), INTENT(IN)                        :: fh
      INTEGER(kind=file_offset), INTENT(IN)      :: offset

#if defined(__parallel)
      INTEGER                                    :: ierr
#endif

      IF (fh%handle == mp_file_null_handle) &
         CPABORT("Invalid file handle")

#if defined(__parallel)
      CALL MPI_FILE_WRITE_AT(fh%handle, offset, msg, LEN(msg), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      IF (ierr .NE. 0) &
         CPABORT("mpi_file_write_at_ch @ mp_file_write_at_ch")
#else
      WRITE (UNIT=fh%handle, POS=offset + 1) msg
#endif
   END SUBROUTINE mp_file_write_at_ch

! **************************************************************************************************
!> \brief (parallel) Blocking collective file write using explicit offsets
!>        (serial) Unformatted stream write
!> \param fh ...
!> \param offset ...
!> \param msg ...
!> \param msglen ...
!> \par MPI-I/O mapping   mpi_file_write_at_all
!> \par STREAM-I/O mapping   WRITE
! **************************************************************************************************
   SUBROUTINE mp_file_write_at_all_chv(fh, offset, msg, msglen)
      CHARACTER, CONTIGUOUS, INTENT(IN)                      :: msg(:)
      CLASS(mp_file_type), INTENT(IN)                        :: fh
      INTEGER, INTENT(IN), OPTIONAL              :: msglen
      INTEGER(kind=file_offset), INTENT(IN)      :: offset

#if defined(__parallel)
      INTEGER                                    :: ierr, msg_len
#endif

      IF (fh%handle == mp_file_null_handle) &
         CPABORT("Invalid file handle")

#if defined(__parallel)
      msg_len = SIZE(msg)
      IF (PRESENT(msglen)) msg_len = msglen
      CALL MPI_FILE_WRITE_AT_ALL(fh%handle, offset, msg, msg_len, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      IF (ierr .NE. 0) &
         CPABORT("mpi_file_write_at_all_chv @ mp_file_write_at_all_chv")
#else
      MARK_USED(msglen)
      WRITE (UNIT=fh%handle, POS=offset + 1) msg
#endif
   END SUBROUTINE mp_file_write_at_all_chv

! **************************************************************************************************
!> \brief ...
!> \param fh ...
!> \param offset ...
!> \param msg ...
! **************************************************************************************************
   SUBROUTINE mp_file_write_at_all_ch(fh, offset, msg)
      CHARACTER(LEN=*), INTENT(IN)               :: msg
      CLASS(mp_file_type), INTENT(IN)                        :: fh
      INTEGER(kind=file_offset), INTENT(IN)      :: offset

#if defined(__parallel)
      INTEGER                                    :: ierr
#endif

      IF (fh%handle == mp_file_null_handle) &
         CPABORT("Invalid file handle")

#if defined(__parallel)
      CALL MPI_FILE_WRITE_AT_ALL(fh%handle, offset, msg, LEN(msg), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      IF (ierr .NE. 0) &
         CPABORT("mpi_file_write_at_all_ch @ mp_file_write_at_all_ch")
#else
      WRITE (UNIT=fh%handle, POS=offset + 1) msg
#endif
   END SUBROUTINE mp_file_write_at_all_ch

! **************************************************************************************************
!> \brief (parallel) Blocking individual file read using explicit offsets
!>        (serial) Unformatted stream read
!> \param[in] fh     file handle (file storage unit)
!> \param[in] offset file offset (position)
!> \param[out] msg   data to be read from the file
!> \param msglen ...
!> \par MPI-I/O mapping   mpi_file_read_at
!> \par STREAM-I/O mapping   READ
!> \param[in](optional) msglen  number of elements of data
! **************************************************************************************************
   SUBROUTINE mp_file_read_at_chv(fh, offset, msg, msglen)
      CHARACTER, CONTIGUOUS, INTENT(OUT)                     :: msg(:)
      CLASS(mp_file_type), INTENT(IN)                        :: fh
      INTEGER, INTENT(IN), OPTIONAL              :: msglen
      INTEGER(kind=file_offset), INTENT(IN)      :: offset

#if defined(__parallel)
      INTEGER                                    :: ierr, msg_len
#endif

      IF (fh%handle == mp_file_null_handle) &
         CPABORT("Invalid file handle")

#if defined(__parallel)
      msg_len = SIZE(msg)
      IF (PRESENT(msglen)) msg_len = msglen
      CALL MPI_FILE_READ_AT(fh%handle, offset, msg, msg_len, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      IF (ierr .NE. 0) &
         CPABORT("mpi_file_read_at_chv @ mp_file_read_at_chv")
#else
      MARK_USED(msglen)
      READ (UNIT=fh%handle, POS=offset + 1) msg
#endif
   END SUBROUTINE mp_file_read_at_chv

! **************************************************************************************************
!> \brief ...
!> \param fh ...
!> \param offset ...
!> \param msg ...
! **************************************************************************************************
   SUBROUTINE mp_file_read_at_ch(fh, offset, msg)
      CHARACTER(LEN=*), INTENT(OUT)              :: msg
      CLASS(mp_file_type), INTENT(IN)                        :: fh
      INTEGER(kind=file_offset), INTENT(IN)      :: offset

#if defined(__parallel)
      INTEGER                                    :: ierr
#endif

      IF (fh%handle == mp_file_null_handle) &
         CPABORT("Invalid file handle")

#if defined(__parallel)
      CALL MPI_FILE_READ_AT(fh%handle, offset, msg, LEN(msg), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      IF (ierr .NE. 0) &
         CPABORT("mpi_file_read_at_ch @ mp_file_read_at_ch")
#else
      READ (UNIT=fh%handle, POS=offset + 1) msg
#endif
   END SUBROUTINE mp_file_read_at_ch

! **************************************************************************************************
!> \brief (parallel) Blocking collective file read using explicit offsets
!>        (serial) Unformatted stream read
!> \param fh ...
!> \param offset ...
!> \param msg ...
!> \param msglen ...
!> \par MPI-I/O mapping    mpi_file_read_at_all
!> \par STREAM-I/O mapping   READ
! **************************************************************************************************
   SUBROUTINE mp_file_read_at_all_chv(fh, offset, msg, msglen)
      CHARACTER, INTENT(OUT)                     :: msg(:)
      CLASS(mp_file_type), INTENT(IN)                        :: fh
      INTEGER, INTENT(IN), OPTIONAL              :: msglen
      INTEGER(kind=file_offset), INTENT(IN)      :: offset

#if defined(__parallel)
      INTEGER                                    :: ierr, msg_len
#endif

      IF (fh%handle == mp_file_null_handle) &
         CPABORT("Invalid file handle")

#if defined(__parallel)
      msg_len = SIZE(msg)
      IF (PRESENT(msglen)) msg_len = msglen
      CALL MPI_FILE_READ_AT_ALL(fh%handle, offset, msg, msg_len, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      IF (ierr .NE. 0) &
         CPABORT("mpi_file_read_at_all_chv @ mp_file_read_at_all_chv")
#else
      MARK_USED(msglen)
      READ (UNIT=fh%handle, POS=offset + 1) msg
#endif
   END SUBROUTINE mp_file_read_at_all_chv

! **************************************************************************************************
!> \brief ...
!> \param fh ...
!> \param offset ...
!> \param msg ...
! **************************************************************************************************
   SUBROUTINE mp_file_read_at_all_ch(fh, offset, msg)
      CHARACTER(LEN=*), INTENT(OUT)              :: msg
      CLASS(mp_file_type), INTENT(IN)                        :: fh
      INTEGER(kind=file_offset), INTENT(IN)      :: offset

#if defined(__parallel)
      INTEGER                                    :: ierr
#endif

      IF (fh%handle == mp_file_null_handle) &
         CPABORT("Invalid file handle")

#if defined(__parallel)
      CALL MPI_FILE_READ_AT_ALL(fh%handle, offset, msg, LEN(msg), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      IF (ierr .NE. 0) &
         CPABORT("mpi_file_read_at_all_ch @ mp_file_read_at_all_ch")
#else
      READ (UNIT=fh%handle, POS=offset + 1) msg
#endif
   END SUBROUTINE mp_file_read_at_all_ch

! **************************************************************************************************
!> \brief Returns the size of a data type in bytes
!> \param[in] type_descriptor  data type
!> \param[out] type_size       size of the data type
!> \par MPI mapping
!>      mpi_type_size
!>
! **************************************************************************************************
   SUBROUTINE mp_type_size(type_descriptor, type_size)
      TYPE(mp_type_descriptor_type), INTENT(IN)          :: type_descriptor
      INTEGER, INTENT(OUT)                               :: type_size

#if defined(__parallel)
      INTEGER                                            :: ierr

      ierr = 0
      CALL MPI_TYPE_SIZE(type_descriptor%type_handle, type_size, ierr)
      IF (ierr .NE. 0) &
         CPABORT("mpi_type_size failed @ mp_type_size")
#else
      SELECT CASE (type_descriptor%type_handle)
      CASE (1)
         type_size = real_4_size
      CASE (3)
         type_size = real_8_size
      CASE (5)
         type_size = 2*real_4_size
      CASE (7)
         type_size = 2*real_8_size
      END SELECT
#endif
   END SUBROUTINE mp_type_size

! **************************************************************************************************
!> \brief ...
!> \param subtypes ...
!> \param vector_descriptor ...
!> \param index_descriptor ...
!> \return ...
! **************************************************************************************************
   FUNCTION mp_type_make_struct(subtypes, &
                                vector_descriptor, index_descriptor) &
      RESULT(type_descriptor)
      TYPE(mp_type_descriptor_type), &
         DIMENSION(:), INTENT(IN)               :: subtypes
      INTEGER, DIMENSION(2), INTENT(IN), &
         OPTIONAL                               :: vector_descriptor
      TYPE(mp_indexing_meta_type), &
         INTENT(IN), OPTIONAL                   :: index_descriptor
      TYPE(mp_type_descriptor_type)            :: type_descriptor

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_type_make_struct'

      INTEGER                                  :: i, n
      INTEGER, ALLOCATABLE, DIMENSION(:)       :: lengths
#if defined(__parallel)
      INTEGER :: ierr
      INTEGER(kind=mpi_address_kind), &
         ALLOCATABLE, DIMENSION(:)              :: displacements
#endif
      MPI_DATA_TYPE, ALLOCATABLE, DIMENSION(:) :: old_types

      n = SIZE(subtypes)
      !type_descriptor%mpi_type_handle = MPI_DATATYPE_NULL
      type_descriptor%length = 1
#if defined(__parallel)
      ierr = 0
      CALL mpi_get_address(MPI_BOTTOM, type_descriptor%base, ierr)
      IF (ierr /= 0) &
         CPABORT("MPI_get_address @ "//routineN)
      ALLOCATE (displacements(n))
#endif
      type_descriptor%vector_descriptor(1:2) = 1
      type_descriptor%has_indexing = .FALSE.
      ALLOCATE (type_descriptor%subtype(n))
      type_descriptor%subtype(:) = subtypes(:)
      ALLOCATE (lengths(n), old_types(n))
      DO i = 1, SIZE(subtypes)
#if defined(__parallel)
         displacements(i) = subtypes(i)%base
#endif
         old_types(i) = subtypes(i)%type_handle
         lengths(i) = subtypes(i)%length
      END DO
#if defined(__parallel)
      CALL MPI_Type_create_struct(n, &
                                  lengths, displacements, old_types, &
                                  type_descriptor%type_handle, ierr)
      IF (ierr /= 0) &
         CPABORT("MPI_Type_create_struct @ "//routineN)
      CALL MPI_Type_commit(type_descriptor%type_handle, ierr)
      IF (ierr /= 0) &
         CPABORT("MPI_Type_commit @ "//routineN)
#endif
      IF (PRESENT(vector_descriptor) .OR. PRESENT(index_descriptor)) THEN
         CPABORT(routineN//" Vectors and indices NYI")
      END IF
   END FUNCTION mp_type_make_struct

! **************************************************************************************************
!> \brief ...
!> \param type_descriptor ...
! **************************************************************************************************
   RECURSIVE SUBROUTINE mp_type_free_m(type_descriptor)
      TYPE(mp_type_descriptor_type), INTENT(inout)       :: type_descriptor

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_type_free_m'

      INTEGER                                            :: handle, i
#if defined(__parallel)
      INTEGER :: ierr
#endif

      CALL mp_timeset(routineN, handle)

      ! If the subtype is associated, then it's a user-defined data type.

      IF (ASSOCIATED(type_descriptor%subtype)) THEN
         DO i = 1, SIZE(type_descriptor%subtype)
            CALL mp_type_free_m(type_descriptor%subtype(i))
         END DO
         DEALLOCATE (type_descriptor%subtype)
      END IF
#if defined(__parallel)
      ierr = 0
      CALL MPI_Type_free(type_descriptor%type_handle, ierr)
      IF (ierr /= 0) &
         CPABORT("MPI_Type_free @ "//routineN)
#endif

      CALL mp_timestop(handle)

   END SUBROUTINE mp_type_free_m

! **************************************************************************************************
!> \brief ...
!> \param type_descriptors ...
! **************************************************************************************************
   SUBROUTINE mp_type_free_v(type_descriptors)
      TYPE(mp_type_descriptor_type), DIMENSION(:), &
         INTENT(inout)                                   :: type_descriptors

      INTEGER                                            :: i

      DO i = 1, SIZE(type_descriptors)
         CALL mp_type_free(type_descriptors(i))
      END DO

   END SUBROUTINE mp_type_free_v

! **************************************************************************************************
!> \brief Creates an indexed MPI type for arrays of strings using bytes for spacing (hindexed type)
!> \param count   number of array blocks to read
!> \param lengths lengths of each array block
!> \param displs  byte offsets for array blocks
!> \return container holding the created type
!> \author Nico Holmberg [05.2017]
! **************************************************************************************************
   FUNCTION mp_file_type_hindexed_make_chv(count, lengths, displs) &
      RESULT(type_descriptor)
      INTEGER, INTENT(IN)                       :: count
      INTEGER, DIMENSION(1:count), &
         INTENT(IN), TARGET                     :: lengths
      INTEGER(kind=file_offset), &
         DIMENSION(1:count), INTENT(in), TARGET :: displs
      TYPE(mp_file_descriptor_type)             :: type_descriptor

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_hindexed_make_chv'

      INTEGER :: ierr, handle

      ierr = 0
      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      CALL MPI_Type_create_hindexed(count, lengths, INT(displs, KIND=address_kind), MPI_CHARACTER, &
                                    type_descriptor%type_handle, ierr)
      IF (ierr /= 0) &
         CPABORT("MPI_Type_create_hindexed @ "//routineN)
      CALL MPI_Type_commit(type_descriptor%type_handle, ierr)
      IF (ierr /= 0) &
         CPABORT("MPI_Type_commit @ "//routineN)
#else
      type_descriptor%type_handle = 68
#endif
      type_descriptor%length = count
      type_descriptor%has_indexing = .TRUE.
      type_descriptor%index_descriptor%index => lengths
      type_descriptor%index_descriptor%chunks => displs

      CALL mp_timestop(handle)

   END FUNCTION mp_file_type_hindexed_make_chv

! **************************************************************************************************
!> \brief Uses a previously created indexed MPI character type to tell the MPI processes
!>        how to partition (set_view) an opened file
!> \param fh      the file handle associated with the input file
!> \param offset  global offset determining where the relevant data begins
!> \param type_descriptor container for the MPI type
!> \author Nico Holmberg [05.2017]
! **************************************************************************************************
   SUBROUTINE mp_file_type_set_view_chv(fh, offset, type_descriptor)
      TYPE(mp_file_type), INTENT(IN)                      :: fh
      INTEGER(kind=file_offset), INTENT(IN)    :: offset
      TYPE(mp_file_descriptor_type)            :: type_descriptor

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_set_view_chv'

      INTEGER                                   :: handle
#if defined(__parallel)
      INTEGER :: ierr
#endif

      CALL mp_timeset(routineN, handle)

      IF (fh%handle == mp_file_null_handle) &
         CPABORT("Invalid file handle")

#if defined(__parallel)
      ierr = 0
      CALL mpi_file_set_errhandler(fh%handle, MPI_ERRORS_RETURN, ierr)
      CALL MPI_File_set_view(fh%handle, offset, MPI_CHARACTER, &
                             type_descriptor%type_handle, "native", MPI_INFO_NULL, ierr)
      IF (ierr .NE. 0) CALL mp_stop(ierr, "mpi_file_set_errhandler @ MPI_File_set_view")
#else
      ! Uses absolute offsets stored in mp_file_descriptor_type
      MARK_USED(fh)
      MARK_USED(offset)
      MARK_USED(type_descriptor)
#endif

      CALL mp_timestop(handle)

   END SUBROUTINE mp_file_type_set_view_chv

! **************************************************************************************************
!> \brief (parallel) Collective, blocking read of a character array from a file. File access pattern
!                    determined by a previously set file view.
!>        (serial)   Unformatted stream read using explicit offsets
!> \param fh     the file handle associated with the input file
!> \param msglen the message length of an individual vector component
!> \param ndims  the number of vector components
!> \param buffer the buffer where the data is placed
!> \param type_descriptor container for the MPI type
!> \author Nico Holmberg [05.2017]
! **************************************************************************************************
   SUBROUTINE mp_file_read_all_chv(fh, msglen, ndims, buffer, type_descriptor)
      CLASS(mp_file_type), INTENT(IN)                       :: fh
      INTEGER, INTENT(IN)                       :: msglen
      INTEGER, INTENT(IN)                       :: ndims
      CHARACTER(LEN=msglen), DIMENSION(ndims), INTENT(INOUT)   :: buffer
      TYPE(mp_file_descriptor_type), &
         INTENT(IN), OPTIONAL                   :: type_descriptor

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_read_all_chv'

      INTEGER                                   :: handle
#if defined(__parallel)
      INTEGER:: ierr
#else
      INTEGER :: i
#endif

      CALL mp_timeset(routineN, handle)

      IF (fh%handle == mp_file_null_handle) &
         CPABORT("Invalid file handle")

#if defined(__parallel)
      ierr = 0
      MARK_USED(type_descriptor)
      CALL MPI_File_read_all(fh%handle, buffer, ndims*msglen, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      IF (ierr .NE. 0) CALL mp_stop(ierr, "mpi_file_set_errhandler @ MPI_File_read_all")
      CALL add_perf(perf_id=28, count=1, msg_size=ndims*msglen)
#else
      MARK_USED(msglen)
      MARK_USED(ndims)
      IF (.NOT. PRESENT(type_descriptor)) &
         CALL cp_abort(__LOCATION__, &
                       "Container for mp_file_descriptor_type must be present in serial call.")
      IF (.NOT. type_descriptor%has_indexing) &
         CALL cp_abort(__LOCATION__, &
                       "File view has not been set in mp_file_descriptor_type.")
      ! Use explicit offsets
      DO i = 1, ndims
         READ (fh%handle, POS=type_descriptor%index_descriptor%chunks(i)) buffer(i)
      END DO
#endif

      CALL mp_timestop(handle)

   END SUBROUTINE mp_file_read_all_chv

! **************************************************************************************************
!> \brief (parallel) Collective, blocking write of a character array to a file. File access pattern
!                    determined by a previously set file view.
!>        (serial)   Unformatted stream write using explicit offsets
!> \param fh     the file handle associated with the output file
!> \param msglen the message length of an individual vector component
!> \param ndims  the number of vector components
!> \param buffer the buffer where the data is placed
!> \param type_descriptor container for the MPI type
!> \author Nico Holmberg [05.2017]
! **************************************************************************************************
   SUBROUTINE mp_file_write_all_chv(fh, msglen, ndims, buffer, type_descriptor)
      CLASS(mp_file_type), INTENT(IN)                       :: fh
      INTEGER, INTENT(IN)                       :: msglen
      INTEGER, INTENT(IN)                       :: ndims
      CHARACTER(LEN=msglen), DIMENSION(ndims), INTENT(IN)   :: buffer
      TYPE(mp_file_descriptor_type), &
         INTENT(IN), OPTIONAL                   :: type_descriptor

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_write_all_chv'

      INTEGER                                   :: handle
#if defined(__parallel)
      INTEGER :: ierr
#else
      INTEGER :: i
#endif

      CALL mp_timeset(routineN, handle)

      IF (fh%handle == mp_file_null_handle) &
         CPABORT("Invalid file handle")

#if defined(__parallel)
      MARK_USED(type_descriptor)
      CALL mpi_file_set_errhandler(fh%handle, MPI_ERRORS_RETURN, ierr)
      CALL MPI_File_write_all(fh%handle, buffer, ndims*msglen, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      IF (ierr .NE. 0) CALL mp_stop(ierr, "mpi_file_set_errhandler @ MPI_File_write_all")
      CALL add_perf(perf_id=28, count=1, msg_size=ndims*msglen)
#else
      MARK_USED(msglen)
      MARK_USED(ndims)
      IF (.NOT. PRESENT(type_descriptor)) &
         CALL cp_abort(__LOCATION__, &
                       "Container for mp_file_descriptor_type must be present in serial call.")
      IF (.NOT. type_descriptor%has_indexing) &
         CALL cp_abort(__LOCATION__, &
                       "File view has not been set in mp_file_descriptor_type.")
      ! Use explicit offsets
      DO i = 1, ndims
         WRITE (fh%handle, POS=type_descriptor%index_descriptor%chunks(i)) buffer(i)
      END DO
#endif

      CALL mp_timestop(handle)

   END SUBROUTINE mp_file_write_all_chv

! **************************************************************************************************
!> \brief Releases the type used for MPI I/O
!> \param type_descriptor the container for the MPI type
!> \author Nico Holmberg [05.2017]
! **************************************************************************************************
   SUBROUTINE mp_file_type_free(type_descriptor)
      TYPE(mp_file_descriptor_type)             :: type_descriptor

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_file_type_free'

      INTEGER                                   :: handle
#if defined(__parallel)
      INTEGER :: ierr
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      CALL MPI_Type_free(type_descriptor%type_handle, ierr)
      IF (ierr /= 0) &
         CPABORT("MPI_Type_free @ "//routineN)
#endif
#if defined(__parallel) && defined(__MPI_F08)
      type_descriptor%type_handle%mpi_val = -1
#else
      type_descriptor%type_handle = -1
#endif
      type_descriptor%length = -1
      IF (type_descriptor%has_indexing) THEN
         NULLIFY (type_descriptor%index_descriptor%index)
         NULLIFY (type_descriptor%index_descriptor%chunks)
         type_descriptor%has_indexing = .FALSE.
      END IF

      CALL mp_timestop(handle)

   END SUBROUTINE mp_file_type_free

! **************************************************************************************************
!> \brief (parallel) Utility routine to determine MPI file access mode based on variables
!                    that in the serial case would get passed to the intrinsic OPEN
!>        (serial)   No action
!> \param mpi_io     flag that determines if MPI I/O will actually be used
!> \param replace    flag that indicates whether file needs to be deleted prior to opening it
!> \param amode      the MPI I/O access mode
!> \param form       formatted or unformatted data?
!> \param action     the variable that determines what to do with file
!> \param status     the status flag:
!> \param position   should the file be appended or rewound
!> \author Nico Holmberg [11.2017]
! **************************************************************************************************
   SUBROUTINE mp_file_get_amode(mpi_io, replace, amode, form, action, status, position)
      LOGICAL, INTENT(INOUT)                             :: mpi_io, replace
      INTEGER, INTENT(OUT)                               :: amode
      CHARACTER(len=*), INTENT(IN)                       :: form, action, status, position

      amode = -1
#if defined(__parallel)
      ! Disable mpi io for unformatted access
      SELECT CASE (form)
      CASE ("FORMATTED")
         ! Do nothing
      CASE ("UNFORMATTED")
         mpi_io = .FALSE.
      CASE DEFAULT
         CPABORT("Unknown MPI file form requested.")
      END SELECT
      ! Determine file access mode (limited set of allowed choices)
      SELECT CASE (action)
      CASE ("WRITE")
         amode = file_amode_wronly
         SELECT CASE (status)
         CASE ("NEW")
            ! Try to open new file for writing, crash if file already exists
            amode = amode + file_amode_create + file_amode_excl
         CASE ("UNKNOWN")
            ! Open file for writing and create it if file does not exist
            amode = amode + file_amode_create
            SELECT CASE (position)
            CASE ("APPEND")
               ! Append existing file
               amode = amode + file_amode_append
            CASE ("REWIND", "ASIS")
               ! Do nothing
            CASE DEFAULT
               CPABORT("Unknown MPI file position requested.")
            END SELECT
         CASE ("OLD")
            SELECT CASE (position)
            CASE ("APPEND")
               ! Append existing file
               amode = amode + file_amode_append
            CASE ("REWIND", "ASIS")
               ! Do nothing
            CASE DEFAULT
               CPABORT("Unknown MPI file position requested.")
            END SELECT
         CASE ("REPLACE")
            ! Overwrite existing file. Must delete existing file first
            amode = amode + file_amode_create
            replace = .TRUE.
         CASE ("SCRATCH")
            ! Disable
            mpi_io = .FALSE.
         CASE DEFAULT
            CPABORT("Unknown MPI file status requested.")
         END SELECT
      CASE ("READ")
         amode = file_amode_rdonly
         SELECT CASE (status)
         CASE ("NEW")
            CPABORT("Cannot read from 'NEW' file.")
         CASE ("REPLACE")
            CPABORT("Illegal status 'REPLACE' for read.")
         CASE ("UNKNOWN", "OLD")
            ! Do nothing
         CASE ("SCRATCH")
            ! Disable
            mpi_io = .FALSE.
         CASE DEFAULT
            CPABORT("Unknown MPI file status requested.")
         END SELECT
      CASE ("READWRITE")
         amode = file_amode_rdwr
         SELECT CASE (status)
         CASE ("NEW")
            ! Try to open new file, crash if file already exists
            amode = amode + file_amode_create + file_amode_excl
         CASE ("UNKNOWN")
            ! Open file and create it if file does not exist
            amode = amode + file_amode_create
            SELECT CASE (position)
            CASE ("APPEND")
               ! Append existing file
               amode = amode + file_amode_append
            CASE ("REWIND", "ASIS")
               ! Do nothing
            CASE DEFAULT
               CPABORT("Unknown MPI file position requested.")
            END SELECT
         CASE ("OLD")
            SELECT CASE (position)
            CASE ("APPEND")
               ! Append existing file
               amode = amode + file_amode_append
            CASE ("REWIND", "ASIS")
               ! Do nothing
            CASE DEFAULT
               CPABORT("Unknown MPI file position requested.")
            END SELECT
         CASE ("REPLACE")
            ! Overwrite existing file. Must delete existing file first
            amode = amode + file_amode_create
            replace = .TRUE.
         CASE ("SCRATCH")
            ! Disable
            mpi_io = .FALSE.
         CASE DEFAULT
            CPABORT("Unknown MPI file status requested.")
         END SELECT
      CASE DEFAULT
         CPABORT("Unknown MPI file action requested.")
      END SELECT
#else
      MARK_USED(replace)
      MARK_USED(form)
      MARK_USED(position)
      MARK_USED(status)
      MARK_USED(action)
      mpi_io = .FALSE.
#endif

   END SUBROUTINE mp_file_get_amode

! **************************************************************************************************
!> \brief Non-blocking send of custom type
!> \param msgin ...
!> \param dest ...
!> \param comm ...
!> \param request ...
!> \param tag ...
! **************************************************************************************************
   SUBROUTINE mp_isend_custom(msgin, dest, comm, request, tag)
      TYPE(mp_type_descriptor_type), INTENT(IN)          :: msgin
      INTEGER, INTENT(IN)                                :: dest
      CLASS(mp_comm_type), INTENT(IN) :: comm
      TYPE(mp_request_type), INTENT(out)                               :: request
      INTEGER, INTENT(in), OPTIONAL                      :: tag

      INTEGER                                            :: ierr, my_tag

      ierr = 0
      my_tag = 0

#if defined(__parallel)
      IF (PRESENT(tag)) my_tag = tag

      CALL mpi_isend(MPI_BOTTOM, 1, msgin%type_handle, dest, my_tag, &
                     comm%handle, request%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_isend @ mp_isend_custom")
#else
      MARK_USED(msgin)
      MARK_USED(dest)
      MARK_USED(comm)
      MARK_USED(tag)
      ierr = 1
      request = mp_request_null
      CALL mp_stop(ierr, "mp_isend called in non parallel case")
#endif
   END SUBROUTINE mp_isend_custom

! **************************************************************************************************
!> \brief Non-blocking receive of vector data
!> \param msgout ...
!> \param source ...
!> \param comm ...
!> \param request ...
!> \param tag ...
! **************************************************************************************************
   SUBROUTINE mp_irecv_custom(msgout, source, comm, request, tag)
      TYPE(mp_type_descriptor_type), INTENT(INOUT)       :: msgout
      INTEGER, INTENT(IN)                                :: source
      CLASS(mp_comm_type), INTENT(IN) :: comm
      TYPE(mp_request_type), INTENT(out)                               :: request
      INTEGER, INTENT(in), OPTIONAL                      :: tag

      INTEGER                                            :: ierr, my_tag

      ierr = 0
      my_tag = 0

#if defined(__parallel)
      IF (PRESENT(tag)) my_tag = tag

      CALL mpi_irecv(MPI_BOTTOM, 1, msgout%type_handle, source, my_tag, &
                     comm%handle, request%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_irecv @ mp_irecv_custom")
#else
      MARK_USED(msgout)
      MARK_USED(source)
      MARK_USED(comm)
      MARK_USED(tag)
      ierr = 1
      request = mp_request_null
      CPABORT("mp_irecv called in non parallel case")
#endif
   END SUBROUTINE mp_irecv_custom

! **************************************************************************************************
!> \brief Window free
!> \param win ...
! **************************************************************************************************
   SUBROUTINE mp_win_free(win)
      CLASS(mp_win_type), INTENT(INOUT)                             :: win

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_win_free'

      INTEGER                                            :: handle
#if defined(__parallel)
      INTEGER :: ierr
#endif

      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      ierr = 0
      CALL mpi_win_free(win%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_win_free @ "//routineN)

      CALL add_perf(perf_id=21, count=1)
#else
      win%handle = mp_win_null_handle
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_win_free

   SUBROUTINE mp_win_assign(win_new, win_old)
      CLASS(mp_win_type), INTENT(OUT) :: win_new
      CLASS(mp_win_type), INTENT(IN) :: win_old

      win_new%handle = win_old%handle

   END SUBROUTINE mp_win_assign

! **************************************************************************************************
!> \brief Window flush
!> \param win ...
! **************************************************************************************************
   SUBROUTINE mp_win_flush_all(win)
      CLASS(mp_win_type), INTENT(IN)                                :: win

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_win_flush_all'

      INTEGER                                            :: handle, ierr

      ierr = 0
      CALL mp_timeset(routineN, handle)

#if defined(__parallel)
      CALL mpi_win_flush_all(win%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_win_flush_all @ "//routineN)
#else
      MARK_USED(win)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_win_flush_all

! **************************************************************************************************
!> \brief Window lock
!> \param win ...
! **************************************************************************************************
   SUBROUTINE mp_win_lock_all(win)
      CLASS(mp_win_type), INTENT(IN)                             :: win

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_win_lock_all'

      INTEGER                                            :: handle, ierr

      ierr = 0
      CALL mp_timeset(routineN, handle)

#if defined(__parallel)

      CALL mpi_win_lock_all(MPI_MODE_NOCHECK, win%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_win_lock_all @ "//routineN)

      CALL add_perf(perf_id=19, count=1)
#else
      MARK_USED(win)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_win_lock_all

! **************************************************************************************************
!> \brief Window lock
!> \param win ...
! **************************************************************************************************
   SUBROUTINE mp_win_unlock_all(win)
      CLASS(mp_win_type), INTENT(IN)                             :: win

      CHARACTER(len=*), PARAMETER :: routineN = 'mp_win_unlock_all'

      INTEGER                                            :: handle, ierr

      ierr = 0
      CALL mp_timeset(routineN, handle)

#if defined(__parallel)

      CALL mpi_win_unlock_all(win%handle, ierr)
      IF (ierr /= 0) CALL mp_stop(ierr, "mpi_win_unlock_all @ "//routineN)

      CALL add_perf(perf_id=19, count=1)
#else
      MARK_USED(win)
#endif
      CALL mp_timestop(handle)
   END SUBROUTINE mp_win_unlock_all

! **************************************************************************************************
!> \brief Tests the MPI library
!> \param comm the relevant, initialized communicator
!> \param npow number of sizes to test, 10**1 .. 10**npow
!> \param output_unit where to direct output
!> \par History
!>      JGH  6-Feb-2001 : Test and performance code
!> \author JGH  1-JAN-2001
!> \note
!>      quickly adapted benchmark code, will only work on an even number of CPUs.
! **************************************************************************************************
   SUBROUTINE mpi_perf_test(comm, npow, output_unit)
      TYPE(mp_comm_type), INTENT(IN) :: comm
      INTEGER, INTENT(IN)                      :: npow, output_unit

#if defined(__parallel)

      INTEGER :: I, ierr, ierror, itask, itests, J, jtask, left, nbufmax, &
                 ncount, Ngrid, Nloc, nprocs, Ntot, partner, right, taskid
      INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcount, rdispl, scount, sdispl
      LOGICAL                                  :: ionode
      REAL(KIND=dp)                            :: maxdiff, res, res2, res3, t1, &
                                                  t2, t3, t4, t5
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: buffer1, buffer2, buffer3, &
                                                  lgrid, lgrid2, lgrid3
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :)                        :: grid, grid2, grid3, &
                                                   send_timings, send_timings2
      REAL(KIND=dp), PARAMETER :: threshold = 1.0E-8_dp

      ! set system sizes !
      ngrid = 10**npow

      CALL mpi_comm_rank(comm%handle, taskid, ierror)
      CALL mpi_comm_size(comm%handle, Nprocs, ierror)
      ionode = (taskid == 0)
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Running with ", nprocs
         WRITE (output_unit, *) "running messages with npow = ", npow
         WRITE (output_unit, *) "use MPI X in the input for larger (e.g. 6) of smaller (e.g. 3) messages"
         IF (MODULO(nprocs, 2) .NE. 0) WRITE (output_unit, *) "Testing only with an even number of tasks"
      END IF

      IF (MODULO(nprocs, 2) .NE. 0) RETURN

      ! equal loads
      Nloc = Ngrid/nprocs
      Ntot = Nprocs*Nloc
      nbufmax = 10**npow
      !
      ALLOCATE (rcount(nprocs))
      ALLOCATE (scount(nprocs))
      ALLOCATE (sdispl(nprocs))
      ALLOCATE (rdispl(nprocs))
      ALLOCATE (buffer1(nbufmax))
      ALLOCATE (buffer2(nbufmax))
      ALLOCATE (buffer3(nbufmax))
      ALLOCATE (grid(Nloc, Nprocs))
      ALLOCATE (grid2(Nloc, Nprocs))
      ALLOCATE (grid3(Nloc, Nprocs))
      ALLOCATE (lgrid(Nloc))
      ALLOCATE (lgrid2(Nloc))
      ALLOCATE (lgrid3(Nloc))
      ALLOCATE (send_timings(0:nprocs - 1, 0:nprocs - 1))
      ALLOCATE (send_timings2(0:nprocs - 1, 0:nprocs - 1))
      buffer1 = 0.0_dp
      buffer2 = 0.0_dp
      buffer3 = 0.0_dp
      ! timings
      send_timings = 0.0_dp
      send_timings2 = 0.0_dp
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ some in memory tests                   ---------------------
      ! -------------------------------------------------------------------------------------------
      CALL MPI_BARRIER(comm%handle, ierror)
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Testing in memory copies just 1 CPU "
         WRITE (output_unit, *) "  could tell something about the motherboard / cache / compiler "
      END IF
      DO i = 1, npow
         ncount = 10**i
         t2 = 0.0E0_dp
         IF (ncount .GT. nbufmax) CPABORT("")
         DO j = 1, 3**(npow - i)
            CALL MPI_BARRIER(comm%handle, ierror)
            t1 = MPI_WTIME()
            buffer2(1:ncount) = buffer1(1:ncount)
            t2 = t2 + MPI_WTIME() - t1 + threshold
         END DO
         CALL MPI_REDUCE(t2, t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm%handle, ierror)
         IF (ionode .AND. output_unit > 0) THEN
            WRITE (output_unit, '(I9,A,F12.4,A)') 8*ncount, " Bytes ", (3**(npow - i))*ncount*8.0E-6_dp/t1, " MB/s"
         END IF
      END DO
      CALL MPI_BARRIER(comm%handle, ierror)
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ some in memory tests                   ---------------------
      ! -------------------------------------------------------------------------------------------
      CALL MPI_BARRIER(comm%handle, ierror)
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Testing in memory copies all cpus"
         WRITE (output_unit, *) "  is the memory bandwidth affected on an SMP machine ?"
      END IF
      DO i = 1, npow
         ncount = 10**i
         t2 = 0.0E0_dp
         IF (ncount .GT. nbufmax) CPABORT("")
         DO j = 1, 3**(npow - i)
            CALL MPI_BARRIER(comm%handle, ierror)
            t1 = MPI_WTIME()
            buffer2(1:ncount) = buffer1(1:ncount)
            t2 = t2 + MPI_WTIME() - t1 + threshold
         END DO
         CALL MPI_REDUCE(t2, t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm%handle, ierror)
         IF (ionode .AND. output_unit > 0) THEN
            WRITE (output_unit, '(I9,A,F12.4,A)') 8*ncount, " Bytes ", (3**(npow - i))*ncount*8.0E-6_dp/t1, " MB/s"
         END IF
      END DO
      CALL MPI_BARRIER(comm%handle, ierror)
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ first test point to point communication ---------------------
      ! -------------------------------------------------------------------------------------------
      CALL MPI_BARRIER(comm%handle, ierror)
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Testing truly point to point communication (i with j only)"
         WRITE (output_unit, *) "  is there some different connection between i j (e.g. shared memory comm)"
      END IF
      ncount = 10**npow
      IF (ionode .AND. output_unit > 0) WRITE (output_unit, *) "For messages of ", ncount*8, " bytes"
      IF (ncount .GT. nbufmax) CPABORT("")
      DO itask = 0, nprocs - 1
         DO jtask = itask + 1, nprocs - 1
            CALL MPI_BARRIER(comm%handle, ierror)
            t1 = MPI_WTIME()
            IF (taskid .EQ. itask) THEN
               CALL MPI_SEND(buffer1, ncount, MPI_DOUBLE_PRECISION, jtask, itask*jtask, comm%handle, ierror)
            END IF
            IF (taskid .EQ. jtask) THEN
               CALL MPI_RECV(buffer1, ncount, MPI_DOUBLE_PRECISION, itask, itask*jtask, comm%handle, MPI_STATUS_IGNORE, ierror)
            END IF
            send_timings(itask, jtask) = MPI_WTIME() - t1 + threshold
         END DO
      END DO
      send_timings2(:, :) = send_timings
      CALL MPI_REDUCE(send_timings2, send_timings, nprocs**2, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm%handle, ierror)
      IF (ionode .AND. output_unit > 0) THEN
         DO itask = 0, nprocs - 1
            DO jtask = itask + 1, nprocs - 1
               WRITE (output_unit, '(I4,I4,F12.4,A)') itask, jtask, ncount*8.0E-6_dp/send_timings(itask, jtask), " MB/s"
            END DO
         END DO
      END IF
      CALL MPI_BARRIER(comm%handle, ierror)
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ second test point to point communication -------------------
      ! -------------------------------------------------------------------------------------------
      CALL MPI_BARRIER(comm%handle, ierror)
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Testing all nearby point to point communication (0,1)(2,3)..."
         WRITE (output_unit, *) "    these could / should all be on the same shared memory node "
      END IF
      DO i = 1, npow
         ncount = 10**i
         t2 = 0.0E0_dp
         IF (ncount .GT. nbufmax) CPABORT("")
         DO j = 1, 3**(npow - i)
            CALL MPI_BARRIER(comm%handle, ierror)
            t1 = MPI_WTIME()
            IF (MODULO(taskid, 2) == 0) THEN
               CALL MPI_SEND(buffer1, ncount, MPI_DOUBLE_PRECISION, taskid + 1, 0, comm%handle, ierror)
            ELSE
               CALL MPI_RECV(buffer1, ncount, MPI_DOUBLE_PRECISION, taskid - 1, 0, comm%handle, MPI_STATUS_IGNORE, ierror)
            END IF
            t2 = t2 + MPI_WTIME() - t1 + threshold
         END DO
         CALL MPI_REDUCE(t2, t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm%handle, ierror)
         IF (ionode .AND. output_unit > 0) THEN
            WRITE (output_unit, '(I9,A,F12.4,A)') 8*ncount, " Bytes ", (3**(npow - i))*ncount*8.0E-6_dp/t1, " MB/s"
         END IF
      END DO
      CALL MPI_BARRIER(comm%handle, ierror)
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ third test point to point communication -------------------
      ! -------------------------------------------------------------------------------------------
      CALL MPI_BARRIER(comm%handle, ierror)
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Testing all far point to point communication (0,nprocs/2),(1,nprocs/2+1),.."
         WRITE (output_unit, *) "    these could all be going over the network, and stress it a lot"
      END IF
      DO i = 1, npow
         ncount = 10**i
         t2 = 0.0E0_dp
         IF (ncount .GT. nbufmax) CPABORT("")
         DO j = 1, 3**(npow - i)
            CALL MPI_BARRIER(comm%handle, ierror)
            t1 = MPI_WTIME()
            ! first half with partner
            IF (taskid .LT. nprocs/2) THEN
               CALL MPI_SEND(buffer1, ncount, MPI_DOUBLE_PRECISION, taskid + nprocs/2, 0, comm%handle, ierror)
            ELSE
               CALL MPI_RECV(buffer1, ncount, MPI_DOUBLE_PRECISION, taskid - nprocs/2, 0, comm%handle, MPI_STATUS_IGNORE, ierror)
            END IF
            t2 = t2 + MPI_WTIME() - t1 + threshold
         END DO
         CALL MPI_REDUCE(t2, t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm%handle, ierror)
         IF (ionode .AND. output_unit > 0) THEN
            WRITE (output_unit, '(I9,A,F12.4,A)') 8*ncount, " Bytes ", (3**(npow - i))*ncount*8.0E-6_dp/t1, " MB/s"
         END IF
      END DO
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ test root to all broadcast               -------------------
      ! -------------------------------------------------------------------------------------------
      CALL MPI_BARRIER(comm%handle, ierror)
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Testing root to all broadcast "
         WRITE (output_unit, *) "    using trees at least ? "
      END IF
      DO i = 1, npow
         ncount = 10**i
         t2 = 0.0E0_dp
         IF (ncount .GT. nbufmax) CPABORT("")
         DO j = 1, 3**(npow - i)
            CALL MPI_BARRIER(comm%handle, ierror)
            t1 = MPI_WTIME()
            CALL MPI_BCAST(buffer1, ncount, MPI_DOUBLE_PRECISION, 0, comm%handle, ierror)
            t2 = t2 + MPI_WTIME() - t1 + threshold
         END DO
         CALL MPI_REDUCE(t2, t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm%handle, ierror)
         IF (ionode .AND. output_unit > 0) THEN
            WRITE (output_unit, '(I9,A,F12.4,A)') 8*ncount, " Bytes ", (3**(npow - i))*ncount*8.0E-6_dp/t1, " MB/s"
         END IF
      END DO
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ test parallel sum like behavior                -------------------
      ! -------------------------------------------------------------------------------------------
      CALL MPI_BARRIER(comm%handle, ierror)
      IF (ionode .AND. output_unit > 0) WRITE (output_unit, *) "Test global summation (mpi_allreduce) "
      DO i = 1, npow
         ncount = 10**i
         t2 = 0.0E0_dp
         IF (ncount .GT. nbufmax) CPABORT("")
         DO j = 1, 3**(npow - i)
            CALL MPI_BARRIER(comm%handle, ierror)
            t1 = MPI_WTIME()
            CALL MPI_ALLREDUCE(buffer1, buffer2, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, comm%handle, ierr)
            t2 = t2 + MPI_WTIME() - t1 + threshold
         END DO
         CALL MPI_REDUCE(t2, t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm%handle, ierror)
         IF (ionode .AND. output_unit > 0) THEN
            WRITE (output_unit, '(I9,A,F12.4,A)') 8*ncount, " Bytes ", (3**(npow - i))*ncount*8.0E-6_dp/t1, " MB/s"
         END IF
      END DO
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ test all to all communication            -------------------
      ! -------------------------------------------------------------------------------------------
      CALL MPI_BARRIER(comm%handle, ierror)
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Test all to all communication (mpi_alltoallv)"
         WRITE (output_unit, *) "    mpi/network getting confused ? "
      END IF
      DO i = 1, npow
         ncount = 10**i
         t2 = 0.0E0_dp
         IF (ncount .GT. nbufmax) CPABORT("")
         scount = ncount/nprocs
         rcount = ncount/nprocs
         DO j = 1, nprocs
            sdispl(j) = (j - 1)*(ncount/nprocs)
            rdispl(j) = (j - 1)*(ncount/nprocs)
         END DO
         DO j = 1, 3**(npow - i)
            CALL MPI_BARRIER(comm%handle, ierror)
            t1 = MPI_WTIME()
            CALL mpi_alltoallv(buffer1, scount, sdispl, MPI_DOUBLE_PRECISION, &
                               buffer2, rcount, rdispl, MPI_DOUBLE_PRECISION, comm%handle, ierr)
            t2 = t2 + MPI_WTIME() - t1 + threshold
         END DO
         CALL MPI_REDUCE(t2, t1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm%handle, ierror)
         IF (ionode .AND. output_unit > 0) THEN
            WRITE (output_unit, '(I9,A,F12.4,A)') 8*(ncount/nprocs)*nprocs, " Bytes ", &
               (3**(npow - i))*(ncount/nprocs)*nprocs*8.0E-6_dp/t1, " MB/s"
         END IF
      END DO

      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ other stuff                            ---------------------
      ! -------------------------------------------------------------------------------------------
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) " Clean tests completed "
         WRITE (output_unit, *) " Testing MPI_REDUCE scatter"
      END IF
      rcount = Nloc
      DO itests = 1, 3
         IF (ionode .AND. output_unit > 0) &
            WRITE (output_unit, *) "------------------------------- test ", itests, " ------------------------"
         ! *** reference ***
         DO j = 1, Nprocs
            DO i = 1, Nloc
               grid(i, j) = MODULO(i*j*taskid, itests)
            END DO
         END DO
         t1 = MPI_WTIME()
         CALL MPI_REDUCE_SCATTER(grid, lgrid, rcount, MPI_DOUBLE_PRECISION, MPI_SUM, comm%handle, ierr)
         t2 = MPI_WTIME() - t1 + threshold
         CALL mpi_allreduce(t2, res, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm%handle, ierr)
         IF (ionode .AND. output_unit > 0) WRITE (output_unit, *) "MPI_REDUCE_SCATTER    ", res
         ! *** simple shift ***
         DO j = 1, Nprocs
            DO i = 1, Nloc
               grid2(i, j) = MODULO(i*j*taskid, itests)
            END DO
         END DO
         left = MODULO(taskid - 1, Nprocs)
         right = MODULO(taskid + 1, Nprocs)
         t3 = MPI_WTIME()
         lgrid2 = 0.0E0_dp
         DO i = 1, Nprocs
            lgrid2(:) = lgrid2 + grid(:, MODULO(taskid - i, Nprocs) + 1)
            IF (i .EQ. nprocs) EXIT
            CALL MPI_SENDRECV_REPLACE(lgrid2, nloc, MPI_DOUBLE_PRECISION, right, 0, left, 0, comm%handle, MPI_STATUS_IGNORE, ierr)
         END DO
         t4 = MPI_WTIME() - t3 + threshold
         CALL mpi_allreduce(t4, res, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm%handle, ierr)
         maxdiff = MAXVAL(ABS(lgrid2 - lgrid))
         CALL mpi_allreduce(maxdiff, res2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm%handle, ierr)
         IF (ionode .AND. output_unit > 0) WRITE (output_unit, *) "MPI_SENDRECV_REPLACE  ", res, res2
         ! *** involved shift ****
         IF (MODULO(nprocs, 2) /= 0) CPABORT("")
         DO j = 1, Nprocs
            DO i = 1, Nloc
               grid3(i, j) = MODULO(i*j*taskid, itests)
            END DO
         END DO
         t3 = MPI_WTIME()
         ! first sum the grid in pairs (0,1),(2,3) should be within an LPAR and fast XXXXXXXXX
         ! 0 will only need parts 0,2,4,... correctly summed
         ! 1 will only need parts 1,3,5,... correctly summed
         ! *** could nicely be generalised ****
         IF (MODULO(taskid, 2) == 0) THEN
            partner = taskid + 1
            DO i = 1, Nprocs, 2 ! sum the full grid with the partner
               CALL MPI_SENDRECV(grid3(1, i + 1), nloc, MPI_DOUBLE_PRECISION, partner, 17, &
                                 lgrid3, nloc, MPI_DOUBLE_PRECISION, partner, 19, comm%handle, MPI_STATUS_IGNORE, ierr)
               grid3(:, i) = grid3(:, i) + lgrid3(:)
            END DO
         ELSE
            partner = taskid - 1
            DO i = 1, Nprocs, 2
               CALL MPI_SENDRECV(grid3(1, i), nloc, MPI_DOUBLE_PRECISION, partner, 19, &
                                 lgrid3, nloc, MPI_DOUBLE_PRECISION, partner, 17, comm%handle, MPI_STATUS_IGNORE, ierr)
               grid3(:, i + 1) = grid3(:, i + 1) + lgrid3(:)
            END DO
         END IF
         t4 = MPI_WTIME() - t3 + threshold
         ! now send a given buffer from 1 to 3 to 5 .. adding the right part of the data
         ! since we've summed an lgrid does only need to pass by even or odd tasks
         left = MODULO(taskid - 2, Nprocs)
         right = MODULO(taskid + 2, Nprocs)
         t3 = MPI_WTIME()
         lgrid3 = 0.0E0_dp
         DO i = 1, Nprocs, 2
            lgrid3(:) = lgrid3 + grid3(:, MODULO(taskid - i - 1, Nprocs) + 1)
            IF (i .EQ. nprocs - 1) EXIT
            CALL MPI_SENDRECV_REPLACE(lgrid3, nloc, MPI_DOUBLE_PRECISION, right, 0, left, 0, comm%handle, MPI_STATUS_IGNORE, ierr)
         END DO
         t5 = MPI_WTIME() - t3 + threshold
         CALL mpi_allreduce(t4, res, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm%handle, ierr)
         CALL mpi_allreduce(t5, res2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm%handle, ierr)
         maxdiff = MAXVAL(ABS(lgrid3 - lgrid))
         CALL mpi_allreduce(maxdiff, res3, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm%handle, ierr)
       IF (ionode .AND. output_unit > 0) WRITE (output_unit, *) "INVOLVED SHIFT        ", res + res2, "(", res, ",", res2, ")", res3
      END DO
      DEALLOCATE (rcount)
      DEALLOCATE (scount)
      DEALLOCATE (sdispl)
      DEALLOCATE (rdispl)
      DEALLOCATE (buffer1)
      DEALLOCATE (buffer2)
      DEALLOCATE (buffer3)
      DEALLOCATE (grid)
      DEALLOCATE (grid2)
      DEALLOCATE (grid3)
      DEALLOCATE (lgrid)
      DEALLOCATE (lgrid2)
      DEALLOCATE (lgrid3)
      DEALLOCATE (send_timings)
      DEALLOCATE (send_timings2)
#else
      MARK_USED(comm)
      MARK_USED(npow)
      IF (output_unit > 0) WRITE (output_unit, *) "No MPI tests for a serial program"
#endif
   END SUBROUTINE mpi_perf_test

! **************************************************************************************************
!> \brief Starts a timer region
!> \param routineN ...
!> \param handle ...
! **************************************************************************************************
   SUBROUTINE mp_timeset(routineN, handle)
      CHARACTER(len=*), INTENT(IN)                       :: routineN
      INTEGER, INTENT(OUT)                               :: handle

      IF (mp_collect_timings) &
         CALL timeset(routineN, handle)
   END SUBROUTINE mp_timeset

! **************************************************************************************************
!> \brief Ends a timer region
!> \param handle ...
! **************************************************************************************************
   SUBROUTINE mp_timestop(handle)
      INTEGER, INTENT(IN)                                :: handle

      IF (mp_collect_timings) &
         CALL timestop(handle)
   END SUBROUTINE mp_timestop

   #:include 'message_passing.fypp'

END MODULE message_passing
