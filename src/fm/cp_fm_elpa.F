!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper for ELPA
!> \author Ole Schuett
! **************************************************************************************************
MODULE cp_fm_elpa
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
   USE cp_fm_struct,                    ONLY: cp_fm_struct_get
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_comm_free,&
                                              mp_comm_split_direct

#include "../base/base_uses.f90"

#if defined (__ELPA)
#include "elpa_kernel_constants.h"
#if (__ELPA < 201611)
   USE elpa2,                           ONLY: solve_evp_real_2stage
#else
   USE elpa2,                           ONLY: elpa_solve_evp_real_2stage_double
#endif
#endif

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_elpa'

   PUBLIC :: cp_fm_diag_elpa
   PUBLIC :: get_elpa_number_kernels
   PUBLIC :: get_elpa_kernel_names
   PUBLIC :: get_elpa_kernel_descriptions
   PUBLIC :: set_elpa_kernel

   INTEGER, SAVE :: elpa_kernel = -1 ! auto

CONTAINS

! **************************************************************************************************
!> \brief Return the number of available ELPA kernels
!> \retval num ...
! **************************************************************************************************
   PURE FUNCTION get_elpa_number_kernels() RESULT(num)
      INTEGER                                            :: num

      num = 19
   END FUNCTION get_elpa_number_kernels

! **************************************************************************************************
!> \brief Returns the names of the availabel ELPA kernels
!> \retval names ...
! **************************************************************************************************
   FUNCTION get_elpa_kernel_names() RESULT(names)
      CHARACTER(len=default_string_length), DIMENSION(19) :: names

      names( 1) = "AUTO"
      names( 2) = "GENERIC"
      names( 3) = "GENERIC_SIMPLE"
      names( 4) = "BGP"
      names( 5) = "BGQ"
      names( 6) = "SSE"
      names( 7) = "SSE_BLOCK2"
      names( 8) = "SSE_BLOCK4"
      names( 9) = "SSE_BLOCK6"
      names(10) = "AVX_BLOCK2"
      names(11) = "AVX_BLOCK4"
      names(12) = "AVX_BLOCK6"
      names(13) = "AVX2_BLOCK2"
      names(14) = "AVX2_BLOCK4"
      names(15) = "AVX2_BLOCK6"
      names(16) = "AVX512_BLOCK2"
      names(17) = "AVX512_BLOCK4"
      names(18) = "AVX512_BLOCK6"
      names(19) = "GPU"

   END FUNCTION get_elpa_kernel_names

! **************************************************************************************************
!> \brief Returns the description texts of the availabel ELPA kernels
!> \retval descriptions ...
! **************************************************************************************************
   FUNCTION get_elpa_kernel_descriptions() RESULT(descriptions)
      CHARACTER(len=default_string_length), DIMENSION(19) :: descriptions

      descriptions( 1) = "Automatically selected kernel"
      descriptions( 2) = "Generic kernel"
      descriptions( 3) = "Simplified generic kernel"
      descriptions( 4) = "Kernel optimized for IBM BGP"
      descriptions( 5) = "Kernel optimized for IBM BGQ"
      descriptions( 6) = "Kernel optimized for x86_64/SSE"
      descriptions( 7) = "Kernel optimized for x86_64/SSE (block=2)"
      descriptions( 8) = "Kernel optimized for x86_64/SSE (block=4)"
      descriptions( 9) = "Kernel optimized for x86_64/SSE (block=6)"
      descriptions(10) = "Kernel optimized for Intel AVX (block=2)"
      descriptions(11) = "Kernel optimized for Intel AVX (block=4)"
      descriptions(12) = "Kernel optimized for Intel AVX (block=6)"
      descriptions(13) = "Kernel optimized for Intel AVX2 (block=2)"
      descriptions(14) = "Kernel optimized for Intel AVX2 (block=4)"
      descriptions(15) = "Kernel optimized for Intel AVX2 (block=6)"
      descriptions(16) = "Kernel optimized for Intel AVX-512 (block=2)"
      descriptions(17) = "Kernel optimized for Intel AVX-512 (block=4)"
      descriptions(18) = "Kernel optimized for Intel AVX-512 (block=6)"
      descriptions(19) = "Kernel targeting GPUs"

   END FUNCTION get_elpa_kernel_descriptions

! **************************************************************************************************
!> \brief Sets the active ELPA kernel.
!> \param kernel Integer between 1 and get_elpa_number_kernels()
! **************************************************************************************************
   SUBROUTINE set_elpa_kernel(kernel)
      INTEGER, INTENT(IN)                                :: kernel

#:def pick_macro(kernel_const)
#if defined(${kernel_const}$)
   elpa_kernel = ${kernel_const}$
#else
   CPABORT("${kernel_const}$ not available")
#endif
#:enddef

      SELECT CASE(kernel)
        CASE(1)
          elpa_kernel = -1 ! auto
        CASE(2)
@:        pick_macro(ELPA2_REAL_KERNEL_GENERIC)
        CASE(3)
@:        pick_macro(ELPA2_REAL_KERNEL_GENERIC_SIMPLE)
        CASE(4)
@:        pick_macro(ELPA2_REAL_KERNEL_BGP)
        CASE(5)
@:        pick_macro(ELPA2_REAL_KERNEL_BGQ)
        CASE(6)
@:        pick_macro(ELPA2_REAL_KERNEL_SSE)
        CASE(7)
@:        pick_macro(ELPA2_REAL_KERNEL_SSE_BLOCK2)
        CASE(8)
@:        pick_macro(ELPA2_REAL_KERNEL_SSE_BLOCK4)
        CASE(9)
@:        pick_macro(ELPA2_REAL_KERNEL_SSE_BLOCK6)
        CASE(10)
@:        pick_macro(ELPA2_REAL_KERNEL_AVX_BLOCK2)
        CASE(11)
@:        pick_macro(ELPA2_REAL_KERNEL_AVX_BLOCK4)
        CASE(12)
@:        pick_macro(ELPA2_REAL_KERNEL_AVX_BLOCK6)
        CASE(13)
@:        pick_macro(ELPA2_REAL_KERNEL_AVX2_BLOCK2)
        CASE(14)
@:        pick_macro(ELPA2_REAL_KERNEL_AVX2_BLOCK4)
        CASE(15)
@:        pick_macro(ELPA2_REAL_KERNEL_AVX2_BLOCK6)
        CASE(16)
@:        pick_macro(ELPA2_REAL_KERNEL_AVX512_BLOCK2)
        CASE(17)
@:        pick_macro(ELPA2_REAL_KERNEL_AVX512_BLOCK4)
        CASE(18)
@:        pick_macro(ELPA2_REAL_KERNEL_AVX512_BLOCK6)
        CASE(19)
@:        pick_macro(ELPA2_REAL_KERNEL_GPU)
        CASE DEFAULT
          CPABORT("Invalid ELPA kernel selected")
     END SELECT

   END SUBROUTINE set_elpa_kernel

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
! **************************************************************************************************
   SUBROUTINE cp_fm_diag_elpa(matrix, eigenvectors, eigenvalues)

      TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:)              :: eigenvalues

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_diag_elpa'

      INTEGER                                  :: handle
#if defined(__ELPA)
      INTEGER                                  :: comm_col, comm_row, group, &
                                                  mypcol, myprow, n, &
                                                  n_rows, n_cols, &
                                                  nblk, neig
      LOGICAL                                  :: success
      REAL(KIND=dp), DIMENSION(:), POINTER     :: eval
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: m, v
      TYPE(cp_blacs_env_type), POINTER         :: context
#endif

      CALL timeset(routineN, handle)
#if defined(__ELPA)

      n = matrix%matrix_struct%nrow_global
      m => matrix%local_data
      context => matrix%matrix_struct%context
      myprow = context%mepos(1)
      mypcol = context%mepos(2)
      group = matrix%matrix_struct%para_env%group

      ! For ELPA, the MPI communicators along rows/cols are sufficient
      ! mpi communicators are created, needed for communicating within
      ! rows or columns of processes
      CALL mp_comm_split_direct(group, comm_row, mypcol, myprow)
      CALL mp_comm_split_direct(group, comm_col, myprow, mypcol)

      v => eigenvectors%local_data

      ! elpa needs the full matrix
      CALL cp_fm_upper_to_full(matrix, eigenvectors)

      CALL cp_fm_struct_get(matrix%matrix_struct, &
                            local_leading_dimension=n_rows, &
                            ncol_local=n_cols, &
                            nrow_block=nblk)

      neig = SIZE(eigenvalues, 1)

      ! ELPA will fail in 'solve_tridi', with no useful error message
      IF (n_cols == 0) &
         CPABORT("ELPA [pre-fail]: Problem contains processor column with zero width.")

      ! the full eigenvalues vector is needed
      ALLOCATE (eval(n))


      ! Make actual call to ELPA to calculate eigenvalues/eigenvectors

!---------------------------------------------------------------------------------------------------
#if (__ELPA==201112) || (__ELPA==201308) || (__ELPA==201311)
      IF (elpa_kernel /= -1) CPABORT("The only available ELPA kernel is AUTO.")
      CALL solve_evp_real_2stage(n, neig, m, n_rows, eval, v, n_rows, nblk, comm_row, comm_col, group)
      success = .TRUE.

!---------------------------------------------------------------------------------------------------
#elif (__ELPA==201406) || (__ELPA==201502)
      IF (elpa_kernel == -1) THEN ! auto
         success = solve_evp_real_2stage(n, neig, m, n_rows, eval, v, n_rows, nblk, &
                                         comm_row, comm_col, group, &
                                         useQR=.FALSE.)
      ELSE
         success = solve_evp_real_2stage(n, neig, m, n_rows, eval, v, n_rows, nblk, &
                                         comm_row, comm_col, group, &
                                         THIS_REAL_ELPA_KERNEL_API=kernel_type, &
                                         useQR=.FALSE.)
      END IF

!---------------------------------------------------------------------------------------------------
#elif (__ELPA==201511) || (__ELPA==201605)
      IF (elpa_kernel == -1) THEN ! auto
         success = solve_evp_real_2stage(n, neig, m, n_rows, eval, v, n_rows, nblk, n_cols, &
                                         comm_row, comm_col, group, &
                                         useQR=.FALSE.)
      ELSE
         success = solve_evp_real_2stage(n, neig, m, n_rows, eval, v, n_rows, nblk, n_cols, &
                                         comm_row, comm_col, group, &
                                         THIS_REAL_ELPA_KERNEL_API=elpa_kernel, &
                                         useQR=.FALSE.)
      END IF

!---------------------------------------------------------------------------------------------------
#elif (__ELPA==201611)
      IF (elpa_kernel == -1) THEN ! auto
         success = elpa_solve_evp_real_2stage_double(n, neig, m, n_rows, eval, v, n_rows, nblk, n_cols, &
                                                     comm_row, comm_col, group, &
                                                     useQR=.FALSE.)
      ELSE
         success = elpa_solve_evp_real_2stage_double(n, neig, m, n_rows, eval, v, n_rows, nblk, n_cols, &
                                                     comm_row, comm_col, group, &
                                                     THIS_REAL_ELPA_KERNEL_API=kernel_type, &
                                                     useQR=.FALSE.)
      END IF

!---------------------------------------------------------------------------------------------------
#else
   Error: Unknown ELPA version, please specifiy library release date via __ELPA=YYYYMM
#endif
!---------------------------------------------------------------------------------------------------

      IF (.NOT. success) &
         CPABORT("ELPA failed to diagonalize a matrix")

      eigenvalues(1:neig) = eval(1:neig)

      DEALLOCATE (eval)

      ! mpi communicators are freed
      CALL mp_comm_free(comm_row)
      CALL mp_comm_free(comm_col)

#else

      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(eigenvalues)

      CPABORT("CP2K compiled without the ELPA library.")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_diag_elpa

END MODULE cp_fm_elpa
