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

#if defined (__ELPA) || defined(__ELPA2) || defined(__ELPA3)
!#if defined (__ELPA) || defined(__ELPA2)
  USE ELPA2,                            ONLY: solve_evp_real_2stage ! ELPA[1-2] and deprecated ELPA3 interface
!#elif defined (__ELPA3)
!  USE ELPA2,                           ONLY: elpa_solve_evp_real_2stage ! new ELPA3 interface
#endif

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_elpa'

   PUBLIC :: cp_fm_diag_elpa
   PUBLIC :: get_elpa_number_kernels
   PUBLIC :: get_elpa_kernel_names
   PUBLIC :: get_elpa_kernel_descriptions
   PUBLIC :: set_elpa_kernel

   INTEGER, SAVE :: elpa_kernel = -1

CONTAINS

! **************************************************************************************************
!> \brief Sets the active ELPA kernel.
!> \param kernel Integer between 1 and get_elpa_number_kernels()
! **************************************************************************************************
   SUBROUTINE set_elpa_kernel(kernel)
      INTEGER, INTENT(IN)                                :: kernel

      CPASSERT(1 <= kernel .AND. kernel <= get_elpa_number_kernels())
      elpa_kernel = kernel

   END SUBROUTINE set_elpa_kernel

! **************************************************************************************************
!> \brief Return the number of available ELPA kernels
!> \retval num ...
! **************************************************************************************************
   PURE FUNCTION get_elpa_number_kernels() RESULT(num)
      INTEGER                                            :: num

      num = 8
   END FUNCTION get_elpa_number_kernels

! **************************************************************************************************
!> \brief Returns the names of the availabel ELPA kernels
!> \retval names ...
! **************************************************************************************************
   FUNCTION get_elpa_kernel_names() RESULT(names)
      CHARACTER(len=default_string_length), DIMENSION(8) :: names

      names(1) = "GENERIC"
      names(2) = "GENERIC_SIMPLE"
      names(3) = "BGP"
      names(4) = "BGQ"
      names(5) = "SSE"
      names(6) = "AVX_BLOCK2"
      names(7) = "AVX_BLOCK4"
      names(8) = "AVX_BLOCK6"

   END FUNCTION get_elpa_kernel_names

! **************************************************************************************************
!> \brief Returns the description texts of the availabel ELPA kernels
!> \retval descriptions ...
! **************************************************************************************************
   FUNCTION get_elpa_kernel_descriptions() RESULT(descriptions)
      CHARACTER(len=default_string_length), DIMENSION(8) :: descriptions

      descriptions(1) = "Generic kernel"
      descriptions(2) = "Simplified generic kernel"
      descriptions(3) = "Kernel optimized for IBM BGP"
      descriptions(4) = "Kernel optimized for IBM BGQ"
      descriptions(5) = "Kernel wiwth assembler for SSE vectorization"
      descriptions(6) = "Kernel optimized for x86_64 using SSE2/SSE3 (Intel)"
      descriptions(7) = "Kernel optimized for x86_64 using SSE2/SSE4 (AMD Bulldozer)"
      descriptions(8) = "Kernel optimized for x86_64 using SSE2/SSE4 (block6)"

   END FUNCTION get_elpa_kernel_descriptions

! **************************************************************************************************
!> \brief Performs a diagonalization using ELPA
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
! **************************************************************************************************
   SUBROUTINE cp_fm_diag_elpa(matrix, eigenvectors, eigenvalues)

      TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:)              :: eigenvalues

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_diag_elpa'

      INTEGER                                  :: handle
#if defined(__SCALAPACK) && (defined(__ELPA) || defined (__ELPA2) || defined(__ELPA3))
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

#if defined(__SCALAPACK) && (defined(__ELPA) || defined (__ELPA2) || defined (__ELPA3))

      CPASSERT(1 <= elpa_kernel .AND. elpa_kernel <= get_elpa_number_kernels())

      success = .TRUE.

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

      ! Calculate eigenvalues/eigenvectors

#if defined (__ELPA)
      CALL solve_evp_real_2stage(n, neig, m, n_rows, eval, v, n_rows, nblk, comm_row, comm_col, group)
#elif defined (__ELPA2)
      success = solve_evp_real_2stage(n, neig, m, n_rows, eval, v, n_rows, nblk, comm_row, comm_col, group, &
                                      THIS_REAL_ELPA_KERNEL_API=elpa_kernel, useQR=.FALSE.)
#else
! deprecated ELPA3 interface
      success = solve_evp_real_2stage(n, neig, m, n_rows, eval, v, n_rows, nblk, n_cols, comm_row, comm_col, group, &
                                      THIS_REAL_ELPA_KERNEL_API=elpa_kernel, useQR=.FALSE.)
! new ELPA3 interface
!      success = elpa_solve_evp_real_2stage(n, neig, m, n_rows, eval, v, n_rows, nblk, n_cols, &
!                                           comm_row, comm_col, group, &
!                                           THIS_REAL_ELPA_KERNEL_API=elpa_kernel, useQR=.FALSE.)
#endif

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
