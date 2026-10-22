!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Auxiliary tools to redistribute cp_fm_type and cp_cfm_type matrices before and after
!>        diagonalization. Heuristics are used to determine the optimal number of CPUs for
!>        diagonalization and the input matrices are redistributed if necessary
!> \par History
!>      - [01.2018] moved redistribution related code from cp_fm_syevd here
!>      - [08.2026] added complex matrix (cp_cfm_type) counterparts
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
MODULE cp_fm_diag_utils
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: gcd
   USE message_passing,                 ONLY: mp_para_env_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_diag_utils'

   ! Information on redistribution
   TYPE, PUBLIC :: cp_fm_redistribute_info
      INTEGER :: matrix_order = -1
      INTEGER :: num_pe_old = -1 ! number of processes before a potential redistribute
      INTEGER :: num_pe_new = -1 ! number of processes after a potential redistribute
      INTEGER :: num_pe_opt = -1 ! optimal number of processes for the given matrix
      INTEGER :: num_pe_max_nz_col = -1 ! the maximal number of processes s.t. no column has zero width, may be < 0 if ignored
      LOGICAL :: redistribute = .FALSE. ! whether or not the matrix was actually redistributed
   CONTAINS
      PROCEDURE, PASS(self) :: write => cp_fm_redistribute_info_write
   END TYPE cp_fm_redistribute_info

   ! Container for redistribution settings and temporary work structs
   TYPE cp_fm_redistribute_type
      ! Settings
      INTEGER                                  :: a = -1, x = -1
      LOGICAL                                  :: should_print = .FALSE.
      LOGICAL                                  :: elpa_force_redistribute = .FALSE.
      ! Temporaries
      INTEGER, DIMENSION(:), POINTER           :: group_distribution => NULL(), &
                                                  group_partition => NULL()
      TYPE(cp_blacs_env_type), POINTER         :: blacs_env_new => NULL()
      TYPE(mp_para_env_type), POINTER          :: para_env_new => NULL()
   END TYPE cp_fm_redistribute_type

   ! Permanent instance of the redistribute type
   TYPE(cp_fm_redistribute_type), PRIVATE, &
      SAVE                                     :: work_redistribute

   ! Public subroutines

   PUBLIC :: cp_fm_redistribute_start, &
             cp_fm_redistribute_end, &
             cp_fm_redistribute_init, &
             cp_cfm_redistribute_start, &
             cp_cfm_redistribute_end
#if defined(__parallel)
   PUBLIC :: cp_fm_max_ncpu_non_zero_column
#endif

CONTAINS

! **************************************************************************************************
!> \brief Write the redistribute info nicely formatted to the given I/O unit
!> \param self reference to the cp_fm_redistribute_info instance
!> \param io_unit I/O unit to use for writing
! **************************************************************************************************
   SUBROUTINE cp_fm_redistribute_info_write(self, io_unit)
      CLASS(cp_fm_redistribute_info), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: io_unit

      WRITE (UNIT=io_unit, FMT="(A)") ""
      WRITE (UNIT=io_unit, FMT="(T2,A,T71,I10)") &
         "CP_FM_DIAG| Number of processes over which the matrix is distributed ", self%num_pe_old, &
         "CP_FM_DIAG| Matrix order ", self%matrix_order
      WRITE (UNIT=io_unit, FMT="(T2,A,T71,I10)") &
         "CP_FM_DIAG| Optimal number of CPUs ", self%num_pe_opt
      IF (self%num_pe_max_nz_col < 0) THEN
         WRITE (UNIT=io_unit, FMT="(T2,A,T71,A10)") &
            "CP_FM_DIAG| Maximum number of CPUs (with non-zero columns) ", "<N/A>"
      ELSE
         WRITE (UNIT=io_unit, FMT="(T2,A,T71,I10)") &
            "CP_FM_DIAG| Maximum number of CPUs (with non-zero columns): ", self%num_pe_max_nz_col
      END IF
      IF (self%redistribute) THEN
         WRITE (UNIT=io_unit, FMT="(T2,A,T71,I10)") &
            "CP_FM_DIAG| Number of processes for the redistribution ", self%num_pe_new
      ELSE
         WRITE (UNIT=io_unit, FMT="(T2,A)") &
            "CP_FM_DIAG| The matrix will NOT be redistributed"
      END IF
      WRITE (UNIT=io_unit, FMT="(A)") ""

   END SUBROUTINE cp_fm_redistribute_info_write

! **************************************************************************************************
!> \brief  Releases the temporary storage needed when redistributing arrays
!> \param  has_redistributed flag that determines if the processors holds a part of the
!>                           redistributed array
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
   SUBROUTINE cp_fm_redistribute_work_finalize(has_redistributed)
      LOGICAL, INTENT(IN)                                :: has_redistributed

      IF (ASSOCIATED(work_redistribute%group_distribution)) THEN
         IF (has_redistributed) THEN
            CALL cp_blacs_env_release(work_redistribute%blacs_env_new)
         END IF
         CALL work_redistribute%para_env_new%free()
         DEALLOCATE (work_redistribute%para_env_new)
         DEALLOCATE (work_redistribute%group_distribution)
         DEALLOCATE (work_redistribute%group_partition)
      END IF
      ! Return the temporaries to their initial state; the settings written by
      ! cp_fm_redistribute_init have to survive every redistribution
      NULLIFY (work_redistribute%group_distribution)
      NULLIFY (work_redistribute%group_partition)
      NULLIFY (work_redistribute%blacs_env_new)
      NULLIFY (work_redistribute%para_env_new)

   END SUBROUTINE cp_fm_redistribute_work_finalize

! **************************************************************************************************
!> \brief  Initializes the parameters that determine how to calculate the optimal number of CPUs
!>         for diagonalizing a matrix. The parameters are read from the GLOBAL input section.
!> \param a                integer parameter used to define the rule for determining the optimal
!>                         number of CPUs for diagonalization
!> \param x                integer parameter used to define the rule for determining the optimal
!>                         number of CPUs for diagonalization
!> \param should_print     flag that determines if information about the redistribution process
!>                         should be printed
!> \param elpa_force_redistribute  flag that if redistribution should always be performed when
!>                                 the ELPA diagonalization library is in use
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
   SUBROUTINE cp_fm_redistribute_init(a, x, should_print, elpa_force_redistribute)
      INTEGER, INTENT(IN)                                :: a, x
      LOGICAL, INTENT(IN)                                :: should_print, elpa_force_redistribute

      CPASSERT(a > 0 .AND. x > 0)

      ! Init work first: the structure constructor resets every component, so the
      ! settings have to be written after it, not before
      work_redistribute = cp_fm_redistribute_type()
      work_redistribute%a = a
      work_redistribute%x = x
      work_redistribute%should_print = should_print
      work_redistribute%elpa_force_redistribute = elpa_force_redistribute

   END SUBROUTINE cp_fm_redistribute_init

! **************************************************************************************************
!> \brief  Calculates the optimal number of CPUs for diagonalizing a matrix.
!> \param  size  the size of the diagonalized matrix
!> \return the optimal number of CPUs
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
   PURE FUNCTION cp_fm_diag_get_optimal_ncpu(size) RESULT(ncpu)
      INTEGER, INTENT(IN)                                :: size
      INTEGER                                            :: ncpu

      ncpu = ((size + work_redistribute%a*work_redistribute%x - 1)/ &
              (work_redistribute%a*work_redistribute%x))*work_redistribute%a

   END FUNCTION cp_fm_diag_get_optimal_ncpu

#if defined(__parallel)
! **************************************************************************************************
!> \brief  Finds a safe ELPA layout at or below the requested number of CPUs.
!> \param  matrix_struct the original matrix layout
!> \param  ncpu_max the upper limit on the number of CPUs
!> \param  block_size the block size of the redistributed matrix
!> \return a number of CPUs whose layout has no empty process rows or columns
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
   FUNCTION cp_fm_max_ncpu_non_zero_column(matrix_struct, ncpu_max, block_size) RESULT(ncpu)
      TYPE(cp_fm_struct_type), INTENT(IN)                :: matrix_struct
      INTEGER, INTENT(IN)                                :: ncpu_max, block_size
      INTEGER                                            :: ncpu

      INTEGER                                            :: gcd_max, ipe, jpe, npcol, nprow, nzero
      INTEGER, EXTERNAL                                  :: numroc

      ncpu = MAX(1, MIN(ncpu_max, matrix_struct%para_env%num_pe))
      ! Avoid layouts with odd number of CPUs (blacs grid layout will be square).
      IF (ncpu > 2) ncpu = ncpu - MODULO(ncpu, 2)

      ! Keep a valid original layout only when no redistribution is needed.
      IF (ncpu == matrix_struct%para_env%num_pe .AND. &
          block_size == matrix_struct%nrow_block .AND. block_size == matrix_struct%ncol_block) THEN
         nzero = COUNT(matrix_struct%nrow_locals == 0) + COUNT(matrix_struct%ncol_locals == 0)
         IF (nzero == 0) RETURN
         ! A smaller rank count ensures that an invalid original grid is replaced.
         ncpu = ncpu - nzero
         IF (ncpu > 2) ncpu = ncpu - MODULO(ncpu, 2)
      END IF

      ! Fewer CPUs can mean more process columns (e.g. 32: 4x8, 20: 2x10).
      ! Check the actual candidate grid with the final block size.
      DO WHILE (ncpu > 1)
         ! Determine layout of new blacs grid with ncpu CPUs
         ! (snippet copied from cp_blacs_env.F:cp_blacs_env_create)
         gcd_max = -1
         DO ipe = 1, CEILING(SQRT(REAL(ncpu, dp)))
            jpe = ncpu/ipe
            IF (ipe*jpe /= ncpu) THEN
               CYCLE
            END IF
            IF (gcd(ipe, jpe) >= gcd_max) THEN
               nprow = ipe
               npcol = jpe
               gcd_max = gcd(ipe, jpe)
            END IF
         END DO

         ! Count empty process rows and columns.
         ! (snippet copied from cp_fm_struct.F:cp_fm_struct_create)
         nzero = 0
         DO ipe = 0, npcol - 1
            IF (numroc(matrix_struct%ncol_global, block_size, ipe, 0, npcol) == 0) THEN
               nzero = nzero + 1
            END IF
         END DO
         DO ipe = 0, nprow - 1
            IF (numroc(matrix_struct%nrow_global, block_size, ipe, 0, nprow) == 0) THEN
               nzero = nzero + 1
            END IF
         END DO

         IF (nzero == 0) THEN
            EXIT
         END IF

         ncpu = ncpu - nzero

         IF (ncpu > 2) THEN
            ncpu = ncpu - MODULO(ncpu, 2)
         END IF
      END DO

   END FUNCTION cp_fm_max_ncpu_non_zero_column
#endif

! **************************************************************************************************
!> \brief   Determines the optimal number of CPUs for matrix diagonalization and redistributes
!>          the input matrices if necessary
!> \param matrix           the input cp_fm_type matrix to be diagonalized
!> \param eigenvectors     the cp_fm_type matrix that will hold the eigenvectors of the input matrix
!> \param matrix_new       the redistributed input matrix which will subsequently be diagonalized,
!>                         or a pointer to the original matrix if no redistribution is required
!> \param eigenvectors_new the redistributed eigenvectors matrix, or a pointer to the original
!>                         matrix if no redistribution is required
!> \param caller_is_elpa   flag that determines if ELPA is used for diagonalization
!> \param redist_info      get info about the redistribution
!> \par History
!>      - [01.2018] created by moving redistribution related code from cp_fm_syevd here
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
   SUBROUTINE cp_fm_redistribute_start(matrix, eigenvectors, matrix_new, eigenvectors_new, &
                                       caller_is_elpa, redist_info)

      TYPE(cp_fm_type), INTENT(IN)             :: matrix, eigenvectors
      TYPE(cp_fm_type), INTENT(OUT)            :: matrix_new, eigenvectors_new
      LOGICAL, OPTIONAL, INTENT(IN)            :: caller_is_elpa

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_redistribute_start'

      INTEGER                                  :: handle
      LOGICAL                                  :: is_elpa
      TYPE(cp_fm_redistribute_info), OPTIONAL, INTENT(OUT) :: redist_info

#if defined(__parallel)
      REAL(KIND=dp)                            :: fake_local_data(1, 1)
      INTEGER                                  :: fake_descriptor(9), mepos_old, &
                                                  io_unit, ngroups, ncol_block, blksize, nrow_block
      TYPE(cp_fm_struct_type), POINTER         :: fm_struct_new
      TYPE(mp_para_env_type), POINTER          :: para_env
      TYPE(cp_logger_type), POINTER            :: logger
      TYPE(cp_fm_redistribute_info)          :: rdinfo
#endif

      CALL timeset(routineN, handle)
      is_elpa = .FALSE.
      IF (PRESENT(caller_is_elpa)) THEN
#if defined(__ELPA)
         is_elpa = caller_is_elpa
#else
         CPABORT("CP2K compiled without the ELPA library.")
#endif
      END IF

#if defined(__parallel)

      logger => cp_get_default_logger()
      io_unit = cp_logger_get_default_io_unit(logger)

      ! first figure out the optimal number of cpus
      ! this is pure heuristics, the defaults are based on rosa timings
      ! that demonstrate that timings go up sharply if too many tasks are used
      ! we take a multiple of 4, and approximately n/60
      para_env => matrix%matrix_struct%para_env
      mepos_old = para_env%mepos
      ncol_block = -1 ! normally we also want to adjust the block size according to the optimal # of CPUs
      nrow_block = -1
      blksize = -1

      rdinfo%matrix_order = matrix%matrix_struct%nrow_global
      rdinfo%num_pe_old = para_env%num_pe
      rdinfo%num_pe_opt = cp_fm_diag_get_optimal_ncpu(rdinfo%matrix_order)
      rdinfo%num_pe_new = rdinfo%num_pe_opt
      rdinfo%num_pe_max_nz_col = -1
      rdinfo%redistribute = .FALSE.

      IF (is_elpa) THEN
         CALL cp_fm_get_info(matrix, ncol_block=ncol_block, nrow_block=nrow_block)

         ! On GPUs, ELPA requires the block size to be a power of 2.
         blksize = 1
         DO WHILE (2*blksize <= MIN(nrow_block, ncol_block))
            blksize = blksize*2
         END DO
         nrow_block = blksize
         ncol_block = blksize

         ! Only redistribute when necessary unless the user requests the heuristic.
         rdinfo%num_pe_max_nz_col = cp_fm_max_ncpu_non_zero_column( &
                                    matrix%matrix_struct, rdinfo%num_pe_old, blksize)
         rdinfo%num_pe_new = rdinfo%num_pe_max_nz_col
         IF (work_redistribute%elpa_force_redistribute .AND. rdinfo%num_pe_opt < rdinfo%num_pe_new) THEN
            rdinfo%num_pe_new = cp_fm_max_ncpu_non_zero_column( &
                                matrix%matrix_struct, rdinfo%num_pe_opt, blksize)
         END IF
      END IF

      ! finally, only redistribute if we're going to use less CPUs than before or changed the block size
      rdinfo%redistribute = (rdinfo%num_pe_old > rdinfo%num_pe_new) .OR. (blksize >= 0 .AND. &
                                   ((blksize /= matrix%matrix_struct%ncol_block) .OR. (blksize /= matrix%matrix_struct%nrow_block)))

      IF (work_redistribute%should_print .AND. io_unit > 0) THEN
         IF (is_elpa) THEN
            IF (work_redistribute%elpa_force_redistribute) THEN
               WRITE (UNIT=io_unit, FMT="(T2,A,T78,A3)") &
                  "CP_FM_DIAG| Force redistribute (ELPA):", "YES"
            ELSE
               WRITE (UNIT=io_unit, FMT="(T2,A,T79,A2)") &
                  "CP_FM_DIAG| Force redistribute (ELPA):", "NO"
            END IF
         END IF
         CALL rdinfo%write(io_unit)
      END IF
      CALL para_env%sync()

      ! if the optimal is smaller than num_pe, we will redistribute the input matrix
      IF (rdinfo%redistribute) THEN
         ! split comm, the first num_pe_new tasks will do the work
         ALLOCATE (work_redistribute%group_distribution(0:rdinfo%num_pe_old - 1))
         ALLOCATE (work_redistribute%group_partition(0:1))
         work_redistribute%group_partition = [rdinfo%num_pe_new, rdinfo%num_pe_old - rdinfo%num_pe_new]
         ALLOCATE (work_redistribute%para_env_new)
         CALL work_redistribute%para_env_new%from_split( &
            comm=para_env, ngroups=ngroups, group_distribution=work_redistribute%group_distribution, &
            n_subgroups=2, group_partition=work_redistribute%group_partition)

         IF (work_redistribute%group_distribution(mepos_old) == 0) THEN

            ! create blacs, should inherit the preferences for the layout and so on, from the higher level
            NULLIFY (work_redistribute%blacs_env_new)
            CALL cp_blacs_env_create(blacs_env=work_redistribute%blacs_env_new, para_env=work_redistribute%para_env_new)

            ! create new matrix
            NULLIFY (fm_struct_new)
            IF (nrow_block == -1 .OR. ncol_block == -1) THEN
               CALL cp_fm_struct_create(fmstruct=fm_struct_new, &
                                        para_env=work_redistribute%para_env_new, &
                                        context=work_redistribute%blacs_env_new, &
                                        nrow_global=rdinfo%matrix_order, ncol_global=rdinfo%matrix_order, &
                                        ncol_block=ncol_block, nrow_block=nrow_block)
            ELSE
               CALL cp_fm_struct_create(fmstruct=fm_struct_new, &
                                        para_env=work_redistribute%para_env_new, &
                                        context=work_redistribute%blacs_env_new, &
                                        nrow_global=rdinfo%matrix_order, ncol_global=rdinfo%matrix_order, &
                                        ncol_block=ncol_block, nrow_block=nrow_block, force_block=.TRUE.)
            END IF
            CALL cp_fm_create(matrix_new, matrix_struct=fm_struct_new, name="yevd_new_mat")
            CALL cp_fm_create(eigenvectors_new, matrix_struct=fm_struct_new, name="yevd_new_vec")
            CALL cp_fm_struct_release(fm_struct_new)

            ! redistribute old
            CALL pdgemr2d(rdinfo%matrix_order, rdinfo%matrix_order, matrix%local_data(1, 1), 1, 1, &
                          matrix%matrix_struct%descriptor, &
                          matrix_new%local_data(1, 1), 1, 1, matrix_new%matrix_struct%descriptor, &
                          matrix%matrix_struct%context)
         ELSE
            ! these tasks must help redistribute (they own part of the data),
            ! but need fake 'new' data, and their descriptor must indicate this with -1
            ! see also scalapack comments on pdgemr2d
            fake_descriptor = -1
            CALL pdgemr2d(rdinfo%matrix_order, rdinfo%matrix_order, matrix%local_data(1, 1), 1, 1, &
                          matrix%matrix_struct%descriptor, &
                          fake_local_data(1, 1), 1, 1, fake_descriptor, &
                          matrix%matrix_struct%context)
         END IF
      ELSE
         ! No need to redistribute, just return pointers to the original arrays
         matrix_new = matrix
         eigenvectors_new = eigenvectors
      END IF

      IF (PRESENT(redist_info)) THEN
         redist_info = rdinfo
      END IF
#else

      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(matrix_new)
      MARK_USED(eigenvectors_new)
      MARK_USED(redist_info)
      CPABORT("Routine called in non-parallel case.")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_redistribute_start

! **************************************************************************************************
!> \brief Redistributes eigenvectors and eigenvalues  back to the original communicator group
!> \param matrix           the input cp_fm_type matrix to be diagonalized
!> \param eigenvectors     the cp_fm_type matrix that will hold the eigenvectors of the input matrix
!> \param eig              global array holding the eigenvalues of the input matrixmatrix
!> \param matrix_new       the redistributed input matrix which will subsequently be diagonalized,
!>                         or a pointer to the original matrix if no redistribution is required
!> \param eigenvectors_new the redistributed eigenvectors matrix, or a pointer to the original
!>                         matrix if no redistribution is required
!> \par History
!>      - [01.2018] created by moving redistribution related code from cp_fm_syevd here
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
   SUBROUTINE cp_fm_redistribute_end(matrix, eigenvectors, eig, matrix_new, eigenvectors_new)

      TYPE(cp_fm_type), INTENT(IN)             :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: eig
      TYPE(cp_fm_type), INTENT(INOUT)          :: matrix_new, eigenvectors_new

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_redistribute_end'

      INTEGER                                  :: handle
#if defined(__parallel)
      REAL(KIND=dp)                            :: fake_local_data(1, 1)
      INTEGER                                  :: fake_descriptor(9), mepos_old, n
      TYPE(mp_para_env_type), POINTER          :: para_env
#endif

      CALL timeset(routineN, handle)

#if defined(__parallel)

      ! Check if matrix was redistributed
      IF (ASSOCIATED(work_redistribute%group_distribution)) THEN
         n = matrix%matrix_struct%nrow_global
         para_env => matrix%matrix_struct%para_env
         mepos_old = para_env%mepos

         IF (work_redistribute%group_distribution(mepos_old) == 0) THEN
            ! redistribute results on CPUs that hold the redistributed matrix
            CALL pdgemr2d(n, n, eigenvectors_new%local_data(1, 1), 1, 1, eigenvectors_new%matrix_struct%descriptor, &
                          eigenvectors%local_data(1, 1), 1, 1, eigenvectors%matrix_struct%descriptor, &
                          eigenvectors%matrix_struct%context)
            CALL cp_fm_release(matrix_new)
            CALL cp_fm_release(eigenvectors_new)
         ELSE
            ! these tasks must help redistribute (they own part of the data),
            ! but need fake 'new' data, and their descriptor must indicate this with -1
            ! see also scalapack comments on pdgemr2d
            fake_descriptor = -1
            CALL pdgemr2d(n, n, fake_local_data(1, 1), 1, 1, fake_descriptor, &
                          eigenvectors%local_data(1, 1), 1, 1, eigenvectors%matrix_struct%descriptor, &
                          eigenvectors%matrix_struct%context)
         END IF
         ! free work
         CALL cp_fm_redistribute_work_finalize(work_redistribute%group_distribution(mepos_old) == 0)

         ! finally, also the eigenvalues need to end up on the non-group member tasks
         CALL para_env%bcast(eig, 0)
      END IF

#else

      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(eig)
      MARK_USED(matrix_new)
      MARK_USED(eigenvectors_new)
      CPABORT("Routine called in non-parallel case.")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_redistribute_end

! **************************************************************************************************
!> \brief   Determines the optimal number of CPUs for matrix diagonalization and redistributes
!>          the input matrices if necessary
!> \param matrix           the input cp_cfm_type matrix to be diagonalized
!> \param eigenvectors     the cp_cfm_type matrix that will hold the eigenvectors of the input matrix
!> \param matrix_new       the redistributed input matrix which will subsequently be diagonalized,
!>                         or a pointer to the original matrix if no redistribution is required
!> \param eigenvectors_new the redistributed eigenvectors matrix, or a pointer to the original
!>                         matrix if no redistribution is required
!> \param caller_is_elpa   flag that determines if ELPA is used for diagonalization
!> \param redist_info      get info about the redistribution
!> \par History
!>      - [08.2026] created by mirroring cp_fm_redistribute_start for complex matrices
! **************************************************************************************************
   SUBROUTINE cp_cfm_redistribute_start(matrix, eigenvectors, matrix_new, eigenvectors_new, &
                                        caller_is_elpa, redist_info)

      TYPE(cp_cfm_type), INTENT(IN)             :: matrix, eigenvectors
      TYPE(cp_cfm_type), INTENT(OUT)            :: matrix_new, eigenvectors_new
      LOGICAL, OPTIONAL, INTENT(IN)             :: caller_is_elpa

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_redistribute_start'

      INTEGER                                  :: handle
      LOGICAL                                  :: is_elpa
      TYPE(cp_fm_redistribute_info), OPTIONAL, INTENT(OUT) :: redist_info

#if defined(__parallel)
      COMPLEX(KIND=dp)                         :: fake_local_data(1, 1)
      INTEGER                                  :: fake_descriptor(9), mepos_old, &
                                                  io_unit, ngroups, ncol_block, blksize, nrow_block
      TYPE(cp_fm_struct_type), POINTER         :: fm_struct_new
      TYPE(mp_para_env_type), POINTER          :: para_env
      TYPE(cp_logger_type), POINTER            :: logger
      TYPE(cp_fm_redistribute_info)            :: rdinfo
#endif

      CALL timeset(routineN, handle)
      is_elpa = .FALSE.
      IF (PRESENT(caller_is_elpa)) THEN
#if defined(__ELPA)
         is_elpa = caller_is_elpa
#else
         CPABORT("CP2K compiled without the ELPA library.")
#endif
      END IF

#if defined(__parallel)

      logger => cp_get_default_logger()
      io_unit = cp_logger_get_default_io_unit(logger)

      ! first figure out the optimal number of cpus
      ! this is pure heuristics, the defaults are based on rosa timings
      ! that demonstrate that timings go up sharply if too many tasks are used
      ! we take a multiple of 4, and approximately n/60
      para_env => matrix%matrix_struct%para_env
      mepos_old = para_env%mepos
      ncol_block = -1 ! normally we also want to adjust the block size according to the optimal # of CPUs
      nrow_block = -1
      blksize = -1

      rdinfo%matrix_order = matrix%matrix_struct%nrow_global
      rdinfo%num_pe_old = para_env%num_pe
      rdinfo%num_pe_opt = cp_fm_diag_get_optimal_ncpu(rdinfo%matrix_order)
      rdinfo%num_pe_new = rdinfo%num_pe_opt
      rdinfo%num_pe_max_nz_col = -1
      rdinfo%redistribute = .FALSE.

      IF (is_elpa) THEN
         CALL cp_cfm_get_info(matrix, ncol_block=ncol_block, nrow_block=nrow_block)

         ! On GPUs, ELPA requires the block size to be a power of 2.
         blksize = 1
         DO WHILE (2*blksize <= MIN(nrow_block, ncol_block))
            blksize = blksize*2
         END DO
         nrow_block = blksize
         ncol_block = blksize

         ! Only redistribute when necessary unless the user requests the heuristic.
         rdinfo%num_pe_max_nz_col = cp_fm_max_ncpu_non_zero_column( &
                                    matrix%matrix_struct, rdinfo%num_pe_old, blksize)
         rdinfo%num_pe_new = rdinfo%num_pe_max_nz_col
         IF (work_redistribute%elpa_force_redistribute .AND. rdinfo%num_pe_opt < rdinfo%num_pe_new) THEN
            rdinfo%num_pe_new = cp_fm_max_ncpu_non_zero_column( &
                                matrix%matrix_struct, rdinfo%num_pe_opt, blksize)
         END IF
      END IF

      ! finally, only redistribute if we're going to use less CPUs than before or changed the block size
      rdinfo%redistribute = (rdinfo%num_pe_old > rdinfo%num_pe_new) .OR. (blksize >= 0 .AND. &
                                   ((blksize /= matrix%matrix_struct%ncol_block) .OR. (blksize /= matrix%matrix_struct%nrow_block)))

      IF (work_redistribute%should_print .AND. io_unit > 0) THEN
         IF (is_elpa) THEN
            IF (work_redistribute%elpa_force_redistribute) THEN
               WRITE (UNIT=io_unit, FMT="(T2,A,T78,A3)") &
                  "CP_FM_DIAG| Force redistribute (ELPA):", "YES"
            ELSE
               WRITE (UNIT=io_unit, FMT="(T2,A,T79,A2)") &
                  "CP_FM_DIAG| Force redistribute (ELPA):", "NO"
            END IF
         END IF
         CALL rdinfo%write(io_unit)
      END IF
      CALL para_env%sync()

      ! if the optimal is smaller than num_pe, we will redistribute the input matrix
      IF (rdinfo%redistribute) THEN
         ! split comm, the first num_pe_new tasks will do the work
         ALLOCATE (work_redistribute%group_distribution(0:rdinfo%num_pe_old - 1))
         ALLOCATE (work_redistribute%group_partition(0:1))
         work_redistribute%group_partition = [rdinfo%num_pe_new, rdinfo%num_pe_old - rdinfo%num_pe_new]
         ALLOCATE (work_redistribute%para_env_new)
         CALL work_redistribute%para_env_new%from_split( &
            comm=para_env, ngroups=ngroups, group_distribution=work_redistribute%group_distribution, &
            n_subgroups=2, group_partition=work_redistribute%group_partition)

         IF (work_redistribute%group_distribution(mepos_old) == 0) THEN

            ! create blacs, should inherit the preferences for the layout and so on, from the higher level
            NULLIFY (work_redistribute%blacs_env_new)
            CALL cp_blacs_env_create(blacs_env=work_redistribute%blacs_env_new, para_env=work_redistribute%para_env_new)

            ! create new matrix
            NULLIFY (fm_struct_new)
            IF (nrow_block == -1 .OR. ncol_block == -1) THEN
               CALL cp_fm_struct_create(fmstruct=fm_struct_new, &
                                        para_env=work_redistribute%para_env_new, &
                                        context=work_redistribute%blacs_env_new, &
                                        nrow_global=rdinfo%matrix_order, ncol_global=rdinfo%matrix_order, &
                                        ncol_block=ncol_block, nrow_block=nrow_block)
            ELSE
               CALL cp_fm_struct_create(fmstruct=fm_struct_new, &
                                        para_env=work_redistribute%para_env_new, &
                                        context=work_redistribute%blacs_env_new, &
                                        nrow_global=rdinfo%matrix_order, ncol_global=rdinfo%matrix_order, &
                                        ncol_block=ncol_block, nrow_block=nrow_block, force_block=.TRUE.)
            END IF
            CALL cp_cfm_create(matrix_new, matrix_struct=fm_struct_new, name="zheevd_new_mat")
            CALL cp_cfm_create(eigenvectors_new, matrix_struct=fm_struct_new, name="zheevd_new_vec")
            CALL cp_fm_struct_release(fm_struct_new)

            ! redistribute old
            CALL pzgemr2d(rdinfo%matrix_order, rdinfo%matrix_order, matrix%local_data(1, 1), 1, 1, &
                          matrix%matrix_struct%descriptor, &
                          matrix_new%local_data(1, 1), 1, 1, matrix_new%matrix_struct%descriptor, &
                          matrix%matrix_struct%context)
         ELSE
            ! these tasks must help redistribute (they own part of the data),
            ! but need fake 'new' data, and their descriptor must indicate this with -1
            ! see also scalapack comments on pzgemr2d
            fake_descriptor = -1
            CALL pzgemr2d(rdinfo%matrix_order, rdinfo%matrix_order, matrix%local_data(1, 1), 1, 1, &
                          matrix%matrix_struct%descriptor, &
                          fake_local_data(1, 1), 1, 1, fake_descriptor, &
                          matrix%matrix_struct%context)
         END IF
      ELSE
         ! No need to redistribute, just return pointers to the original arrays
         matrix_new = matrix
         eigenvectors_new = eigenvectors
      END IF

      IF (PRESENT(redist_info)) THEN
         redist_info = rdinfo
      END IF
#else

      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(matrix_new)
      MARK_USED(eigenvectors_new)
      MARK_USED(redist_info)
      CPABORT("Routine called in non-parallel case.")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_redistribute_start

! **************************************************************************************************
!> \brief Redistributes eigenvectors and eigenvalues  back to the original communicator group
!> \param matrix           the input cp_cfm_type matrix to be diagonalized
!> \param eigenvectors     the cp_cfm_type matrix that will hold the eigenvectors of the input matrix
!> \param eig              global array holding the eigenvalues of the input matrixmatrix
!> \param matrix_new       the redistributed input matrix which will subsequently be diagonalized,
!>                         or a pointer to the original matrix if no redistribution is required
!> \param eigenvectors_new the redistributed eigenvectors matrix, or a pointer to the original
!>                         matrix if no redistribution is required
!> \par History
!>      - [08.2026] created by mirroring cp_fm_redistribute_end for complex matrices
! **************************************************************************************************
   SUBROUTINE cp_cfm_redistribute_end(matrix, eigenvectors, eig, matrix_new, eigenvectors_new)

      TYPE(cp_cfm_type), INTENT(IN)             :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: eig
      TYPE(cp_cfm_type), INTENT(INOUT)          :: matrix_new, eigenvectors_new

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_redistribute_end'

      INTEGER                                  :: handle
#if defined(__parallel)
      COMPLEX(KIND=dp)                         :: fake_local_data(1, 1)
      INTEGER                                  :: fake_descriptor(9), mepos_old, n
      TYPE(mp_para_env_type), POINTER          :: para_env
#endif

      CALL timeset(routineN, handle)

#if defined(__parallel)

      ! Check if matrix was redistributed
      IF (ASSOCIATED(work_redistribute%group_distribution)) THEN
         n = matrix%matrix_struct%nrow_global
         para_env => matrix%matrix_struct%para_env
         mepos_old = para_env%mepos

         IF (work_redistribute%group_distribution(mepos_old) == 0) THEN
            ! redistribute results on CPUs that hold the redistributed matrix
            CALL pzgemr2d(n, n, eigenvectors_new%local_data(1, 1), 1, 1, eigenvectors_new%matrix_struct%descriptor, &
                          eigenvectors%local_data(1, 1), 1, 1, eigenvectors%matrix_struct%descriptor, &
                          eigenvectors%matrix_struct%context)
            CALL cp_cfm_release(matrix_new)
            CALL cp_cfm_release(eigenvectors_new)
         ELSE
            ! these tasks must help redistribute (they own part of the data),
            ! but need fake 'new' data, and their descriptor must indicate this with -1
            ! see also scalapack comments on pzgemr2d
            fake_descriptor = -1
            CALL pzgemr2d(n, n, fake_local_data(1, 1), 1, 1, fake_descriptor, &
                          eigenvectors%local_data(1, 1), 1, 1, eigenvectors%matrix_struct%descriptor, &
                          eigenvectors%matrix_struct%context)
         END IF
         ! free work
         CALL cp_fm_redistribute_work_finalize(work_redistribute%group_distribution(mepos_old) == 0)

         ! finally, also the eigenvalues need to end up on the non-group member tasks
         CALL para_env%bcast(eig, 0)
      END IF

#else

      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(eig)
      MARK_USED(matrix_new)
      MARK_USED(eigenvectors_new)
      CPABORT("Routine called in non-parallel case.")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_cfm_redistribute_end

END MODULE cp_fm_diag_utils
