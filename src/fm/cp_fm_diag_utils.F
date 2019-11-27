!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Auxiliary tools to redistribute cp_fm_type matrices before and after diagonalization.
!>        Heuristics are used to determine the optimal number of CPUs for diagonalization and the
!>        input matrices are redistributed if necessary
!> \par History
!>      - [01.2018] moved redistribution related code from cp_fm_syevd here
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
MODULE cp_fm_diag_utils
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
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
   USE cp_para_env,                     ONLY: cp_para_env_create,&
                                              cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: gcd
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_comm_free,&
                                              mp_comm_split,&
                                              mp_sync
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_diag_utils'

   ! Information on redistribution
   TYPE, PUBLIC :: cp_fm_redistribute_info
      INTEGER :: matrix_order
      INTEGER :: num_pe_old ! number of processes before a potential redistribute
      INTEGER :: num_pe_new ! number of processes after a potential redistribute
      INTEGER :: num_pe_opt ! optimal number of processes for the given matrix
      INTEGER :: num_pe_max_nz_col ! the maximal number of processes s.t. no column has zero width, may be < 0 if ignored
      LOGICAL :: redistribute ! whether or not the matrix was actually redistributed
   CONTAINS
      PROCEDURE, PASS(self) :: write => cp_fm_redistribute_info_write
   END TYPE

   ! Container for redistribution settings and temporary work structs
   TYPE cp_fm_redistribute_type
      ! Settings
      INTEGER                                  :: a, x
      LOGICAL                                  :: should_print
      LOGICAL                                  :: elpa_force_redistribute
      ! Temporaries
      INTEGER                                  :: subgroup
      INTEGER, DIMENSION(:), POINTER           :: group_distribution, &
                                                  group_partition
      TYPE(cp_blacs_env_type), POINTER         :: blacs_env_new
      TYPE(cp_para_env_type), POINTER          :: para_env_new
   END TYPE cp_fm_redistribute_type

   ! Permanent instance of the redistribute type
   TYPE(cp_fm_redistribute_type), PRIVATE, &
      SAVE                                     :: work_redistribute

   ! Public subroutines

   PUBLIC :: cp_fm_redistribute_start, &
             cp_fm_redistribute_end, &
             cp_fm_redistribute_init

CONTAINS

! **************************************************************************************************
!> \brief Write the redistribute info nicely formatted to the given I/O unit
!> \param self reference to the cp_fm_redistribute_info instance
!> \param io_unit I/O unit to use for writing
! **************************************************************************************************
   SUBROUTINE cp_fm_redistribute_info_write(self, io_unit)
      CLASS(cp_fm_redistribute_info), INTENT(IN) :: self
      INTEGER, INTENT(IN) :: io_unit

      WRITE (io_unit, '(/,A,I10)') "CP_FM_DIAG| Matrix order:           ", self%matrix_order
      WRITE (io_unit, '(A,I6,A)') "CP_FM_DIAG| Matrix distributed on ", self%num_pe_old, " processes"
      WRITE (io_unit, '(A,I5)') "CP_FM_DIAG| Optimal number of CPUs:      ", self%num_pe_opt
      IF (self%num_pe_max_nz_col < 0) THEN
         WRITE (io_unit, '(A,A)') "CP_FM_DIAG| Max number of CPUs (with non-zero columns): ", "<N/A>"
      ELSE
         WRITE (io_unit, '(A,I5)') "CP_FM_DIAG| Max number of CPUs (with non-zero columns): ", self%num_pe_max_nz_col
      END IF
      IF (self%redistribute) THEN
         WRITE (io_unit, '(A,I5,A)') "CP_FM_DIAG| The matrix will be redistributed onto ", self%num_pe_new, " processes"
      ELSE
         WRITE (io_unit, '(A)') "CP_FM_DIAG| The matrix will NOT be redistributed"
      END IF
      WRITE (io_unit, '(A)') " "
   END SUBROUTINE

! **************************************************************************************************
!> \brief  Initializes temporary storage needed when redistributing arrays
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
   SUBROUTINE cp_fm_redistribute_work_init()

      work_redistribute%subgroup = -1
      NULLIFY (work_redistribute%group_distribution)
      NULLIFY (work_redistribute%group_partition)
      NULLIFY (work_redistribute%blacs_env_new)
      NULLIFY (work_redistribute%para_env_new)

   END SUBROUTINE cp_fm_redistribute_work_init

! **************************************************************************************************
!> \brief  Releases the temporary storage needed when redistributing arrays
!> \param  has_redistributed flag that determines if the processors holds a part of the
!>                           redistributed array
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
   SUBROUTINE cp_fm_redistribute_work_finalize(has_redistributed)
      LOGICAL                                            :: has_redistributed

      IF (ASSOCIATED(work_redistribute%group_distribution)) THEN
         IF (has_redistributed) THEN
            CALL cp_blacs_env_release(work_redistribute%blacs_env_new)
            CALL cp_para_env_release(work_redistribute%para_env_new)
         ELSE
            CALL mp_comm_free(work_redistribute%subgroup)
         ENDIF
         DEALLOCATE (work_redistribute%group_distribution)
         DEALLOCATE (work_redistribute%group_partition)
      END IF
      ! Return work to its initial state
      CALL cp_fm_redistribute_work_init()

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
!> \param elpa_force_redistribute  flag that if redistribution should always be perfomed when
!>                                 the ELPA diagonalization library is in use
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
   SUBROUTINE cp_fm_redistribute_init(a, x, should_print, elpa_force_redistribute)
      INTEGER                                            :: a, x
      LOGICAL                                            :: should_print, elpa_force_redistribute

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_redistribute_init', &
         routineP = moduleN//':'//routineN

      work_redistribute%a = a
      work_redistribute%x = x
      work_redistribute%should_print = should_print
      work_redistribute%elpa_force_redistribute = elpa_force_redistribute
      ! Init work
      CALL cp_fm_redistribute_work_init()

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

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_diag_get_optimal_ncpu', &
         routineP = moduleN//':'//routineN

      ncpu = ((size + work_redistribute%a*work_redistribute%x - 1)/ &
              (work_redistribute%a*work_redistribute%x))*work_redistribute%a

   END FUNCTION cp_fm_diag_get_optimal_ncpu

! **************************************************************************************************
!> \brief  Determines the largest number of CPUs a matrix can be distributed on without any of the
!>         processors getting a zero-width column (currently only needed for ELPA).
!> \param  matrix the matrix that will be diagonalized
!> \return the maximum number of CPUs for ELPA
!> \author Nico Holmberg [01.2018]
! **************************************************************************************************
   FUNCTION cp_fm_max_ncpu_non_zero_column(matrix) RESULT(ncpu)
      TYPE(cp_fm_type), POINTER                :: matrix
      INTEGER                                  :: ncpu

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_max_ncpu_non_zero_column', &
                                     routineP = moduleN//':'//routineN

#if defined(__SCALAPACK)
      INTEGER                                  :: num_pe_old, nzero, &
                                                  ncol_global, ncol_block, nrow_block, &
                                                  nrow_global, gcd_max, ipe, jpe, npe, &
                                                  nprow, npcol
      INTEGER, DIMENSION(:), POINTER           :: ncol_locals
      INTEGER, EXTERNAL                        :: numroc
      LOGICAL                                  :: max_cpu_found
#endif

      ncpu = 0
#if defined(__SCALAPACK)
      NULLIFY (ncol_locals)
      ! First check if there are any zero width columns in current layout
      CALL cp_fm_get_info(matrix, ncol_locals=ncol_locals, &
                          nrow_global=nrow_global, ncol_global=ncol_global, &
                          nrow_block=nrow_block, ncol_block=ncol_block)
      nzero = COUNT(ncol_locals == 0)
      num_pe_old = matrix%matrix_struct%para_env%num_pe
      ncpu = num_pe_old - nzero
      ! Avoid layouts with odd number of CPUs (blacs grid layout will be square)
      IF (ncpu /= 1) ncpu = ncpu - MODULO(ncpu, 2)
      ! Iteratively search for the maximum number of CPUs for ELPA
      ! On each step, we test whether the blacs grid created with ncpu processes
      ! contains any columns with zero width
      max_cpu_found = .FALSE.
      IF (ncpu == num_pe_old) max_cpu_found = .TRUE.
      DO WHILE (.NOT. max_cpu_found)
         ! Determine layout of new blacs grid with ncpu CPUs
         ! (snippet copied from cp_blacs_env.F:cp_blacs_env_create)
         npe = ncpu
         gcd_max = -1
         DO ipe = 1, CEILING(SQRT(REAL(npe, dp)))
            jpe = npe/ipe
            IF (ipe*jpe .NE. npe) CYCLE
            IF (gcd(ipe, jpe) >= gcd_max) THEN
               nprow = ipe
               npcol = jpe
               gcd_max = gcd(ipe, jpe)
            ENDIF
         END DO
         ! Test if there are any columns with zero width
         ! (snippet copied from cp_fm_struct.F:cp_fm_struct_create)
         nzero = 0
         DO ipe = 0, npcol - 1
            IF (numroc(ncol_global, ncol_block, ipe, 0, npcol) == 0) &
               nzero = nzero + 1
         END DO
         IF (nzero == 0) THEN
            max_cpu_found = .TRUE.
         ELSE
            ncpu = ncpu - nzero
            IF (ncpu /= 1) ncpu = ncpu - MODULO(ncpu, 2)
         END IF
      END DO
      CPASSERT(ncpu .GT. 0)
#else

      MARK_USED(matrix)
      CPABORT("Routine called in non-parallel case.")
#endif

   END FUNCTION cp_fm_max_ncpu_non_zero_column

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

      TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
      TYPE(cp_fm_type), POINTER, INTENT(OUT)   :: matrix_new, eigenvectors_new
      LOGICAL, OPTIONAL, INTENT(IN)            :: caller_is_elpa

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_redistribute_start', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: handle
      LOGICAL                                  :: is_elpa
      TYPE(cp_fm_redistribute_info), OPTIONAL, INTENT(OUT) :: redist_info

#if defined(__SCALAPACK)
      REAL(KIND=dp)                            :: fake_local_data(1, 1)
      INTEGER                                  :: fake_descriptor(9), mepos_old, &
                                                  io_unit, ngroups
      TYPE(cp_fm_struct_type), POINTER         :: fm_struct_new
      TYPE(cp_para_env_type), POINTER          :: para_env
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

#if defined(__SCALAPACK)

      NULLIFY (matrix_new)
      NULLIFY (eigenvectors_new)
      logger => cp_get_default_logger()
      io_unit = cp_logger_get_default_io_unit(logger)

      ! first figure out the optimal number of cpus
      ! this is pure heuristics, the defaults are based on rosa timings
      ! that demonstrate that timings go up sharply if too many tasks are used
      ! we take a multiple of 4, and approximately n/60
      para_env => matrix%matrix_struct%para_env
      mepos_old = para_env%mepos

      rdinfo%matrix_order = matrix%matrix_struct%nrow_global
      rdinfo%num_pe_old = para_env%num_pe
      rdinfo%num_pe_opt = cp_fm_diag_get_optimal_ncpu(rdinfo%matrix_order)
      rdinfo%num_pe_new = rdinfo%num_pe_opt
      rdinfo%num_pe_max_nz_col = -1
      rdinfo%redistribute = .FALSE.

      IF (.NOT. is_elpa) THEN
         rdinfo%redistribute = (rdinfo%num_pe_new < rdinfo%num_pe_old)
      ELSE
         ! Diagonalization with ELPA fails when a processor column has zero width
         ! Determine the maximum number of CPUs the matrix can be distributed without zero-width columns
         rdinfo%num_pe_max_nz_col = cp_fm_max_ncpu_non_zero_column(matrix)
         IF (rdinfo%num_pe_old > rdinfo%num_pe_max_nz_col) THEN
            ! Must redistribute to avoid crash if we exceed the max number of processors
            rdinfo%num_pe_new = rdinfo%num_pe_max_nz_col
            rdinfo%redistribute = .TRUE.
         END IF

         IF (work_redistribute%elpa_force_redistribute .AND. rdinfo%num_pe_opt < rdinfo%num_pe_max_nz_col) THEN
            ! Use heuristics to determine the need for redistribution (when num_pe_opt is smaller than the safe maximum)
            rdinfo%num_pe_new = rdinfo%num_pe_opt
            rdinfo%redistribute = (rdinfo%num_pe_old > rdinfo%num_pe_new)
         END IF
      END IF

      IF (work_redistribute%should_print .AND. io_unit > 0) &
         CALL rdinfo%write(io_unit)

      ! if the optimal is smaller than num_pe, we will redistribute the input matrix
      IF (rdinfo%redistribute) THEN
         ! split comm, the first num_pe_new tasks will do the work
         ALLOCATE (work_redistribute%group_distribution(0:rdinfo%num_pe_old - 1))
         ALLOCATE (work_redistribute%group_partition(0:1))
         work_redistribute%group_partition = (/rdinfo%num_pe_new, rdinfo%num_pe_old - rdinfo%num_pe_new/)
         CALL mp_comm_split(comm=para_env%group, sub_comm=work_redistribute%subgroup, &
                            ngroups=ngroups, group_distribution=work_redistribute%group_distribution, &
                            n_subgroups=2, group_partition=work_redistribute%group_partition)

         IF (work_redistribute%group_distribution(mepos_old) == 0) THEN

            ! create para_env, might need a proper bound to this para_env
            NULLIFY (work_redistribute%para_env_new)
            CALL cp_para_env_create(work_redistribute%para_env_new, work_redistribute%subgroup)
            ! test a sync
            CALL mp_sync(work_redistribute%para_env_new%group)

            ! create blacs, should inherit the preferences for the layout and so on, from the higher level
            NULLIFY (work_redistribute%blacs_env_new)
            CALL cp_blacs_env_create(blacs_env=work_redistribute%blacs_env_new, para_env=work_redistribute%para_env_new)

            ! create new matrix
            NULLIFY (fm_struct_new)
            CALL cp_fm_struct_create(fmstruct=fm_struct_new, &
                                     para_env=work_redistribute%para_env_new, &
                                     context=work_redistribute%blacs_env_new, &
                                     nrow_global=rdinfo%matrix_order, ncol_global=rdinfo%matrix_order)
            CALL cp_fm_create(matrix_new, matrix_struct=fm_struct_new, name="yevd_new_mat")
            CALL cp_fm_create(eigenvectors_new, matrix_struct=fm_struct_new, name="yevd_new_vec")
            CALL cp_fm_struct_release(fm_struct_new)

            ! redistribute old
            CALL pdgemr2d(rdinfo%matrix_order, rdinfo%matrix_order, matrix%local_data(1, 1), 1, 1, &
                          matrix%matrix_struct%descriptor, &
                          matrix_new%local_data(1, 1), 1, 1, matrix_new%matrix_struct%descriptor, &
                          matrix%matrix_struct%context%group)
         ELSE
            ! these tasks must help redistribute (they own part of the data),
            ! but need fake 'new' data, and their descriptor must indicate this with -1
            ! see also scalapack comments on pdgemr2d
            fake_descriptor = -1
            CALL pdgemr2d(rdinfo%matrix_order, rdinfo%matrix_order, matrix%local_data(1, 1), 1, 1, &
                          matrix%matrix_struct%descriptor, &
                          fake_local_data(1, 1), 1, 1, fake_descriptor, &
                          matrix%matrix_struct%context%group)
         ENDIF
      ELSE
         ! No need to redistribute, just return pointers to the original arrays
         matrix_new => matrix
         eigenvectors_new => eigenvectors
      ENDIF

      IF (PRESENT(redist_info)) &
         redist_info = rdinfo
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

      TYPE(cp_fm_type), POINTER                :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:)              :: eig
      TYPE(cp_fm_type), POINTER, INTENT(OUT)   :: matrix_new, eigenvectors_new

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_redistribute_end', &
                                     routineP = moduleN//':'//routineN

      INTEGER                                  :: handle
#if defined(__SCALAPACK)
      REAL(KIND=dp)                            :: fake_local_data(1, 1)
      INTEGER                                  :: fake_descriptor(9), mepos_old, n
      TYPE(cp_para_env_type), POINTER          :: para_env
#endif

      CALL timeset(routineN, handle)

#if defined(__SCALAPACK)

      ! Check if matrix was redistributed
      IF (ASSOCIATED(work_redistribute%group_distribution)) THEN
         n = matrix%matrix_struct%nrow_global
         para_env => matrix%matrix_struct%para_env
         mepos_old = para_env%mepos

         IF (work_redistribute%group_distribution(mepos_old) == 0) THEN
            ! redistribute results on CPUs that hold the redistributed matrix
            CALL pdgemr2d(n, n, eigenvectors_new%local_data(1, 1), 1, 1, eigenvectors_new%matrix_struct%descriptor, &
                          eigenvectors%local_data(1, 1), 1, 1, eigenvectors%matrix_struct%descriptor, &
                          eigenvectors%matrix_struct%context%group)
            CALL cp_fm_release(matrix_new)
            CALL cp_fm_release(eigenvectors_new)
         ELSE
            ! these tasks must help redistribute (they own part of the data),
            ! but need fake 'new' data, and their descriptor must indicate this with -1
            ! see also scalapack comments on pdgemr2d
            fake_descriptor = -1
            CALL pdgemr2d(n, n, fake_local_data(1, 1), 1, 1, fake_descriptor, &
                          eigenvectors%local_data(1, 1), 1, 1, eigenvectors%matrix_struct%descriptor, &
                          eigenvectors%matrix_struct%context%group)
         ENDIF
         ! free work
         CALL cp_fm_redistribute_work_finalize(work_redistribute%group_distribution(mepos_old) == 0)

         ! finally, also the eigenvalues need to end up on the non-group member tasks
         CALL mp_bcast(eig, 0, para_env%group)

      ELSE
         ! Just nullify pointers to the original matrices
         NULLIFY (matrix_new)
         NULLIFY (eigenvectors_new)

      ENDIF

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

END MODULE cp_fm_diag_utils
