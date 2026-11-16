!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for writing PAO restart files from MAO.
! **************************************************************************************************
MODULE mao_io
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_block_p,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mao_io'

   PUBLIC :: mao_write_pao_restart

   INTEGER, PARAMETER, PRIVATE :: file_format_version = 4

CONTAINS

! **************************************************************************************************
!> \brief Writes restart file
!> \param mao_coef ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE mao_write_pao_restart(mao_coef, qs_env)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coef
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: printkey_section = 'DFT%PRINT%MAO_ANALYSIS', &
         routineN = 'mao_write_pao_restart'

      INTEGER                                            :: handle, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      CALL get_qs_env(qs_env, input=input, para_env=para_env)

      IF (SIZE(mao_coef) == 1) THEN
         ! open file
         unit_nr = cp_print_key_unit_nr(logger, &
                                        input, &
                                        printkey_section, &
                                        extension=".pao", &
                                        file_action="WRITE", &
                                        file_position="REWIND", &
                                        file_status="UNKNOWN", &
                                        on_file=.TRUE.)

         IF (unit_nr > 0) CALL write_restart_header(qs_env, unit_nr)

         CALL write_diagonal_blocks(para_env, mao_coef(1)%matrix, "Xblock", unit_nr)
         IF (unit_nr > 0) WRITE (unit_nr, '(A)') "THE_END"

         CALL cp_print_key_finished_output(unit_nr, logger, input, printkey_section)
      ELSE
         CPWARN("MAO/PAO restart only for restricted case available.")
      END IF

      CALL timestop(handle)

   END SUBROUTINE mao_write_pao_restart

! **************************************************************************************************
!> \brief Write the digonal blocks of given DBCSR matrix into the provided unit_nr
!> \param para_env ...
!> \param matrix ...
!> \param label ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE write_diagonal_blocks(para_env, matrix, label, unit_nr)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type)                                   :: matrix
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      INTEGER, INTENT(IN)                                :: unit_nr

      INTEGER                                            :: i, iatom, m, n, natoms
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, row_blk_sizes
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:, :), POINTER                 :: local_block, mpi_buffer

      !TODO: this is a serial algorithm
      CALL dbcsr_get_info(matrix, row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes)
      CPASSERT(SIZE(row_blk_sizes) == SIZE(col_blk_sizes))
      natoms = SIZE(row_blk_sizes)

      DO iatom = 1, natoms
         n = row_blk_sizes(iatom)
         m = col_blk_sizes(iatom)
         ALLOCATE (mpi_buffer(n, m))
         NULLIFY (local_block)
         CALL dbcsr_get_block_p(matrix=matrix, row=iatom, col=iatom, block=local_block, found=found)
         IF (ASSOCIATED(local_block)) THEN
            IF (SIZE(local_block) > 0) THEN
               ! catch corner-case
               mpi_buffer(:, :) = local_block(:, :)
            END IF
         ELSE
            mpi_buffer(:, :) = 0.0_dp
         END IF

         CALL para_env%sum(mpi_buffer)
         IF (unit_nr > 0) THEN
            ! normalize vectors
            DO i = 1, m
               mpi_buffer(:, i) = mpi_buffer(:, i)/NORM2(mpi_buffer(:, i))
            END DO

            WRITE (unit_nr, fmt="(A,1X,I10,1X)", advance='no') label, iatom
            WRITE (unit_nr, *) mpi_buffer
         END IF
         DEALLOCATE (mpi_buffer)
      END DO

      ! flush
      IF (unit_nr > 0) FLUSH (unit_nr)

   END SUBROUTINE write_diagonal_blocks

! **************************************************************************************************
!> \brief Writes header of restart file
!> \param qs_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE write_restart_header(qs_env, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=default_string_length)               :: kindname
      INTEGER                                            :: iatom, ikind, ipot, istep, nmao, &
                                                            nparams, nsgf, z
      REAL(KIND=dp)                                      :: energy
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set)

      WRITE (unit_nr, "(A,5X,I0)") "Version", file_format_version
      energy = 0.0_dp
      WRITE (unit_nr, "(A,5X,F20.10)") "Energy", energy
      istep = 1
      WRITE (unit_nr, "(A,5X,I0)") "Step", istep
      WRITE (unit_nr, "(A,5X,A)") "Parametrization", "EQUIVARIANT"

      ! write kinds
      WRITE (unit_nr, "(A,5X,I0)") "Nkinds", SIZE(atomic_kind_set)
      DO ikind = 1, SIZE(atomic_kind_set)
         CALL get_atomic_kind(atomic_kind_set(ikind), name=kindname, z=z)
         CALL get_qs_kind(qs_kind_set(ikind), mao=nmao, nsgf=nsgf, basis_set=basis_set)
         nparams = nmao*nsgf
         WRITE (unit_nr, "(A,5X,I10,1X,A,1X,I3)") "Kind", ikind, TRIM(kindname), z
         WRITE (unit_nr, "(A,5X,I10,1X,I3)") "NParams", ikind, nparams
         WRITE (unit_nr, "(A,5X,I10,1X,I10,1X,A)") "PrimBasis", ikind, nsgf, TRIM(basis_set%name)
         WRITE (unit_nr, "(A,5X,I10,1X,I3)") "PaoBasis", ikind, nmao
         ipot = 0
         WRITE (unit_nr, "(A,5X,I10,1X,I3)") "NPaoPotentials", ikind, ipot
      END DO

      ! write cell
      WRITE (unit_nr, fmt="(A,5X)", advance='no') "Cell"
      WRITE (unit_nr, *) cell%hmat*angstrom

      ! write atoms
      WRITE (unit_nr, "(A,5X,I0)") "Natoms", SIZE(particle_set)
      DO iatom = 1, SIZE(particle_set)
         kindname = particle_set(iatom)%atomic_kind%name
         WRITE (unit_nr, fmt="(A,5X,I10,5X,A,1X)", advance='no') "Atom ", iatom, TRIM(kindname)
         WRITE (unit_nr, *) particle_set(iatom)%r*angstrom
      END DO

   END SUBROUTINE write_restart_header

END MODULE mao_io
