!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief methods related to the blacs parallel environment
!> \par History
!>      08.2002 created [fawzi]
!>      02.2004 modified to associate a blacs_env with a given para_env
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE cp_blacs_env
   USE cp_array_utils,                  ONLY: cp_2d_i_write
   USE cp_blacs_types,                  ONLY: cp_blacs_type
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush
   USE mathlib,                         ONLY: gcd
   USE message_passing,                 ONLY: mp_para_env_release,&
                                              mp_para_env_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_blacs_env'

   ! Blacs type of distribution
   INTEGER, PARAMETER, PUBLIC               :: BLACS_GRID_SQUARE = 1, &
                                               BLACS_GRID_ROW = 2, &
                                               BLACS_GRID_COL = 3

   PUBLIC :: cp_blacs_env_type
   PUBLIC :: cp_blacs_env_create, cp_blacs_env_release

! **************************************************************************************************
!> \brief represent a blacs multidimensional parallel environment
!>      (for the mpi corrispective see cp_paratypes/mp_para_cart_type)
!> \param ref_count the reference count, when it is zero this object gets
!>        deallocated
!> \param my_pid process id of the actual processor
!> \param n_pid number of process ids
!> \param the para_env associated (and compatible) with this blacs_env
!> \param blacs2mpi: maps mepos(1)-mepos(2) of blacs to its mpi rank
!> \param mpi2blacs(i,rank): maps the mpi rank to the mepos(i)
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   TYPE, EXTENDS(cp_blacs_type) :: cp_blacs_env_type
      INTEGER :: my_pid = -1, n_pid = -1, ref_count = -1
      TYPE(mp_para_env_type), POINTER :: para_env => NULL()
      INTEGER, DIMENSION(:, :), POINTER :: blacs2mpi => NULL()
      INTEGER, DIMENSION(:, :), POINTER :: mpi2blacs => NULL()
      LOGICAL :: repeatable = .FALSE.
   CONTAINS
      PROCEDURE, PUBLIC, PASS, NON_OVERRIDABLE :: create => cp_blacs_env_create_low
      PROCEDURE, PUBLIC, PASS, NON_OVERRIDABLE :: retain => cp_blacs_env_retain
      PROCEDURE, PUBLIC, PASS, NON_OVERRIDABLE :: release => cp_blacs_env_release_low
      PROCEDURE, PUBLIC, PASS, NON_OVERRIDABLE :: get => get_blacs_info
      PROCEDURE, PUBLIC, PASS, NON_OVERRIDABLE :: write => cp_blacs_env_write
   END TYPE cp_blacs_env_type

!***
CONTAINS

! **************************************************************************************************
!> \brief   Return informations about the specified BLACS context.
!> \param blacs_env ...
!> \param my_process_row ...
!> \param my_process_column ...
!> \param my_process_number ...
!> \param number_of_process_rows ...
!> \param number_of_process_columns ...
!> \param number_of_processes ...
!> \param para_env ...
!> \param blacs2mpi ...
!> \param mpi2blacs ...
!> \date    19.06.2001
!> \par     History
!>          MM.YYYY moved here from qs_blacs (Joost VandeVondele)
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE get_blacs_info(blacs_env, my_process_row, my_process_column, &
                             my_process_number, number_of_process_rows, &
                             number_of_process_columns, number_of_processes, &
                             para_env, blacs2mpi, mpi2blacs)
      CLASS(cp_blacs_env_type), INTENT(IN)                :: blacs_env
      INTEGER, INTENT(OUT), OPTIONAL :: my_process_row, my_process_column, my_process_number, &
                                        number_of_process_rows, number_of_process_columns, number_of_processes
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env
      INTEGER, DIMENSION(:, :), OPTIONAL, POINTER        :: blacs2mpi, mpi2blacs

      IF (PRESENT(my_process_row)) my_process_row = blacs_env%mepos(1)
      IF (PRESENT(my_process_column)) my_process_column = blacs_env%mepos(2)
      IF (PRESENT(my_process_number)) my_process_number = blacs_env%my_pid
      IF (PRESENT(number_of_process_rows)) number_of_process_rows = blacs_env%num_pe(1)
      IF (PRESENT(number_of_process_columns)) number_of_process_columns = blacs_env%num_pe(2)
      IF (PRESENT(number_of_processes)) number_of_processes = blacs_env%n_pid
      IF (PRESENT(para_env)) para_env => blacs_env%para_env
      IF (PRESENT(blacs2mpi)) blacs2mpi => blacs_env%blacs2mpi
      IF (PRESENT(mpi2blacs)) mpi2blacs => blacs_env%mpi2blacs

   END SUBROUTINE get_blacs_info

! **************************************************************************************************
!> \brief allocates and initializes a type that represent a blacs context
!> \param blacs_env the type to initialize
!> \param para_env the para_env for which a blacs env should be created
!> \param blacs_grid_layout ...
!> \param blacs_repeatable ...
!> \param row_major ...
!> \param grid_2d ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cp_blacs_env_create(blacs_env, para_env, blacs_grid_layout, blacs_repeatable, row_major, grid_2d)
      TYPE(cp_blacs_env_type), INTENT(OUT), POINTER      :: blacs_env
      TYPE(mp_para_env_type), INTENT(INOUT), TARGET      :: para_env
      INTEGER, INTENT(IN), OPTIONAL                      :: blacs_grid_layout
      LOGICAL, INTENT(IN), OPTIONAL                      :: blacs_repeatable, row_major
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: grid_2d

      ALLOCATE (blacs_env)
      CALL blacs_env%create(para_env, blacs_grid_layout, blacs_repeatable, row_major, grid_2d)

   END SUBROUTINE

! **************************************************************************************************
!> \brief allocates and initializes a type that represent a blacs context
!> \param blacs_env the type to initialize
!> \param para_env the para_env for which a blacs env should be created
!> \param blacs_grid_layout ...
!> \param blacs_repeatable ...
!> \param row_major ...
!> \param grid_2d ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cp_blacs_env_create_low(blacs_env, para_env, blacs_grid_layout, blacs_repeatable, row_major, grid_2d)
      CLASS(cp_blacs_env_type), INTENT(OUT)    :: blacs_env
      TYPE(mp_para_env_type), TARGET, INTENT(INOUT) :: para_env
      INTEGER, INTENT(IN), OPTIONAL            :: blacs_grid_layout
      LOGICAL, INTENT(IN), OPTIONAL            :: blacs_repeatable, row_major
      INTEGER, DIMENSION(:), INTENT(IN), &
         OPTIONAL                               :: grid_2d

      INTEGER                                  :: ipcol, iprow
#if defined(__parallel)
      INTEGER                                  :: gcd_max, ipe, jpe, &
                                                  my_blacs_grid_layout, &
                                                  npcol, npe, nprow
      LOGICAL                                  :: my_blacs_repeatable, &
                                                  my_row_major
#endif

#ifdef __parallel
      ! get the number of cpus for this blacs grid
      nprow = 1
      npcol = 1
      npe = para_env%num_pe
      ! get the layout of this grid

      IF (PRESENT(grid_2d)) THEN
         nprow = grid_2d(1)
         npcol = grid_2d(2)
      END IF

      IF (nprow*npcol .NE. npe) THEN
         ! hard code for the time being the grid layout
         my_blacs_grid_layout = BLACS_GRID_SQUARE
         IF (PRESENT(blacs_grid_layout)) my_blacs_grid_layout = blacs_grid_layout
         ! XXXXXX
         SELECT CASE (my_blacs_grid_layout)
         CASE (BLACS_GRID_SQUARE)
            ! make the grid as 'square' as possible, where square is defined as nprow and npcol
            ! having the largest possible gcd
            gcd_max = -1
            DO ipe = 1, CEILING(SQRT(REAL(npe, dp)))
               jpe = npe/ipe
               IF (ipe*jpe .NE. npe) CYCLE
               IF (gcd(ipe, jpe) >= gcd_max) THEN
                  nprow = ipe
                  npcol = jpe
                  gcd_max = gcd(ipe, jpe)
               END IF
            END DO
         CASE (BLACS_GRID_ROW)
            nprow = 1
            npcol = npe
         CASE (BLACS_GRID_COL)
            nprow = npe
            npcol = 1
         END SELECT
      END IF

      my_row_major = .TRUE.
      IF (PRESENT(row_major)) my_row_major = row_major
      IF (my_row_major) THEN
         CALL blacs_env%gridinit(para_env, "Row-major", nprow, npcol)
      ELSE
         CALL blacs_env%gridinit(para_env, "Col-major", nprow, npcol)
      END IF

      ! We set the components of blacs_env here such that we can still use INTENT(OUT) with gridinit
      blacs_env%my_pid = para_env%mepos
      blacs_env%n_pid = para_env%num_pe
      blacs_env%ref_count = 1

      my_blacs_repeatable = .FALSE.
      IF (PRESENT(blacs_repeatable)) my_blacs_repeatable = blacs_repeatable
      blacs_env%repeatable = my_blacs_repeatable
      IF (blacs_env%repeatable) CALL blacs_env%set(15, 1)

#else
      ! In serial mode, we just have to setup the object
      CALL blacs_env%gridinit(para_env, "Row-major", 1, 1)

      blacs_env%ref_count = 1
      blacs_env%my_pid = 0
      blacs_env%n_pid = 1
      MARK_USED(blacs_grid_layout)
      MARK_USED(blacs_repeatable)
      MARK_USED(grid_2d)
      MARK_USED(row_major)
#endif

      CALL para_env%retain()
      blacs_env%para_env => para_env

      ! generate the mappings blacs2mpi and mpi2blacs
      ALLOCATE (blacs_env%blacs2mpi(0:blacs_env%num_pe(1) - 1, 0:blacs_env%num_pe(2) - 1))
      blacs_env%blacs2mpi = 0
      blacs_env%blacs2mpi(blacs_env%mepos(1), blacs_env%mepos(2)) = para_env%mepos
      CALL para_env%sum(blacs_env%blacs2mpi)
      ALLOCATE (blacs_env%mpi2blacs(2, 0:para_env%num_pe - 1))
      blacs_env%mpi2blacs = -1
      DO ipcol = 0, blacs_env%num_pe(2) - 1
         DO iprow = 0, blacs_env%num_pe(1) - 1
            blacs_env%mpi2blacs(1, blacs_env%blacs2mpi(iprow, ipcol)) = iprow
            blacs_env%mpi2blacs(2, blacs_env%blacs2mpi(iprow, ipcol)) = ipcol
         END DO
      END DO
   END SUBROUTINE cp_blacs_env_create_low

! **************************************************************************************************
!> \brief retains the given blacs env
!> \param blacs_env the blacs env to retain
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cp_blacs_env_retain(blacs_env)
      CLASS(cp_blacs_env_type), INTENT(INOUT)            :: blacs_env

      CPASSERT(blacs_env%ref_count > 0)
      blacs_env%ref_count = blacs_env%ref_count + 1
   END SUBROUTINE cp_blacs_env_retain

! **************************************************************************************************
!> \brief releases the given blacs_env
!> \param blacs_env the blacs env to release
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cp_blacs_env_release(blacs_env)
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      IF (ASSOCIATED(blacs_env)) THEN
         CPASSERT(blacs_env%ref_count > 0)
         blacs_env%ref_count = blacs_env%ref_count - 1
         IF (blacs_env%ref_count < 1) THEN
            CALL blacs_env%release()
            DEALLOCATE (blacs_env)
         END IF
      END IF
      NULLIFY (blacs_env)
   END SUBROUTINE cp_blacs_env_release

! **************************************************************************************************
!> \brief releases the given blacs_env
!> \param blacs_env the blacs env to release
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cp_blacs_env_release_low(blacs_env)
      CLASS(cp_blacs_env_type), INTENT(INOUT)                   :: blacs_env

      CALL blacs_env%gridexit()
      CALL mp_para_env_release(blacs_env%para_env)
      DEALLOCATE (blacs_env%mpi2blacs)
      DEALLOCATE (blacs_env%blacs2mpi)

   END SUBROUTINE cp_blacs_env_release_low

! **************************************************************************************************
!> \brief writes the description of the given blacs env
!> \param blacs_env the blacs environment to write
!> \param unit_nr the unit number where to write the description of the
!>        blacs environment
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE cp_blacs_env_write(blacs_env, unit_nr)
      CLASS(cp_blacs_env_type), INTENT(IN)                :: blacs_env
      INTEGER, INTENT(in)                                :: unit_nr

      WRITE (unit=unit_nr, fmt="('  group=',i10,', ref_count=',i10,',')") &
         blacs_env%get_handle(), blacs_env%ref_count
      WRITE (unit=unit_nr, fmt="('  mepos=(',i8,',',i8,'),')") &
         blacs_env%mepos(1), blacs_env%mepos(2)
      WRITE (unit=unit_nr, fmt="('  num_pe=(',i8,',',i8,'),')") &
         blacs_env%num_pe(1), blacs_env%num_pe(2)
      IF (ASSOCIATED(blacs_env%blacs2mpi)) THEN
         WRITE (unit=unit_nr, fmt="('  blacs2mpi=')", advance="no")
         CALL cp_2d_i_write(blacs_env%blacs2mpi, unit_nr=unit_nr)
      ELSE
         WRITE (unit=unit_nr, fmt="('  blacs2mpi=*null*')")
      END IF
      IF (ASSOCIATED(blacs_env%para_env)) THEN
         WRITE (unit=unit_nr, fmt="('  para_env=<cp_para_env id=',i6,'>,')") &
            blacs_env%para_env%get_handle()
      ELSE
         WRITE (unit=unit_nr, fmt="('  para_env=*null*')")
      END IF
      WRITE (unit=unit_nr, fmt="('  my_pid=',i10,', n_pid=',i10,' }')") &
         blacs_env%my_pid, blacs_env%n_pid
      CALL m_flush(unit_nr)
   END SUBROUTINE cp_blacs_env_write

END MODULE cp_blacs_env
