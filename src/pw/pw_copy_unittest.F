!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM pw_copy_unittest
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE pw_copy_all,                     ONLY: pw_copy_match
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_types,                        ONLY: pw_c1d_gs_type

   IMPLICIT NONE

   INTEGER                                            :: mode
   TYPE(mp_para_env_type), POINTER                     :: para_env

   ALLOCATE (para_env)
   CALL mp_world_init(para_env)
   DO mode = 0, 8
      CALL check_copy(mode)
   END DO
   DEALLOCATE (para_env)
   CALL mp_world_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Compare identical, permuted, redistributed, partial and empty PW grids.
!> \param mode Layout variant; mode two changes only one rank's ordering.
!> \author Thomas D. Kuehne (tkuehne@cp2k.org)
! **************************************************************************************************
   SUBROUTINE check_copy(mode)
      INTEGER, INTENT(IN)                                :: mode

      COMPLEX(dp), PARAMETER :: untouched = (-999.0_dp, 17.0_dp)

      COMPLEX(dp)                                        :: expected
      INTEGER                                            :: dimensions(1), i, j, nsource, ntarget, &
                                                            owner, rank
      TYPE(pw_c1d_gs_type)                               :: destination, source
      TYPE(pw_grid_type), TARGET                         :: source_grid, target_grid

      rank = para_env%mepos
      nsource = 7
      IF (mode == 6 .OR. (mode >= 7 .AND. rank == 0)) nsource = 0
      owner = rank
      IF (mode == 3 .OR. mode == 8) owner = MOD(rank + 1, para_env%num_pe)
      ntarget = nsource
      IF (mode == 4) ntarget = nsource - 1
      IF (mode == 5) ntarget = nsource + 1
      IF (mode == 8) ntarget = MERGE(0, 7, owner == 0)
      ALLOCATE (source_grid%g_hat(3, nsource), target_grid%g_hat(3, ntarget))
      ALLOCATE (source%array(nsource), destination%array(ntarget))
      source%pw_grid => source_grid
      destination%pw_grid => target_grid
      dimensions = [para_env%num_pe]
      CALL source_grid%para%group%create(para_env, 1, dimensions)
      target_grid%para%group = source_grid%para%group
      source_grid%ngpts_cut_local = nsource
      target_grid%ngpts_cut_local = ntarget
      DO i = 1, nsource
         source_grid%g_hat(:, i) = [i, rank, 0]
         source%array(i) = CMPLX(i + 11*rank, -i - 3*rank, dp)
      END DO
      DO i = 1, ntarget
         j = i
         IF (mode == 1 .OR. (mode == 2 .AND. rank == 0)) j = ntarget + 1 - i
         target_grid%g_hat(:, i) = [j, owner, 0]
      END DO
      destination%array(:) = untouched
      CALL pw_copy_match(source, destination)
      DO i = 1, ntarget
         j = target_grid%g_hat(1, i)
         expected = CMPLX(j + 11*owner, -j - 3*owner, dp)
         IF (mode == 5 .AND. j > nsource) expected = untouched
         IF (ABS(destination%array(i) - expected) > 1.e-14_dp) ERROR STOP 'Incorrect PW copy'
      END DO
      DEALLOCATE (source%array, destination%array)
      CALL source_grid%para%group%free()
   END SUBROUTINE check_copy

END PROGRAM pw_copy_unittest
