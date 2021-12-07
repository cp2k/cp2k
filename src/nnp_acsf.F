!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Functionality for atom centered symmetry functions
!>         for neural network potentials
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
MODULE nnp_acsf
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: pi
   USE nnp_environment_types,           ONLY: nnp_acsf_ang_type,&
                                              nnp_acsf_rad_type,&
                                              nnp_cut_cos,&
                                              nnp_cut_tanh,&
                                              nnp_env_get,&
                                              nnp_neighbor_type,&
                                              nnp_type
   USE periodic_table,                  ONLY: get_ptable_info
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nnp_acsf'
   ! Public subroutines ***
   PUBLIC :: nnp_calc_acsf, &
             nnp_init_acsf_groups, &
             nnp_sort_acsf, &
             nnp_sort_ele, &
             nnp_write_acsf

CONTAINS

! **************************************************************************************************
!> \brief Calculate atom centered symmetry functions for given atom i
!> \param nnp ...
!> \param i ...
!> \param dsymdxyz ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_calc_acsf(nnp, i, dsymdxyz)

      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      INTEGER, INTENT(IN)                                :: i
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dsymdxyz

      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_calc_acsf', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, handle_sf, ind, j, jj, k, kk, l, &
                                                            m, off, s, sf
      REAL(KIND=dp)                                      :: r1, r2, r3
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: symtmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: forcetmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: force3tmp
      REAL(KIND=dp), DIMENSION(3)                        :: rvect1, rvect2, rvect3
      TYPE(nnp_neighbor_type)                            :: neighbor

      CALL timeset(routineN, handle)

      !determine index of atom type
      ind = nnp%ele_ind(i)

      ! compute neighbors of atom i
      CALL nnp_neighbor_create(nnp, ind, nnp%num_atoms, neighbor)
      CALL nnp_compute_neighbors(nnp, neighbor, i)

      ! Reset y:
      nnp%rad(ind)%y = 0.0_dp
      nnp%ang(ind)%y = 0.0_dp

      !calc forces
      IF (PRESENT(dsymdxyz)) THEN
         !loop over radial sym fnct grps
         CALL timeset('nnp_acsf_radial', handle_sf)
         DO s = 1, nnp%rad(ind)%n_symfgrp
            ALLOCATE (symtmp(nnp%rad(ind)%symfgrp(s)%n_symf))
            ALLOCATE (forcetmp(3, nnp%rad(ind)%symfgrp(s)%n_symf))
            !loop over associated neighbors
            DO j = 1, neighbor%n_rad(s)
               rvect1 = neighbor%dist_rad(1:3, j, s)
               r1 = neighbor%dist_rad(4, j, s)
               CALL nnp_calc_rad(nnp, ind, s, rvect1, r1, symtmp, forcetmp)
               jj = neighbor%ind_rad(j, s)
               DO sf = 1, nnp%rad(ind)%symfgrp(s)%n_symf
                  m = nnp%rad(ind)%symfgrp(s)%symf(sf)
                  ! update forces into dsymdxyz
                  DO l = 1, 3
                     dsymdxyz(l, m, i) = dsymdxyz(l, m, i) + forcetmp(l, sf)
                     dsymdxyz(l, m, jj) = dsymdxyz(l, m, jj) - forcetmp(l, sf)
                  END DO
                  nnp%rad(ind)%y(m) = nnp%rad(ind)%y(m) + symtmp(sf)
               END DO
            END DO
            DEALLOCATE (symtmp)
            DEALLOCATE (forcetmp)
         END DO
         CALL timestop(handle_sf)
         !loop over angular sym fnct grps
         CALL timeset('nnp_acsf_angular', handle_sf)
         off = nnp%n_rad(ind)
         DO s = 1, nnp%ang(ind)%n_symfgrp
            ALLOCATE (symtmp(nnp%ang(ind)%symfgrp(s)%n_symf))
            ALLOCATE (force3tmp(3, 3, nnp%ang(ind)%symfgrp(s)%n_symf))
            !loop over associated neighbors
            IF (nnp%ang(ind)%symfgrp(s)%ele(1) == nnp%ang(ind)%symfgrp(s)%ele(2)) THEN
               DO j = 1, neighbor%n_ang1(s)
                  rvect1 = neighbor%dist_ang1(1:3, j, s)
                  r1 = neighbor%dist_ang1(4, j, s)
                  DO k = j + 1, neighbor%n_ang1(s)
                     rvect2 = neighbor%dist_ang1(1:3, k, s)
                     r2 = neighbor%dist_ang1(4, k, s)
                     rvect3 = rvect2 - rvect1
                     r3 = NORM2(rvect3(:))
                     IF (r3 < nnp%ang(ind)%symfgrp(s)%cutoff) THEN
                        jj = neighbor%ind_ang1(j, s)
                        kk = neighbor%ind_ang1(k, s)
                        CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, &
                                          r1, r2, r3, symtmp, force3tmp)
                        DO sf = 1, nnp%ang(ind)%symfgrp(s)%n_symf
                           m = off + nnp%ang(ind)%symfgrp(s)%symf(sf)
                           ! update forces into dsymdxy
                           DO l = 1, 3
                              dsymdxyz(l, m, i) = dsymdxyz(l, m, i) &
                                                  + force3tmp(l, 1, sf)
                              dsymdxyz(l, m, jj) = dsymdxyz(l, m, jj) &
                                                   + force3tmp(l, 2, sf)
                              dsymdxyz(l, m, kk) = dsymdxyz(l, m, kk) &
                                                   + force3tmp(l, 3, sf)
                           END DO
                           nnp%ang(ind)%y(m - off) = nnp%ang(ind)%y(m - off) + symtmp(sf)
                        END DO
                     END IF
                  END DO
               END DO
            ELSE
               DO j = 1, neighbor%n_ang1(s)
                  rvect1 = neighbor%dist_ang1(1:3, j, s)
                  r1 = neighbor%dist_ang1(4, j, s)
                  DO k = 1, neighbor%n_ang2(s)
                     rvect2 = neighbor%dist_ang2(1:3, k, s)
                     r2 = neighbor%dist_ang2(4, k, s)
                     rvect3 = rvect2 - rvect1
                     r3 = NORM2(rvect3(:))
                     IF (r3 < nnp%ang(ind)%symfgrp(s)%cutoff) THEN
                        jj = neighbor%ind_ang1(j, s)
                        kk = neighbor%ind_ang1(k, s)
                        CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, &
                                          r1, r2, r3, symtmp, force3tmp)
                        !loop over associated sym fncts
                        DO sf = 1, nnp%ang(ind)%symfgrp(s)%n_symf
                           m = off + nnp%ang(ind)%symfgrp(s)%symf(sf)
                           jj = neighbor%ind_ang1(j, s)
                           kk = neighbor%ind_ang2(k, s)
                           ! update forces into dsymdxy
                           DO l = 1, 3
                              dsymdxyz(l, m, i) = dsymdxyz(l, m, i) &
                                                  + force3tmp(l, 1, sf)
                              dsymdxyz(l, m, jj) = dsymdxyz(l, m, jj) &
                                                   + force3tmp(l, 2, sf)
                              dsymdxyz(l, m, kk) = dsymdxyz(l, m, kk) &
                                                   + force3tmp(l, 3, sf)
                           END DO
                           nnp%ang(ind)%y(m - off) = nnp%ang(ind)%y(m - off) + symtmp(sf)
                        END DO
                     END IF
                  END DO
               END DO
            END IF
            DEALLOCATE (symtmp)
            DEALLOCATE (force3tmp)
         END DO
         CALL timestop(handle_sf)
         !no forces
      ELSE
         !loop over radial sym fnct grps
         CALL timeset('nnp_acsf_radial', handle_sf)
         DO s = 1, nnp%rad(ind)%n_symfgrp
            ALLOCATE (symtmp(nnp%rad(ind)%symfgrp(s)%n_symf))
            !loop over associated neighbors
            DO j = 1, neighbor%n_rad(s)
               rvect1 = neighbor%dist_rad(1:3, j, s)
               r1 = neighbor%dist_rad(4, j, s)
               CALL nnp_calc_rad(nnp, ind, s, rvect1, r1, symtmp)
               DO sf = 1, nnp%rad(ind)%symfgrp(s)%n_symf
                  m = nnp%rad(ind)%symfgrp(s)%symf(sf)
                  nnp%rad(ind)%y(m) = nnp%rad(ind)%y(m) + symtmp(sf)
               END DO
            END DO
            DEALLOCATE (symtmp)
         END DO
         CALL timestop(handle_sf)
         !loop over angular sym fnct grps
         CALL timeset('nnp_acsf_angular', handle_sf)
         off = nnp%n_rad(ind)
         DO s = 1, nnp%ang(ind)%n_symfgrp
            ALLOCATE (symtmp(nnp%ang(ind)%symfgrp(s)%n_symf))
            !loop over associated neighbors
            IF (nnp%ang(ind)%symfgrp(s)%ele(1) == nnp%ang(ind)%symfgrp(s)%ele(2)) THEN
               DO j = 1, neighbor%n_ang1(s)
                  rvect1 = neighbor%dist_ang1(1:3, j, s)
                  r1 = neighbor%dist_ang1(4, j, s)
                  DO k = j + 1, neighbor%n_ang1(s)
                     rvect2 = neighbor%dist_ang1(1:3, k, s)
                     r2 = neighbor%dist_ang1(4, k, s)
                     rvect3 = rvect2 - rvect1
                     r3 = NORM2(rvect3(:))
                     IF (r3 < nnp%ang(ind)%symfgrp(s)%cutoff) THEN
                        CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, r1, r2, r3, symtmp)
                        DO sf = 1, nnp%ang(ind)%symfgrp(s)%n_symf
                           m = off + nnp%ang(ind)%symfgrp(s)%symf(sf)
                           nnp%ang(ind)%y(m - off) = nnp%ang(ind)%y(m - off) + symtmp(sf)
                        END DO
                     END IF
                  END DO
               END DO
            ELSE
               DO j = 1, neighbor%n_ang1(s)
                  rvect1 = neighbor%dist_ang1(1:3, j, s)
                  r1 = neighbor%dist_ang1(4, j, s)
                  DO k = 1, neighbor%n_ang2(s)
                     rvect2 = neighbor%dist_ang2(1:3, k, s)
                     r2 = neighbor%dist_ang2(4, k, s)
                     rvect3 = rvect2 - rvect1
                     r3 = NORM2(rvect3(:))
                     IF (r3 < nnp%ang(ind)%symfgrp(s)%cutoff) THEN
                        CALL nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, r1, r2, r3, symtmp)
                        !loop over associated sym fncts
                        DO sf = 1, nnp%ang(ind)%symfgrp(s)%n_symf
                           m = off + nnp%ang(ind)%symfgrp(s)%symf(sf)
                           nnp%ang(ind)%y(m - off) = nnp%ang(ind)%y(m - off) + symtmp(sf)
                        END DO
                     END IF
                  END DO
               END DO
            END IF
            DEALLOCATE (symtmp)
         END DO
         CALL timestop(handle_sf)
      END IF

      !check extrapolation
      CALL nnp_check_extrapolation(nnp, ind)

      IF (PRESENT(dsymdxyz)) THEN
         CALL nnp_scale_acsf(nnp, ind, dsymdxyz)
      ELSE
         CALL nnp_scale_acsf(nnp, ind)
      END IF

      CALL nnp_neighbor_release(neighbor)
      CALL timestop(handle)

   END SUBROUTINE nnp_calc_acsf

! **************************************************************************************************
!> \brief Check if the nnp is extrapolating
!> \param nnp ...
!> \param ind ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_check_extrapolation(nnp, ind)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: ind

      REAL(KIND=dp), PARAMETER                           :: threshold = 0.0001_dp

      INTEGER                                            :: j
      LOGICAL                                            :: extrapolate

      extrapolate = .FALSE.

      DO j = 1, nnp%n_rad(ind)
         IF (nnp%rad(ind)%y(j) - &
             nnp%rad(ind)%loc_max(j) > threshold) THEN
            extrapolate = .TRUE.
         ELSE IF (-nnp%rad(ind)%y(j) + &
                  nnp%rad(ind)%loc_min(j) > threshold) THEN
            extrapolate = .TRUE.
         END IF
      END DO
      DO j = 1, nnp%n_ang(ind)
         IF (nnp%ang(ind)%y(j) - &
             nnp%ang(ind)%loc_max(j) > threshold) THEN
            extrapolate = .TRUE.
         ELSE IF (-nnp%ang(ind)%y(j) + &
                  nnp%ang(ind)%loc_min(j) > threshold) THEN
            extrapolate = .TRUE.
         END IF
      END DO

      nnp%output_expol = extrapolate

   END SUBROUTINE nnp_check_extrapolation

! **************************************************************************************************
!> \brief Scale and center symetry functions (and gradients)
!> \param nnp ...
!> \param ind ...
!> \param dsymdxyz ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_scale_acsf(nnp, ind, dsymdxyz)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: ind
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT), &
         OPTIONAL                                        :: dsymdxyz

      INTEGER                                            :: j, k, off

      IF (nnp%center_acsf) THEN
         DO j = 1, nnp%n_rad(ind)
            nnp%arc(ind)%layer(1)%node(j) = &
               (nnp%rad(ind)%y(j) - nnp%rad(ind)%loc_av(j))
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            nnp%arc(ind)%layer(1)%node(j + off) = &
               (nnp%ang(ind)%y(j) - nnp%ang(ind)%loc_av(j))
         END DO

         IF (nnp%scale_acsf) THEN
            DO j = 1, nnp%n_rad(ind)
               nnp%arc(ind)%layer(1)%node(j) = &
                  nnp%arc(ind)%layer(1)%node(j)/ &
                  (nnp%rad(ind)%loc_max(j) - nnp%rad(ind)%loc_min(j))* &
                  (nnp%scmax - nnp%scmin) + nnp%scmin
            END DO
            off = nnp%n_rad(ind)
            DO j = 1, nnp%n_ang(ind)
               nnp%arc(ind)%layer(1)%node(j + off) = &
                  nnp%arc(ind)%layer(1)%node(j + off)/ &
                  (nnp%ang(ind)%loc_max(j) - nnp%ang(ind)%loc_min(j))* &
                  (nnp%scmax - nnp%scmin) + nnp%scmin
            END DO
         END IF
      ELSE IF (nnp%scale_acsf) THEN
         DO j = 1, nnp%n_rad(ind)
            nnp%arc(ind)%layer(1)%node(j) = &
               (nnp%rad(ind)%y(j) - nnp%rad(ind)%loc_min(j))/ &
               (nnp%rad(ind)%loc_max(j) - nnp%rad(ind)%loc_min(j))* &
               (nnp%scmax - nnp%scmin) + nnp%scmin
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            nnp%arc(ind)%layer(1)%node(j + off) = &
               (nnp%ang(ind)%y(j) - nnp%ang(ind)%loc_min(j))/ &
               (nnp%ang(ind)%loc_max(j) - nnp%ang(ind)%loc_min(j))* &
               (nnp%scmax - nnp%scmin) + nnp%scmin
         END DO
      ELSE IF (nnp%scale_sigma_acsf) THEN
         DO j = 1, nnp%n_rad(ind)
            nnp%arc(ind)%layer(1)%node(j) = &
               (nnp%rad(ind)%y(j) - nnp%rad(ind)%loc_av(j))/ &
               nnp%rad(ind)%sigma(j)* &
               (nnp%scmax - nnp%scmin) + nnp%scmin
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            nnp%arc(ind)%layer(1)%node(j + off) = &
               (nnp%ang(ind)%y(j) - nnp%ang(ind)%loc_av(j))/ &
               nnp%ang(ind)%sigma(j)* &
               (nnp%scmax - nnp%scmin) + nnp%scmin
         END DO
      ELSE
         DO j = 1, nnp%n_rad(ind)
            nnp%arc(ind)%layer(1)%node(j) = nnp%rad(ind)%y(j)
         END DO
         off = nnp%n_rad(ind)
         DO j = 1, nnp%n_ang(ind)
            nnp%arc(ind)%layer(1)%node(j + off) = nnp%ang(ind)%y(j)
         END DO
      END IF

      IF (PRESENT(dsymdxyz)) THEN
         IF (nnp%scale_acsf) THEN
            DO k = 1, nnp%num_atoms
               DO j = 1, nnp%n_rad(ind)
                  dsymdxyz(:, j, k) = dsymdxyz(:, j, k)/ &
                                      (nnp%rad(ind)%loc_max(j) - &
                                       nnp%rad(ind)%loc_min(j))* &
                                      (nnp%scmax - nnp%scmin)
               END DO
            END DO
            off = nnp%n_rad(ind)
            DO k = 1, nnp%num_atoms
               DO j = 1, nnp%n_ang(ind)
                  dsymdxyz(:, j + off, k) = dsymdxyz(:, j + off, k)/ &
                                            (nnp%ang(ind)%loc_max(j) - &
                                             nnp%ang(ind)%loc_min(j))* &
                                            (nnp%scmax - nnp%scmin)
               END DO
            END DO
         ELSE IF (nnp%scale_sigma_acsf) THEN
            DO k = 1, nnp%num_atoms
               DO j = 1, nnp%n_rad(ind)
                  dsymdxyz(:, j, k) = dsymdxyz(:, j, k)/ &
                                      nnp%rad(ind)%sigma(j)* &
                                      (nnp%scmax - nnp%scmin)
               END DO
            END DO
            off = nnp%n_rad(ind)
            DO k = 1, nnp%num_atoms
               DO j = 1, nnp%n_ang(ind)
                  dsymdxyz(:, j + off, k) = dsymdxyz(:, j + off, k)/ &
                                            nnp%ang(ind)%sigma(j)* &
                                            (nnp%scmax - nnp%scmin)
               END DO
            END DO
         END IF
      END IF

   END SUBROUTINE nnp_scale_acsf

! **************************************************************************************************
!> \brief Calculate radial symmetry function and gradient (optinal)
!>        for given displacment vecotr rvect of atom i and j
!> \param nnp ...
!> \param ind ...
!> \param s ...
!> \param rvect ...
!> \param r ...
!> \param sym ...
!> \param force ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_calc_rad(nnp, ind, s, rvect, r, sym, force)

      TYPE(nnp_type), INTENT(IN)                         :: nnp
      INTEGER, INTENT(IN)                                :: ind, s
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rvect
      REAL(KIND=dp), INTENT(IN)                          :: r
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: sym
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: force

      INTEGER                                            :: k, sf
      REAL(KIND=dp)                                      :: dfcutdr, dsymdr, eta, fcut, funccut, rs, &
                                                            tmp
      REAL(KIND=dp), DIMENSION(3)                        :: drdx

      !init
      drdx = 0.0_dp
      fcut = 0.0_dp
      dfcutdr = 0.0_dp

      !Calculate cutoff function and partial derivative
      funccut = nnp%rad(ind)%symfgrp(s)%cutoff !cutoff
      SELECT CASE (nnp%cut_type)
      CASE (nnp_cut_cos)
         tmp = pi*r/funccut
         fcut = 0.5_dp*(COS(tmp) + 1.0_dp)
         IF (PRESENT(force)) THEN
            dfcutdr = 0.5_dp*(-SIN(tmp))*(pi/funccut)
         END IF
      CASE (nnp_cut_tanh)
         tmp = TANH(1.0_dp - r/funccut)
         fcut = tmp**3
         IF (PRESENT(force)) THEN
            dfcutdr = (-3.0_dp/funccut)*(tmp**2 - tmp**4)
         END IF
      CASE DEFAULT
         CPABORT("NNP| Cutoff function unknown")
      END SELECT

      IF (PRESENT(force)) drdx(:) = rvect(:)/r

      !loop over sym fncts of sym fnct group s
      DO sf = 1, nnp%rad(ind)%symfgrp(s)%n_symf
         k = nnp%rad(ind)%symfgrp(s)%symf(sf) !symf indice
         eta = nnp%rad(ind)%eta(k) !eta
         rs = nnp%rad(ind)%rs(k) !rshift

         ! Calculate radial symmetry function
         sym(sf) = EXP(-eta*(r - rs)**2)

         ! Calculate partial derivatives of symmetry function and distance
         ! and combine them to obtain force
         IF (PRESENT(force)) THEN
            dsymdr = sym(sf)*(-2.0_dp*eta*(r - rs))
            force(:, sf) = fcut*dsymdr*drdx(:) + sym(sf)*dfcutdr*drdx(:)
         END IF

         ! combine radial symmetry function and cutoff function
         sym(sf) = sym(sf)*fcut
      END DO

   END SUBROUTINE nnp_calc_rad

! **************************************************************************************************
!> \brief Calculate angular symmetry function and gradient (optinal)
!>        for given displacment vectors rvect1,2,3 of atom i,j and k
!> \param nnp ...
!> \param ind ...
!> \param s ...
!> \param rvect1 ...
!> \param rvect2 ...
!> \param rvect3 ...
!> \param r1 ...
!> \param r2 ...
!> \param r3 ...
!> \param sym ...
!> \param force ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_calc_ang(nnp, ind, s, rvect1, rvect2, rvect3, r1, r2, r3, sym, force)

      TYPE(nnp_type), INTENT(IN)                         :: nnp
      INTEGER, INTENT(IN)                                :: ind, s
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rvect1, rvect2, rvect3
      REAL(KIND=dp), INTENT(IN)                          :: r1, r2, r3
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: sym
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT), &
         OPTIONAL                                        :: force

      INTEGER                                            :: i, m, sf
      REAL(KIND=dp) :: angular, costheta, dfcutdr1, dfcutdr2, dfcutdr3, dsymdr1, dsymdr2, dsymdr3, &
         eta, f, fcut, fcut1, fcut2, fcut3, ftot, g, lam, pref, prefzeta, rsqr1, rsqr2, rsqr3, &
         symtmp, tmp, tmp1, tmp2, tmp3, tmpzeta, zeta
      REAL(KIND=dp), DIMENSION(3) :: dangulardx1, dangulardx2, dangulardx3, dcosthetadx1, &
         dcosthetadx2, dcosthetadx3, dfdx1, dfdx2, dfdx3, dgdx1, dgdx2, dgdx3, dr1dx, dr2dx, dr3dx

      rsqr1 = r1**2
      rsqr2 = r2**2
      rsqr3 = r3**2

      !init
      dangulardx1 = 0.0_dp
      dangulardx2 = 0.0_dp
      dangulardx3 = 0.0_dp
      dr1dx = 0.0_dp
      dr2dx = 0.0_dp
      dr3dx = 0.0_dp
      ftot = 0.0_dp
      dfcutdr1 = 0.0_dp
      dfcutdr2 = 0.0_dp
      dfcutdr3 = 0.0_dp

      ! Calculate cos(theta)
      ! use law of cosine for theta
      ! cos(ang(r1,r2)) = (r3**2 - r1**2 - r2**2)/(-2*r1*r2)
      !                   |          f           |    g    |
      f = (rsqr3 - rsqr1 - rsqr2)
      g = -2.0_dp*r1*r2
      costheta = f/g

      ! Calculate cutoff function and partial derivatives
      fcut = nnp%ang(ind)%symfgrp(s)%cutoff !cutoff
      SELECT CASE (nnp%cut_type)
      CASE (nnp_cut_cos)
         tmp1 = pi*r1/fcut
         tmp2 = pi*r2/fcut
         tmp3 = pi*r3/fcut
         fcut1 = 0.5_dp*(COS(tmp1) + 1.0_dp)
         fcut2 = 0.5_dp*(COS(tmp2) + 1.0_dp)
         fcut3 = 0.5_dp*(COS(tmp3) + 1.0_dp)
         ftot = fcut1*fcut2*fcut3
         IF (PRESENT(force)) THEN
            pref = 0.5_dp*(pi/fcut)
            dfcutdr1 = pref*(-SIN(tmp1))*fcut2*fcut3
            dfcutdr2 = pref*(-SIN(tmp2))*fcut1*fcut3
            dfcutdr3 = pref*(-SIN(tmp3))*fcut1*fcut2
         END IF
      CASE (nnp_cut_tanh)
         tmp1 = TANH(1.0_dp - r1/fcut)
         tmp2 = TANH(1.0_dp - r2/fcut)
         tmp3 = TANH(1.0_dp - r3/fcut)
         fcut1 = tmp1**3
         fcut2 = tmp2**3
         fcut3 = tmp3**3
         ftot = fcut1*fcut2*fcut3
         IF (PRESENT(force)) THEN
            pref = -3.0_dp/fcut
            dfcutdr1 = pref*(tmp1**2 - tmp1**4)*fcut2*fcut3
            dfcutdr2 = pref*(tmp2**2 - tmp2**4)*fcut1*fcut3
            dfcutdr3 = pref*(tmp3**2 - tmp3**4)*fcut1*fcut2
         END IF
      CASE DEFAULT
         CPABORT("NNP| Cutoff function unknown")
      END SELECT

      IF (PRESENT(force)) THEN
         dr1dx(:) = rvect1(:)/r1
         dr2dx(:) = rvect2(:)/r2
         dr3dx(:) = rvect3(:)/r3
      END IF

      !loop over associated sym fncts
      DO sf = 1, nnp%ang(ind)%symfgrp(s)%n_symf
         m = nnp%ang(ind)%symfgrp(s)%symf(sf)
         lam = nnp%ang(ind)%lam(m) !lambda
         zeta = nnp%ang(ind)%zeta(m) !zeta
         prefzeta = nnp%ang(ind)%prefzeta(m) ! 2**(1-zeta)
         eta = nnp%ang(ind)%eta(m) !eta

         tmp = (1.0_dp + lam*costheta)

         IF (tmp <= 0.0_dp) THEN
            sym(sf) = 0.0_dp
            IF (PRESENT(force)) force(:, :, sf) = 0.0_dp
            CYCLE
         END IF

         ! Calculate symmetry function
         ! Calculate angular symmetry function
         ! ang = (1+lam*cos(theta_ijk))**zeta
         i = NINT(zeta)
         IF (1.0_dp*i == zeta) THEN
            tmpzeta = tmp**(i - 1)   ! integer power is a LOT faster
         ELSE
            tmpzeta = tmp**(zeta - 1.0_dp)
         END IF
         angular = tmpzeta*tmp
         ! exponential part:
         ! exp(-eta*(r1**2+r2**2+r3**2))
         symtmp = EXP(-eta*(rsqr1 + rsqr2 + rsqr3))

         ! Partial derivatives (f/g)' = (f'g - fg')/g^2
         IF (PRESENT(force)) THEN
            pref = zeta*(tmpzeta)/g**2
            DO i = 1, 3
               dfdx1(i) = -2.0_dp*lam*(rvect1(i) + rvect2(i))
               dfdx2(i) = 2.0_dp*lam*(rvect3(i) + rvect1(i))
               dfdx3(i) = 2.0_dp*lam*(rvect2(i) - rvect3(i))

               tmp1 = 2.0_dp*r2*dr1dx(i)
               tmp2 = 2.0_dp*r1*dr2dx(i)
               dgdx1(i) = -(tmp1 + tmp2)
               dgdx2(i) = tmp1
               dgdx3(i) = tmp2

               dcosthetadx1(i) = dfdx1(i)*g - lam*f*dgdx1(i)
               dcosthetadx2(i) = dfdx2(i)*g - lam*f*dgdx2(i)
               dcosthetadx3(i) = dfdx3(i)*g - lam*f*dgdx3(i)

               dangulardx1(i) = pref*dcosthetadx1(i)
               dangulardx2(i) = pref*dcosthetadx2(i)
               dangulardx3(i) = pref*dcosthetadx3(i)
            END DO

            ! Calculate partial derivatives of exponential part and distance
            ! and combine partial derivatives to obtain force
            pref = prefzeta
            tmp = -2.0_dp*symtmp*eta
            dsymdr1 = tmp*r1
            dsymdr2 = tmp*r2
            dsymdr3 = tmp*r3

            ! G(r1,r2,r3) = pref*angular(r1,r2,r3)*exp(r1,r2,r3)*fcut(r1,r2,r3)
            ! dG/dx1 = (dangular/dx1*  exp    *    fcut   +
            !            angular    * dexp/dx1*    fcut   +
            !            angular    *  exp    *  dfcut/dx1)*pref
            ! dr1/dx1 = -dr1/dx2
            tmp = pref*symtmp*ftot
            tmp1 = pref*angular*(ftot*dsymdr1 + dfcutdr1*symtmp)
            tmp2 = pref*angular*(ftot*dsymdr2 + dfcutdr2*symtmp)
            tmp3 = pref*angular*(ftot*dsymdr3 + dfcutdr3*symtmp)
            DO i = 1, 3
               force(i, 1, sf) = tmp*dangulardx1(i) + tmp1*dr1dx(i) + tmp2*dr2dx(i)
               force(i, 2, sf) = tmp*dangulardx2(i) - tmp1*dr1dx(i) + tmp3*dr3dx(i)
               force(i, 3, sf) = tmp*dangulardx3(i) - tmp2*dr2dx(i) - tmp3*dr3dx(i)
            END DO
         END IF

         ! Don't forget prefactor: 2**(1-ang%zeta)
         pref = prefzeta
         ! combine angular and exponential part with cutoff function
         sym(sf) = pref*angular*symtmp*ftot
      END DO

   END SUBROUTINE nnp_calc_ang

! **************************************************************************************************
!> \brief Sort element array according to atomic number
!> \param ele ...
!> \param nuc_ele ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_sort_ele(ele, nuc_ele)
      CHARACTER(len=2), DIMENSION(:), INTENT(INOUT)      :: ele
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: nuc_ele

      CHARACTER(len=2)                                   :: tmp_ele
      INTEGER                                            :: i, j, loc, minimum, tmp_nuc_ele

      ! Determine atomic number
      DO i = 1, SIZE(ele)
         CALL get_ptable_info(ele(i), number=nuc_ele(i))
      END DO

      ! Sort both arrays
      DO i = 1, SIZE(ele) - 1
         minimum = nuc_ele(i)
         loc = i
         DO j = i + 1, SIZE(ele)
            IF (nuc_ele(j) .LT. minimum) THEN
               loc = j
               minimum = nuc_ele(j)
            END IF
         END DO
         tmp_nuc_ele = nuc_ele(i)
         nuc_ele(i) = nuc_ele(loc)
         nuc_ele(loc) = tmp_nuc_ele

         tmp_ele = ele(i)
         ele(i) = ele(loc)
         ele(loc) = tmp_ele

      END DO

   END SUBROUTINE nnp_sort_ele

! **************************************************************************************************
!> \brief Sort symmetry functions according to different criteria
!> \param nnp ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_sort_acsf(nnp)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp

      INTEGER                                            :: i, j, k, loc

      ! First sort is according to symmetry function type
      ! This is done manually, since data structures are separate
      ! Note: Bubble sort is OK here, since small sort + special
      DO i = 1, nnp%n_ele
         ! sort by cutoff
         ! rad
         DO j = 1, nnp%n_rad(i) - 1
            loc = j
            DO k = j + 1, nnp%n_rad(i)
               IF (nnp%rad(i)%funccut(loc) .GT. nnp%rad(i)%funccut(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swaprad(nnp%rad(i), j, loc)
         END DO

         ! sort by eta
         ! rad
         DO j = 1, nnp%n_rad(i) - 1
            loc = j
            DO k = j + 1, nnp%n_rad(i)
               IF (nnp%rad(i)%funccut(loc) .EQ. nnp%rad(i)%funccut(k) .AND. &
                   nnp%rad(i)%eta(loc) .GT. nnp%rad(i)%eta(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swaprad(nnp%rad(i), j, loc)
         END DO

         ! sort by rshift
         ! rad
         DO j = 1, nnp%n_rad(i) - 1
            loc = j
            DO k = j + 1, nnp%n_rad(i)
               IF (nnp%rad(i)%funccut(loc) .EQ. nnp%rad(i)%funccut(k) .AND. &
                   nnp%rad(i)%eta(loc) .EQ. nnp%rad(i)%eta(k) .AND. &
                   nnp%rad(i)%rs(loc) .GT. nnp%rad(i)%rs(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swaprad(nnp%rad(i), j, loc)
         END DO

         ! sort by ele
         ! rad
         DO j = 1, nnp%n_rad(i) - 1
            loc = j
            DO k = j + 1, nnp%n_rad(i)
               IF (nnp%rad(i)%funccut(loc) .EQ. nnp%rad(i)%funccut(k) .AND. &
                   nnp%rad(i)%eta(loc) .EQ. nnp%rad(i)%eta(k) .AND. &
                   nnp%rad(i)%rs(loc) .EQ. nnp%rad(i)%rs(k) .AND. &
                   nnp%rad(i)%nuc_ele(loc) .GT. nnp%rad(i)%nuc_ele(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swaprad(nnp%rad(i), j, loc)
         END DO

         ! ang
         ! sort by cutoff
         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) .GT. nnp%ang(i)%funccut(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO

         ! sort by eta
         ! ang
         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) .EQ. nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) .GT. nnp%ang(i)%eta(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO

         ! sort by zeta
         ! ang
         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) .EQ. nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) .EQ. nnp%ang(i)%eta(k) .AND. &
                   nnp%ang(i)%zeta(loc) .GT. nnp%ang(i)%zeta(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO

         ! sort by lambda
         ! ang
         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) .EQ. nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) .EQ. nnp%ang(i)%eta(k) .AND. &
                   nnp%ang(i)%zeta(loc) .EQ. nnp%ang(i)%zeta(k) .AND. &
                   nnp%ang(i)%lam(loc) .GT. nnp%ang(i)%lam(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO

         ! sort by ele
         ! ang, ele1
         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) .EQ. nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) .EQ. nnp%ang(i)%eta(k) .AND. &
                   nnp%ang(i)%zeta(loc) .EQ. nnp%ang(i)%zeta(k) .AND. &
                   nnp%ang(i)%lam(loc) .EQ. nnp%ang(i)%lam(k) .AND. &
                   nnp%ang(i)%nuc_ele1(loc) .GT. nnp%ang(i)%nuc_ele1(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO
         ! ang, ele2
         DO j = 1, nnp%n_ang(i) - 1
            loc = j
            DO k = j + 1, nnp%n_ang(i)
               IF (nnp%ang(i)%funccut(loc) .EQ. nnp%ang(i)%funccut(k) .AND. &
                   nnp%ang(i)%eta(loc) .EQ. nnp%ang(i)%eta(k) .AND. &
                   nnp%ang(i)%zeta(loc) .EQ. nnp%ang(i)%zeta(k) .AND. &
                   nnp%ang(i)%lam(loc) .EQ. nnp%ang(i)%lam(k) .AND. &
                   nnp%ang(i)%nuc_ele1(loc) .EQ. nnp%ang(i)%nuc_ele1(k) .AND. &
                   nnp%ang(i)%nuc_ele2(loc) .GT. nnp%ang(i)%nuc_ele2(k)) THEN
                  loc = k
               END IF
            END DO
            ! swap symfnct
            CALL nnp_swapang(nnp%ang(i), j, loc)
         END DO
      END DO

   END SUBROUTINE nnp_sort_acsf

! **************************************************************************************************
!> \brief Swap two radial symmetry functions
!> \param rad ...
!> \param i ...
!> \param j ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_swaprad(rad, i, j)
      TYPE(nnp_acsf_rad_type), INTENT(INOUT)             :: rad
      INTEGER, INTENT(IN)                                :: i, j

      CHARACTER(len=2)                                   :: tmpc
      INTEGER                                            :: tmpi
      REAL(KIND=dp)                                      :: tmpr

      tmpr = rad%funccut(i)
      rad%funccut(i) = rad%funccut(j)
      rad%funccut(j) = tmpr

      tmpr = rad%eta(i)
      rad%eta(i) = rad%eta(j)
      rad%eta(j) = tmpr

      tmpr = rad%rs(i)
      rad%rs(i) = rad%rs(j)
      rad%rs(j) = tmpr

      tmpc = rad%ele(i)
      rad%ele(i) = rad%ele(j)
      rad%ele(j) = tmpc

      tmpi = rad%nuc_ele(i)
      rad%nuc_ele(i) = rad%nuc_ele(j)
      rad%nuc_ele(j) = tmpi

   END SUBROUTINE nnp_swaprad

! **************************************************************************************************
!> \brief Swap two angular symmetry functions
!> \param ang ...
!> \param i ...
!> \param j ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_swapang(ang, i, j)
      TYPE(nnp_acsf_ang_type), INTENT(INOUT)             :: ang
      INTEGER, INTENT(IN)                                :: i, j

      CHARACTER(len=2)                                   :: tmpc
      INTEGER                                            :: tmpi
      REAL(KIND=dp)                                      :: tmpr

      tmpr = ang%funccut(i)
      ang%funccut(i) = ang%funccut(j)
      ang%funccut(j) = tmpr

      tmpr = ang%eta(i)
      ang%eta(i) = ang%eta(j)
      ang%eta(j) = tmpr

      tmpr = ang%zeta(i)
      ang%zeta(i) = ang%zeta(j)
      ang%zeta(j) = tmpr

      tmpr = ang%prefzeta(i)
      ang%prefzeta(i) = ang%prefzeta(j)
      ang%prefzeta(j) = tmpr

      tmpr = ang%lam(i)
      ang%lam(i) = ang%lam(j)
      ang%lam(j) = tmpr

      tmpc = ang%ele1(i)
      ang%ele1(i) = ang%ele1(j)
      ang%ele1(j) = tmpc

      tmpi = ang%nuc_ele1(i)
      ang%nuc_ele1(i) = ang%nuc_ele1(j)
      ang%nuc_ele1(j) = tmpi

      tmpc = ang%ele2(i)
      ang%ele2(i) = ang%ele2(j)
      ang%ele2(j) = tmpc

      tmpi = ang%nuc_ele2(i)
      ang%nuc_ele2(i) = ang%nuc_ele2(j)
      ang%nuc_ele2(j) = tmpi

   END SUBROUTINE nnp_swapang

! **************************************************************************************************
!> \brief Initialize symmetry function groups
!> \param nnp ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_init_acsf_groups(nnp)

      TYPE(nnp_type), INTENT(INOUT)                      :: nnp

      INTEGER                                            :: ang, i, j, k, rad, s
      REAL(KIND=dp)                                      :: funccut

      !find out how many symmetry functions groups are needed
      DO i = 1, nnp%n_ele
         nnp%rad(i)%n_symfgrp = 0
         nnp%ang(i)%n_symfgrp = 0
         !search radial symmetry functions
         DO j = 1, nnp%n_ele
            funccut = -1.0_dp
            DO s = 1, nnp%n_rad(i)
               IF (nnp%rad(i)%ele(s) == nnp%ele(j)) THEN
                  IF (ABS(nnp%rad(i)%funccut(s) - funccut) > 1.0e-5_dp) THEN
                     nnp%rad(i)%n_symfgrp = nnp%rad(i)%n_symfgrp + 1
                     funccut = nnp%rad(i)%funccut(s)
                  END IF
               END IF
            END DO
         END DO
         !search angular symmetry functions
         DO j = 1, nnp%n_ele
            DO k = j, nnp%n_ele
               funccut = -1.0_dp
               DO s = 1, nnp%n_ang(i)
                  IF ((nnp%ang(i)%ele1(s) == nnp%ele(j) .AND. &
                       nnp%ang(i)%ele2(s) == nnp%ele(k)) .OR. &
                      (nnp%ang(i)%ele1(s) == nnp%ele(k) .AND. &
                       nnp%ang(i)%ele2(s) == nnp%ele(j))) THEN
                     IF (ABS(nnp%ang(i)%funccut(s) - funccut) > 1.0e-5_dp) THEN
                        nnp%ang(i)%n_symfgrp = nnp%ang(i)%n_symfgrp + 1
                        funccut = nnp%ang(i)%funccut(s)
                     END IF
                  END IF
               END DO
            END DO
         END DO
      END DO

      !allocate memory for symmetry functions group
      DO i = 1, nnp%n_ele
         ALLOCATE (nnp%rad(i)%symfgrp(nnp%rad(i)%n_symfgrp))
         ALLOCATE (nnp%ang(i)%symfgrp(nnp%ang(i)%n_symfgrp))
         DO j = 1, nnp%rad(i)%n_symfgrp
            nnp%rad(i)%symfgrp(j)%n_symf = 0
         END DO
         DO j = 1, nnp%ang(i)%n_symfgrp
            nnp%ang(i)%symfgrp(j)%n_symf = 0
         END DO
      END DO

      !init symmetry functions group
      DO i = 1, nnp%n_ele
         rad = 0
         ang = 0
         DO j = 1, nnp%n_ele
            funccut = -1.0_dp
            DO s = 1, nnp%n_rad(i)
               IF (nnp%rad(i)%ele(s) == nnp%ele(j)) THEN
                  IF (ABS(nnp%rad(i)%funccut(s) - funccut) > 1.0e-5_dp) THEN
                     rad = rad + 1
                     funccut = nnp%rad(i)%funccut(s)
                     nnp%rad(i)%symfgrp(rad)%cutoff = funccut
                     ALLOCATE (nnp%rad(i)%symfgrp(rad)%ele(1))
                     ALLOCATE (nnp%rad(i)%symfgrp(rad)%ele_ind(1))
                     nnp%rad(i)%symfgrp(rad)%ele(1) = nnp%ele(j)
                     nnp%rad(i)%symfgrp(rad)%ele_ind(1) = j
                  END IF
                  nnp%rad(i)%symfgrp(rad)%n_symf = nnp%rad(i)%symfgrp(rad)%n_symf + 1
               END IF
            END DO
         END DO
         DO j = 1, nnp%n_ele
            DO k = j, nnp%n_ele
               funccut = -1.0_dp
               DO s = 1, nnp%n_ang(i)
                  IF ((nnp%ang(i)%ele1(s) == nnp%ele(j) .AND. &
                       nnp%ang(i)%ele2(s) == nnp%ele(k)) .OR. &
                      (nnp%ang(i)%ele1(s) == nnp%ele(k) .AND. &
                       nnp%ang(i)%ele2(s) == nnp%ele(j))) THEN
                     IF (ABS(nnp%ang(i)%funccut(s) - funccut) > 1.0e-5_dp) THEN
                        ang = ang + 1
                        funccut = nnp%ang(i)%funccut(s)
                        nnp%ang(i)%symfgrp(ang)%cutoff = funccut
                        ALLOCATE (nnp%ang(i)%symfgrp(ang)%ele(2))
                        ALLOCATE (nnp%ang(i)%symfgrp(ang)%ele_ind(2))
                        nnp%ang(i)%symfgrp(ang)%ele(1) = nnp%ele(j)
                        nnp%ang(i)%symfgrp(ang)%ele(2) = nnp%ele(k)
                        nnp%ang(i)%symfgrp(ang)%ele_ind(1) = j
                        nnp%ang(i)%symfgrp(ang)%ele_ind(2) = k
                     END IF
                     nnp%ang(i)%symfgrp(ang)%n_symf = nnp%ang(i)%symfgrp(ang)%n_symf + 1
                  END IF
               END DO
            END DO
         END DO
      END DO

      !add symmetry functions to associated groups
      DO i = 1, nnp%n_ele
         DO j = 1, nnp%rad(i)%n_symfgrp
            ALLOCATE (nnp%rad(i)%symfgrp(j)%symf(nnp%rad(i)%symfgrp(j)%n_symf))
            rad = 0
            DO s = 1, nnp%n_rad(i)
               IF (nnp%rad(i)%ele(s) == nnp%rad(i)%symfgrp(j)%ele(1)) THEN
                  IF (ABS(nnp%rad(i)%funccut(s) - nnp%rad(i)%symfgrp(j)%cutoff) < 1.0e-7_dp) THEN
                     rad = rad + 1
                     nnp%rad(i)%symfgrp(j)%symf(rad) = s
                  END IF
               END IF
            END DO
         END DO
         DO j = 1, nnp%ang(i)%n_symfgrp
            ALLOCATE (nnp%ang(i)%symfgrp(j)%symf(nnp%ang(i)%symfgrp(j)%n_symf))
            ang = 0
            DO s = 1, nnp%n_ang(i)
               IF ((nnp%ang(i)%ele1(s) == nnp%ang(i)%symfgrp(j)%ele(1) .AND. &
                    nnp%ang(i)%ele2(s) == nnp%ang(i)%symfgrp(j)%ele(2)) .OR. &
                   (nnp%ang(i)%ele1(s) == nnp%ang(i)%symfgrp(j)%ele(2) .AND. &
                    nnp%ang(i)%ele2(s) == nnp%ang(i)%symfgrp(j)%ele(1))) THEN
                  IF (ABS(nnp%ang(i)%funccut(s) - nnp%ang(i)%symfgrp(j)%cutoff) < 1.0e-7_dp) THEN
                     ang = ang + 1
                     nnp%ang(i)%symfgrp(j)%symf(ang) = s
                  END IF
               END IF
            END DO
         END DO
      END DO

   END SUBROUTINE nnp_init_acsf_groups

! **************************************************************************************************
!> \brief Write symmetry function information
!> \param nnp ...
!> \param para_env ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_write_acsf(nnp, para_env)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=default_string_length)               :: my_label
      INTEGER                                            :: i, j, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      my_label = "NNP| "
      IF (para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger)
         WRITE (unit_nr, '(1X,A,1X,10(I2,1X))') TRIM(my_label)//" Activation functions:", nnp%actfnct(:)
         DO i = 1, nnp%n_ele
            WRITE (unit_nr, *) TRIM(my_label)//" short range atomic symmetry functions element "// &
               nnp%ele(i)//":"
            DO j = 1, nnp%n_rad(i)
               WRITE (unit_nr, '(1X,A,1X,I3,1X,A2,1X,I2,1X,A2,11X,3(F6.3,1X))') TRIM(my_label), j, nnp%ele(i), 2, &
                  nnp%rad(i)%ele(j), nnp%rad(i)%eta(j), &
                  nnp%rad(i)%rs(j), nnp%rad(i)%funccut(j)
            END DO
            DO j = 1, nnp%n_ang(i)
               WRITE (unit_nr, '(1X,A,1X,I3,1X,A2,1X,I2,2(1X,A2),1X,4(F6.3,1X))') &
                  TRIM(my_label), j, nnp%ele(i), 3, &
                  nnp%ang(i)%ele1(j), nnp%ang(i)%ele2(j), &
                  nnp%ang(i)%eta(j), nnp%ang(i)%lam(j), &
                  nnp%ang(i)%zeta(j), nnp%ang(i)%funccut(j)
            END DO
         END DO
      END IF

   END SUBROUTINE nnp_write_acsf

! **************************************************************************************************
!> \brief Create neighbor object
!> \param nnp ...
!> \param ind ...
!> \param nat ...
!> \param neighbor ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_neighbor_create(nnp, ind, nat, neighbor)

      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      INTEGER, INTENT(IN)                                :: ind, nat
      TYPE(nnp_neighbor_type), INTENT(INOUT)             :: neighbor

      INTEGER                                            :: n
      TYPE(cell_type), POINTER                           :: cell

      NULLIFY (cell)
      CALL nnp_env_get(nnp_env=nnp, cell=cell)

      CALL nnp_compute_pbc_copies(neighbor%pbc_copies, cell, nnp%max_cut)
      n = (SUM(neighbor%pbc_copies) + 1)*nat
      ALLOCATE (neighbor%dist_rad(4, n, nnp%rad(ind)%n_symfgrp))
      ALLOCATE (neighbor%dist_ang1(4, n, nnp%ang(ind)%n_symfgrp))
      ALLOCATE (neighbor%dist_ang2(4, n, nnp%ang(ind)%n_symfgrp))
      ALLOCATE (neighbor%ind_rad(n, nnp%rad(ind)%n_symfgrp))
      ALLOCATE (neighbor%ind_ang1(n, nnp%ang(ind)%n_symfgrp))
      ALLOCATE (neighbor%ind_ang2(n, nnp%ang(ind)%n_symfgrp))
      ALLOCATE (neighbor%n_rad(nnp%rad(ind)%n_symfgrp))
      ALLOCATE (neighbor%n_ang1(nnp%ang(ind)%n_symfgrp))
      ALLOCATE (neighbor%n_ang2(nnp%ang(ind)%n_symfgrp))
      neighbor%n_rad(:) = 0
      neighbor%n_ang1(:) = 0
      neighbor%n_ang2(:) = 0

   END SUBROUTINE nnp_neighbor_create

! **************************************************************************************************
!> \brief Destroy neighbor object
!> \param neighbor ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_neighbor_release(neighbor)

      TYPE(nnp_neighbor_type), INTENT(INOUT)             :: neighbor

      DEALLOCATE (neighbor%dist_rad)
      DEALLOCATE (neighbor%dist_ang1)
      DEALLOCATE (neighbor%dist_ang2)
      DEALLOCATE (neighbor%ind_rad)
      DEALLOCATE (neighbor%ind_ang1)
      DEALLOCATE (neighbor%ind_ang2)
      DEALLOCATE (neighbor%n_rad)
      DEALLOCATE (neighbor%n_ang1)
      DEALLOCATE (neighbor%n_ang2)

   END SUBROUTINE nnp_neighbor_release

! **************************************************************************************************
!> \brief Generate neighboring list for an atomic configuration
!> \param nnp ...
!> \param neighbor ...
!> \param i ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_compute_neighbors(nnp, neighbor, i)

      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      TYPE(nnp_neighbor_type), INTENT(INOUT)             :: neighbor
      INTEGER, INTENT(IN)                                :: i

      INTEGER                                            :: c1, c2, c3, ind, j, s
      INTEGER, DIMENSION(3)                              :: nl
      REAL(KIND=dp)                                      :: norm
      REAL(KIND=dp), DIMENSION(3)                        :: dr
      TYPE(cell_type), POINTER                           :: cell

      NULLIFY (cell)
      CALL nnp_env_get(nnp_env=nnp, cell=cell)

      ind = nnp%ele_ind(i)

      DO j = 1, nnp%num_atoms
         DO c1 = 1, 2*neighbor%pbc_copies(1) + 1
            nl(1) = -neighbor%pbc_copies(1) + c1 - 1
            DO c2 = 1, 2*neighbor%pbc_copies(2) + 1
               nl(2) = -neighbor%pbc_copies(2) + c2 - 1
               DO c3 = 1, 2*neighbor%pbc_copies(3) + 1
                  nl(3) = -neighbor%pbc_copies(3) + c3 - 1
                  IF (j == i .AND. nl(1) == 0 .AND. nl(2) == 0 .AND. nl(3) == 0) CYCLE
                  dr(:) = nnp%coord(:, i) - nnp%coord(:, j)
                  !Apply pbc, but subtract nl boxes from periodic image
                  dr = pbc(dr, cell, nl)
                  norm = NORM2(dr(:))
                  DO s = 1, nnp%rad(ind)%n_symfgrp
                     IF (nnp%ele_ind(j) == nnp%rad(ind)%symfgrp(s)%ele_ind(1)) THEN
                        IF (norm < nnp%rad(ind)%symfgrp(s)%cutoff) THEN
                           neighbor%n_rad(s) = neighbor%n_rad(s) + 1
                           neighbor%ind_rad(neighbor%n_rad(s), s) = j
                           neighbor%dist_rad(1:3, neighbor%n_rad(s), s) = dr(:)
                           neighbor%dist_rad(4, neighbor%n_rad(s), s) = norm
                        END IF
                     END IF
                  END DO
                  DO s = 1, nnp%ang(ind)%n_symfgrp
                     IF (norm < nnp%ang(ind)%symfgrp(s)%cutoff) THEN
                        IF (nnp%ele_ind(j) == nnp%ang(ind)%symfgrp(s)%ele_ind(1)) THEN
                           neighbor%n_ang1(s) = neighbor%n_ang1(s) + 1
                           neighbor%ind_ang1(neighbor%n_ang1(s), s) = j
                           neighbor%dist_ang1(1:3, neighbor%n_ang1(s), s) = dr(:)
                           neighbor%dist_ang1(4, neighbor%n_ang1(s), s) = norm
                        END IF
                        IF (nnp%ele_ind(j) == nnp%ang(ind)%symfgrp(s)%ele_ind(2)) THEN
                           neighbor%n_ang2(s) = neighbor%n_ang2(s) + 1
                           neighbor%ind_ang2(neighbor%n_ang2(s), s) = j
                           neighbor%dist_ang2(1:3, neighbor%n_ang2(s), s) = dr(:)
                           neighbor%dist_ang2(4, neighbor%n_ang2(s), s) = norm
                        END IF
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE nnp_compute_neighbors

! **************************************************************************************************
!> \brief Determine required pbc copies for small cells
!> \param pbc_copies ...
!> \param cell ...
!> \param cutoff ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_compute_pbc_copies(pbc_copies, cell, cutoff)
      INTEGER, DIMENSION(3), INTENT(INOUT)               :: pbc_copies
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      REAL(KIND=dp), INTENT(IN)                          :: cutoff

      REAL(KIND=dp)                                      :: proja, projb, projc
      REAL(KIND=dp), DIMENSION(3)                        :: axb, axc, bxc

      axb(1) = cell%hmat(2, 1)*cell%hmat(3, 2) - cell%hmat(3, 1)*cell%hmat(2, 2)
      axb(2) = cell%hmat(3, 1)*cell%hmat(1, 2) - cell%hmat(1, 1)*cell%hmat(3, 2)
      axb(3) = cell%hmat(1, 1)*cell%hmat(2, 2) - cell%hmat(2, 1)*cell%hmat(1, 2)
      axb(:) = axb(:)/NORM2(axb(:))

      axc(1) = cell%hmat(2, 1)*cell%hmat(3, 3) - cell%hmat(3, 1)*cell%hmat(2, 3)
      axc(2) = cell%hmat(3, 1)*cell%hmat(1, 3) - cell%hmat(1, 1)*cell%hmat(3, 3)
      axc(3) = cell%hmat(1, 1)*cell%hmat(2, 3) - cell%hmat(2, 1)*cell%hmat(1, 3)
      axc(:) = axc(:)/NORM2(axc(:))

      bxc(1) = cell%hmat(2, 2)*cell%hmat(3, 3) - cell%hmat(3, 2)*cell%hmat(2, 3)
      bxc(2) = cell%hmat(3, 2)*cell%hmat(1, 3) - cell%hmat(1, 2)*cell%hmat(3, 3)
      bxc(3) = cell%hmat(1, 2)*cell%hmat(2, 3) - cell%hmat(2, 2)*cell%hmat(1, 3)
      bxc(:) = bxc(:)/NORM2(bxc(:))

      proja = ABS(SUM(cell%hmat(:, 1)*bxc(:)))*0.5_dp
      projb = ABS(SUM(cell%hmat(:, 2)*axc(:)))*0.5_dp
      projc = ABS(SUM(cell%hmat(:, 3)*axb(:)))*0.5_dp

      pbc_copies(:) = 0
      DO WHILE ((pbc_copies(1) + 1)*proja <= cutoff)
         pbc_copies(1) = pbc_copies(1) + 1
      END DO
      DO WHILE ((pbc_copies(2) + 1)*projb <= cutoff)
         pbc_copies(2) = pbc_copies(2) + 1
      END DO
      DO WHILE ((pbc_copies(3) + 1)*projc <= cutoff)
         pbc_copies(3) = pbc_copies(3) + 1
      END DO
      ! Apply non periodic setting
      pbc_copies(:) = pbc_copies(:)*cell%perd(:)

   END SUBROUTINE nnp_compute_pbc_copies

END MODULE nnp_acsf
