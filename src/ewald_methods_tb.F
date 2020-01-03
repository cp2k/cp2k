!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of Ewald contributions in DFTB
!> \author JGH
! **************************************************************************************************
MODULE ewald_methods_tb
   USE atprop_types,                    ONLY: atprop_type
   USE cell_types,                      ONLY: cell_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dgs,                             ONLY: dg_sum_patch,&
                                              dg_sum_patch_force_1d,&
                                              dg_sum_patch_force_3d
   USE ewald_environment_types,         ONLY: ewald_env_get,&
                                              ewald_environment_type
   USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                              ewald_pw_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi,&
                                              oorootpi
   USE particle_types,                  ONLY: particle_type
   USE pme_tools,                       ONLY: get_center,&
                                              set_list
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_grids,                        ONLY: get_pw_grid_info
   USE pw_methods,                      ONLY: pw_copy,&
                                              pw_derive,&
                                              pw_integral_a2b,&
                                              pw_transfer
   USE pw_poisson_methods,              ONLY: pw_poisson_rebuild,&
                                              pw_poisson_solve
   USE pw_poisson_types,                ONLY: greens_fn_type,&
                                              pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE realspace_grid_types,            ONLY: &
        pw2rs, realspace_grid_desc_type, realspace_grid_p_type, realspace_grid_type, rs2pw, &
        rs_grid_create, rs_grid_release, rs_grid_set_box, rs_grid_zero, rs_pw_transfer
   USE spme,                            ONLY: get_patch
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ewald_methods_tb'

   PUBLIC :: tb_spme_evaluate, tb_ewald_overlap

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param particle_set ...
!> \param box ...
!> \param gmcharge ...
!> \param mcharge ...
!> \param calculate_forces ...
!> \param virial ...
!> \param use_virial ...
!> \param atprop ...
! **************************************************************************************************
   SUBROUTINE tb_spme_evaluate(ewald_env, ewald_pw, particle_set, box, &
                               gmcharge, mcharge, calculate_forces, virial, use_virial, atprop)

      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      TYPE(cell_type), POINTER                           :: box
      REAL(KIND=dp), DIMENSION(:, :), INTENT(inout)      :: gmcharge
      REAL(KIND=dp), DIMENSION(:), INTENT(in)            :: mcharge
      LOGICAL, INTENT(in)                                :: calculate_forces
      TYPE(virial_type), POINTER                         :: virial
      LOGICAL, INTENT(in)                                :: use_virial
      TYPE(atprop_type), POINTER                         :: atprop

      CHARACTER(len=*), PARAMETER :: routineN = 'tb_spme_evaluate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: group, handle, i, ig, ipart, j, n, &
                                                            npart, o_spline, p1
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: center
      INTEGER, DIMENSION(3)                              :: nd, npts
      REAL(KIND=dp)                                      :: alpha, dvols, fat(3), ffa, ffb, fint, vgc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: rhos
      REAL(KIND=dp), DIMENSION(3, 3)                     :: f_stress, h_stress
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(greens_fn_type), POINTER                      :: green
      TYPE(pw_grid_type), POINTER                        :: grid_spme
      TYPE(pw_p_type), DIMENSION(3)                      :: dphi_g
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(pw_type), POINTER                             :: phi_g, phib_g, rhob_g, rhob_r
      TYPE(realspace_grid_desc_type), POINTER            :: rs_desc
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: drpot
      TYPE(realspace_grid_type), POINTER                 :: rden, rpot

      CALL timeset(routineN, handle)
      !-------------- INITIALISATION ---------------------
      CALL ewald_env_get(ewald_env, alpha=alpha, o_spline=o_spline, group=group, &
                         para_env=para_env)
      NULLIFY (green, poisson_env, pw_pool)
      CALL ewald_pw_get(ewald_pw, pw_big_pool=pw_pool, rs_desc=rs_desc, &
                        poisson_env=poisson_env)
      CALL pw_poisson_rebuild(poisson_env)
      green => poisson_env%green_fft
      grid_spme => pw_pool%pw_grid

      CALL get_pw_grid_info(grid_spme, dvol=dvols, npts=npts)

      npart = SIZE(particle_set)

      n = o_spline
      ALLOCATE (rhos(n, n, n))

      CALL rs_grid_create(rden, rs_desc)
      CALL rs_grid_set_box(grid_spme, rs=rden)
      CALL rs_grid_zero(rden)

      ALLOCATE (center(3, npart))
      CALL get_center(particle_set, box, center, npts, n)

      !-------------- DENSITY CALCULATION ----------------
      ipart = 0
      DO
         CALL set_list(particle_set, npart, center, p1, rden, ipart)
         IF (p1 == 0) EXIT

         ! calculate function on small boxes
         CALL get_patch(particle_set, box, green, npts, p1, rhos, is_core=.FALSE., &
                        is_shell=.FALSE., unit_charge=.TRUE.)
         rhos(:, :, :) = rhos(:, :, :)*mcharge(p1)

         ! add boxes to real space grid (big box)
         CALL dg_sum_patch(rden, rhos, center(:, p1))
      END DO

      CALL pw_pool_create_pw(pw_pool, rhob_r, use_data=REALDATA3D, &
                             in_space=REALSPACE)

      CALL rs_pw_transfer(rden, rhob_r, rs2pw)

      ! transform density to G space and add charge function
      CALL pw_pool_create_pw(pw_pool, rhob_g, use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_transfer(rhob_r, rhob_g)
      ! update charge function
      rhob_g%cc = rhob_g%cc*green%p3m_charge%cr

      !-------------- ELECTROSTATIC CALCULATION -----------

      ! allocate intermediate arrays
      DO i = 1, 3
         NULLIFY (dphi_g(i)%pw)
         CALL pw_pool_create_pw(pw_pool, dphi_g(i)%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      END DO
      CALL pw_pool_create_pw(pw_pool, phi_g, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      IF (use_virial) THEN
         CALL pw_poisson_solve(poisson_env, rhob_g, vgc, phi_g, dphi_g, h_stress=h_stress)
      ELSE
         CALL pw_poisson_solve(poisson_env, rhob_g, vgc, phi_g, dphi_g)
      END IF

      CALL rs_grid_create(rpot, rs_desc)
      CALL rs_grid_set_box(grid_spme, rs=rpot)

      ! Atomic Stress
      IF (ASSOCIATED(atprop)) THEN
         IF (atprop%stress .AND. use_virial) THEN

            CALL rs_grid_zero(rpot)
            rhob_g%cc = phi_g%cc*green%p3m_charge%cr
            CALL pw_transfer(rhob_g, rhob_r)
            CALL rs_pw_transfer(rpot, rhob_r, pw2rs)
            ipart = 0
            DO
               CALL set_list(particle_set, npart, center, p1, rden, ipart)
               IF (p1 == 0) EXIT
               ! calculate function on small boxes
               CALL get_patch(particle_set, box, green, npts, p1, rhos, is_core=.FALSE., &
                              is_shell=.FALSE., unit_charge=.TRUE.)
               CALL dg_sum_patch_force_1d(rpot, rhos, center(:, p1), fint)
               atprop%atstress(1, 1, p1) = atprop%atstress(1, 1, p1) + 0.5_dp*mcharge(p1)*fint*dvols
               atprop%atstress(2, 2, p1) = atprop%atstress(2, 2, p1) + 0.5_dp*mcharge(p1)*fint*dvols
               atprop%atstress(3, 3, p1) = atprop%atstress(3, 3, p1) + 0.5_dp*mcharge(p1)*fint*dvols
            END DO

            CALL pw_pool_create_pw(pw_pool, phib_g, &
                                   use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            ffa = (0.5_dp/alpha)**2
            ffb = 1.0_dp/fourpi
            DO i = 1, 3
               DO ig = grid_spme%first_gne0, grid_spme%ngpts_cut_local
                  phib_g%cc(ig) = ffb*dphi_g(i)%pw%cc(ig)*(1.0_dp + ffa*grid_spme%gsq(ig))
                  phib_g%cc(ig) = phib_g%cc(ig)*green%influence_fn%cc(ig)
               END DO
               IF (grid_spme%have_g0) phib_g%cc(1) = 0.0_dp
               DO j = 1, i
                  nd = 0
                  nd(j) = 1
                  CALL pw_copy(phib_g, rhob_g)
                  CALL pw_derive(rhob_g, nd)
                  rhob_g%cc = rhob_g%cc*green%p3m_charge%cr
                  CALL pw_transfer(rhob_g, rhob_r)
                  CALL rs_pw_transfer(rpot, rhob_r, pw2rs)

                  ipart = 0
                  DO
                     CALL set_list(particle_set, npart, center, p1, rden, ipart)
                     IF (p1 == 0) EXIT
                     ! calculate function on small boxes
                     CALL get_patch(particle_set, box, green, npts, p1, rhos, &
                                    is_core=.FALSE., is_shell=.FALSE., unit_charge=.TRUE.)
                     ! integrate box and potential
                     CALL dg_sum_patch_force_1d(rpot, rhos, center(:, p1), fint)
                     atprop%atstress(i, j, p1) = atprop%atstress(i, j, p1) + fint*dvols*mcharge(p1)
                     IF (i /= j) atprop%atstress(j, i, p1) = atprop%atstress(j, i, p1) + fint*dvols*mcharge(p1)
                  END DO

               END DO
            END DO
            CALL pw_pool_give_back_pw(pw_pool, phib_g)

         END IF
      END IF

      CALL pw_pool_give_back_pw(pw_pool, rhob_g)

      CALL rs_grid_zero(rpot)
      phi_g%cc = phi_g%cc*green%p3m_charge%cr
      CALL pw_transfer(phi_g, rhob_r)
      CALL pw_pool_give_back_pw(pw_pool, phi_g)
      CALL rs_pw_transfer(rpot, rhob_r, pw2rs)

      !---------- END OF ELECTROSTATIC CALCULATION --------

      !------------- STRESS TENSOR CALCULATION ------------

      IF (use_virial) THEN
         DO i = 1, 3
            DO j = i, 3
               f_stress(i, j) = pw_integral_a2b(dphi_g(i)%pw, dphi_g(j)%pw)
               f_stress(j, i) = f_stress(i, j)
            END DO
         END DO
         ffa = (1.0_dp/fourpi)*(0.5_dp/alpha)**2
         virial%pv_virial = virial%pv_virial - (ffa*f_stress - h_stress)/REAL(para_env%num_pe, dp)
      END IF

      !--------END OF STRESS TENSOR CALCULATION -----------

      IF (calculate_forces) THEN
         ! move derivative of potential to real space grid and
         ! multiply by charge function in g-space
         ALLOCATE (drpot(1:3))
         DO i = 1, 3
            CALL rs_grid_create(drpot(i)%rs_grid, rs_desc)
            CALL rs_grid_set_box(grid_spme, rs=drpot(i)%rs_grid)
            dphi_g(i)%pw%cc = dphi_g(i)%pw%cc*green%p3m_charge%cr
            CALL pw_transfer(dphi_g(i)%pw, rhob_r)
            CALL pw_pool_give_back_pw(pw_pool, dphi_g(i)%pw)
            CALL rs_pw_transfer(drpot(i)%rs_grid, rhob_r, pw2rs)
         END DO
      ELSE
         DO i = 1, 3
            CALL pw_pool_give_back_pw(pw_pool, dphi_g(i)%pw)
         END DO
      END IF
      CALL pw_pool_give_back_pw(pw_pool, rhob_r)

      !----------------- FORCE CALCULATION ----------------

      ipart = 0
      DO

         CALL set_list(particle_set, npart, center, p1, rden, ipart)
         IF (p1 == 0) EXIT

         ! calculate function on small boxes
         CALL get_patch(particle_set, box, green, npts, p1, rhos, is_core=.FALSE., &
                        is_shell=.FALSE., unit_charge=.TRUE.)

         CALL dg_sum_patch_force_1d(rpot, rhos, center(:, p1), fint)
         gmcharge(p1, 1) = gmcharge(p1, 1) + fint*dvols

         IF (calculate_forces) THEN
            CALL dg_sum_patch_force_3d(drpot, rhos, center(:, p1), fat)
            gmcharge(p1, 2) = gmcharge(p1, 2) - fat(1)*dvols
            gmcharge(p1, 3) = gmcharge(p1, 3) - fat(2)*dvols
            gmcharge(p1, 4) = gmcharge(p1, 4) - fat(3)*dvols
         END IF

      END DO

      !--------------END OF FORCE CALCULATION -------------

      !------------------CLEANING UP ----------------------

      CALL rs_grid_release(rden)
      CALL rs_grid_release(rpot)
      IF (calculate_forces) THEN
         IF (ASSOCIATED(drpot)) THEN
            DO i = 1, 3
               CALL rs_grid_release(drpot(i)%rs_grid)
            END DO
            DEALLOCATE (drpot)
         END IF
      END IF
      DEALLOCATE (rhos)
      DEALLOCATE (center)

      CALL timestop(handle)

   END SUBROUTINE tb_spme_evaluate

! **************************************************************************************************
!> \brief ...
!> \param gmcharge ...
!> \param mcharge ...
!> \param alpha ...
!> \param n_list ...
!> \param virial ...
!> \param use_virial ...
!> \param atprop ...
! **************************************************************************************************
   SUBROUTINE tb_ewald_overlap(gmcharge, mcharge, alpha, n_list, virial, use_virial, atprop)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(inout)      :: gmcharge
      REAL(KIND=dp), DIMENSION(:), INTENT(in)            :: mcharge
      REAL(KIND=dp), INTENT(in)                          :: alpha
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list
      TYPE(virial_type), POINTER                         :: virial
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(atprop_type), POINTER                         :: atprop

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tb_ewald_overlap', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iatom, jatom, nmat
      REAL(KIND=dp)                                      :: dfr, dr, fr, pfr, rij(3)
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      CALL timeset(routineN, handle)

      nmat = SIZE(gmcharge, 2)

      CALL neighbor_list_iterator_create(nl_iterator, n_list)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, r=rij)

         dr = SQRT(SUM(rij(:)**2))
         IF (dr > 1.e-10) THEN
            fr = erfc(alpha*dr)/dr
            gmcharge(iatom, 1) = gmcharge(iatom, 1) + mcharge(jatom)*fr
            gmcharge(jatom, 1) = gmcharge(jatom, 1) + mcharge(iatom)*fr
            IF (nmat > 1) THEN
               dfr = -2._dp*alpha*EXP(-alpha*alpha*dr*dr)*oorootpi/dr - fr/dr
               dfr = -dfr/dr
               DO i = 2, nmat
                  gmcharge(iatom, i) = gmcharge(iatom, i) - rij(i - 1)*mcharge(jatom)*dfr
                  gmcharge(jatom, i) = gmcharge(jatom, i) + rij(i - 1)*mcharge(iatom)*dfr
               END DO
            END IF
            IF (use_virial) THEN
               IF (iatom == jatom) THEN
                  pfr = -0.5_dp*dfr*mcharge(iatom)*mcharge(jatom)
               ELSE
                  pfr = -dfr*mcharge(iatom)*mcharge(jatom)
               END IF
               CALL virial_pair_force(virial%pv_virial, -pfr, rij, rij)
               IF (ASSOCIATED(atprop)) THEN
                  IF (atprop%stress) THEN
                     CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp*pfr, rij, rij)
                     CALL virial_pair_force(atprop%atstress(:, :, jatom), -0.5_dp*pfr, rij, rij)
                  END IF
               END IF
            END IF
         END IF

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      CALL timestop(handle)

   END SUBROUTINE tb_ewald_overlap

END MODULE ewald_methods_tb

