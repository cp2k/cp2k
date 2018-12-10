!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of Coulomb contributions in xTB
!> \author JGH
! **************************************************************************************************
MODULE xtb_coulomb
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_array_init,&
                                              atprop_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_get_block_p,&
                                              dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop,&
                                              dbcsr_iterator_type,&
                                              dbcsr_p_type
   USE dgs,                             ONLY: dg_sum_patch,&
                                              dg_sum_patch_force_1d,&
                                              dg_sum_patch_force_3d
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE erf_fn,                          ONLY: erfc
   USE ewald_environment_types,         ONLY: ewald_env_get,&
                                              ewald_environment_type
   USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                              ewald_pw_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE mathconstants,                   ONLY: fourpi,&
                                              oorootpi,&
                                              pi
   USE message_passing,                 ONLY: mp_sum
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
   USE pw_poisson_types,                ONLY: do_ewald_ewald,&
                                              do_ewald_none,&
                                              do_ewald_pme,&
                                              do_ewald_spme,&
                                              greens_fn_type,&
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
   USE qs_dftb3_methods,                ONLY: build_dftb3_diagonal
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE realspace_grid_types,            ONLY: &
        pw2rs, realspace_grid_desc_type, realspace_grid_p_type, realspace_grid_type, rs2pw, &
        rs_grid_create, rs_grid_release, rs_grid_set_box, rs_grid_zero, rs_pw_transfer
   USE spme,                            ONLY: get_patch
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_coulomb'

   PUBLIC :: build_xtb_coulomb

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ks_matrix ...
!> \param rho ...
!> \param charges ...
!> \param mcharge ...
!> \param energy ...
!> \param calculate_forces ...
!> \param just_energy ...
! **************************************************************************************************
   SUBROUTINE build_xtb_coulomb(qs_env, ks_matrix, rho, charges, mcharge, energy, &
                                calculate_forces, just_energy)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(qs_rho_type), POINTER                         :: rho
      REAL(dp), DIMENSION(:, :)                          :: charges
      REAL(dp), DIMENSION(:)                             :: mcharge
      TYPE(qs_energy_type), POINTER                      :: energy
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy

      CHARACTER(len=*), PARAMETER :: routineN = 'build_xtb_coulomb', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_i, atom_j, blk, ewald_type, handle, &
                                                            i, ia, iatom, ic, icol, ikind, img, &
                                                            irow, is, j, jatom, jkind, natom, ni, &
                                                            nimg, nj, nkind, nmat
      INTEGER, DIMENSION(3)                              :: cellind, periodic
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind, kind_of
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_ewald, found, use_virial
      REAL(KIND=dp)                                      :: alpha, deth, dr, ecsr, fi, gmij, zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: xgamma, zeffk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gcij, gmcharge
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: gchrg
      REAL(KIND=dp), DIMENSION(25)                       :: gcint
      REAL(KIND=dp), DIMENSION(3)                        :: fij, rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dsblock, gblock, ksblock, pblock, sblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gamma_matrix
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)

      NULLIFY (gamma_matrix, matrix_p, matrix_s, virial, atprop, dft_control)

      CALL get_qs_env(qs_env, &
                      gamma_matrix=gamma_matrix, &
                      particle_set=particle_set, &
                      cell=cell, &
                      virial=virial, &
                      atprop=atprop, &
                      dft_control=dft_control)

      use_virial = .FALSE.
      IF (calculate_forces) THEN
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      END IF

      IF (atprop%energy) THEN
         CALL get_qs_env(qs_env=qs_env, particle_set=particle_set)
         natom = SIZE(particle_set)
         CALL atprop_array_init(atprop%atecoul, natom)
      END IF

      IF (calculate_forces) THEN
         nmat = 4
      ELSE
         nmat = 1
      END IF

      natom = SIZE(mcharge)
      ALLOCATE (gchrg(natom, 5, nmat))
      gchrg = 0._dp
      ALLOCATE (gmcharge(natom, nmat))
      gmcharge = 0._dp

      ! short range contribution (gamma)
      CALL dbcsr_iterator_start(iter, gamma_matrix(1)%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         NULLIFY (gblock)
         CALL dbcsr_iterator_next_block(iter, iatom, jatom, gblock, blk)
         ni = SIZE(gblock, 1)
         nj = SIZE(gblock, 2)
         gcint(1:ni) = MATMUL(gblock, charges(jatom, 1:nj))
         gchrg(iatom, 1:ni, 1) = gchrg(iatom, 1:ni, 1)+gcint(1:ni)
         IF (calculate_forces) THEN
            DO i = 2, 4
               NULLIFY (gblock)
               CALL dbcsr_get_block_p(matrix=gamma_matrix(i)%matrix, &
                                      row=iatom, col=jatom, block=gblock, found=found)
               ni = SIZE(gblock, 1)
               nj = SIZE(gblock, 2)
               gcint(1:ni) = MATMUL(gblock, charges(jatom, 1:nj))
               gchrg(iatom, 1:ni, i) = gchrg(iatom, 1:ni, i)+gcint(1:ni)
            END DO
         END IF
      ENDDO
      CALL dbcsr_iterator_stop(iter)

      energy%hartree = energy%hartree+ecsr

      IF (calculate_forces .AND. use_virial) THEN
         CALL dbcsr_iterator_start(iter, gamma_matrix(1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            NULLIFY (gblock)
            CALL dbcsr_iterator_next_block(iter, iatom, jatom, gblock, blk)
            DO i = 1, 3
               NULLIFY (gblock)
               CALL dbcsr_get_block_p(matrix=gamma_matrix(i+1)%matrix, &
                                      row=iatom, col=jatom, block=gblock, found=found)
               ni = SIZE(gblock, 1)
               nj = SIZE(gblock, 2)
               gcint(1:ni) = MATMUL(gblock, charges(jatom, 1:nj))
               fij(i) = SUM(charges(iatom, 1:ni)*gcint(1:ni))
            END DO
            rij = particle_set(iatom)%r-particle_set(jatom)%r
            rij = pbc(rij, cell)
            CALL virial_pair_force(virial%pv_virial, 1._dp, fij, rij)
            IF (atprop%stress) THEN
               CALL virial_pair_force(atprop%atstress(:, :, iatom), 0.5_dp, fij, rij)
               CALL virial_pair_force(atprop%atstress(:, :, jatom), 0.5_dp, fij, rij)
            END IF
         ENDDO
         CALL dbcsr_iterator_stop(iter)
      END IF

      ! 1/R contribution

      do_ewald = dft_control%qs_control%xtb_control%do_ewald
      IF (do_ewald) THEN
         ! Ewald sum
         NULLIFY (ewald_env, ewald_pw)
         CALL get_qs_env(qs_env=qs_env, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw)
         CALL get_cell(cell=cell, periodic=periodic, deth=deth)
         CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type)
         CALL get_qs_env(qs_env=qs_env, sab_tbe=n_list)
         CALL ewald_overlap(gmcharge, mcharge, alpha, n_list, virial, use_virial, atprop)
         SELECT CASE (ewald_type)
         CASE DEFAULT
            CPABORT("Invalid Ewald type")
         CASE (do_ewald_none)
            CPABORT("Not allowed with DFTB")
         CASE (do_ewald_ewald)
            CPABORT("Standard Ewald not implemented in DFTB")
         CASE (do_ewald_pme)
            CPABORT("PME not implemented in DFTB")
         CASE (do_ewald_spme)
            CALL spme_evaluate(ewald_env, ewald_pw, particle_set, cell, &
                               gmcharge, mcharge, calculate_forces, virial, use_virial, atprop)
         END SELECT
      ELSE
         ! direct sum
         CALL get_qs_env(qs_env=qs_env, &
                         local_particles=local_particles)
         DO ikind = 1, SIZE(local_particles%n_el)
            DO ia = 1, local_particles%n_el(ikind)
               iatom = local_particles%list(ikind)%array(ia)
               DO jatom = 1, iatom-1
                  rij = particle_set(iatom)%r-particle_set(jatom)%r
                  rij = pbc(rij, cell)
                  dr = SQRT(SUM(rij(:)**2))
                  gmcharge(iatom, 1) = gmcharge(iatom, 1)+mcharge(jatom)/dr
                  gmcharge(jatom, 1) = gmcharge(jatom, 1)+mcharge(iatom)/dr
                  DO i = 2, nmat
                     gmcharge(iatom, i) = gmcharge(iatom, i)+rij(i-1)*mcharge(jatom)/dr**3
                     gmcharge(jatom, i) = gmcharge(jatom, i)-rij(i-1)*mcharge(iatom)/dr**3
                  END DO
               END DO
            END DO
         END DO
         CPASSERT(.NOT. use_virial)
      END IF

      ! global sum of gamma*p arrays
      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      force=force, para_env=para_env)
      CALL mp_sum(gmcharge(:, 1), para_env%group)
      CALL mp_sum(gchrg(:, :, 1), para_env%group)

      IF (do_ewald) THEN
         ! add self charge interaction and background charge contribution
         gmcharge(:, 1) = gmcharge(:, 1)-2._dp*alpha*oorootpi*mcharge(:)
         IF (ANY(periodic(:) == 1)) THEN
            gmcharge(:, 1) = gmcharge(:, 1)-pi/alpha**2/deth
         END IF
      END IF

      ! energy
      ALLOCATE (atom_of_kind(natom), kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               kind_of=kind_of, &
                               atom_of_kind=atom_of_kind)
      ecsr = 0.0_dp
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
         CALL get_xtb_atom_param(xtb_kind, nshell=ni)
         ecsr = ecsr+SUM(charges(iatom, 1:ni)*gchrg(iatom, 1:ni, 1))
      END DO
      energy%hartree = energy%hartree+0.5_dp*ecsr
      energy%hartree = energy%hartree+0.5_dp*SUM(mcharge(:)*gmcharge(:, 1))
      IF (atprop%energy) THEN
         CALL get_qs_env(qs_env=qs_env, local_particles=local_particles)
         DO ikind = 1, SIZE(local_particles%n_el)
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, nshell=ni)
            DO ia = 1, local_particles%n_el(ikind)
               iatom = local_particles%list(ikind)%array(ia)
               atprop%atecoul(iatom) = atprop%atecoul(iatom)+ &
                                       0.5_dp*SUM(charges(iatom, 1:ni)*gchrg(iatom, 1:ni, 1))
               atprop%atecoul(iatom) = atprop%atecoul(iatom)+ &
                                       0.5_dp*mcharge(iatom)*gmcharge(iatom, 1)
            END DO
         END DO
      END IF

      IF (calculate_forces) THEN
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atom_i = atom_of_kind(iatom)
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, nshell=ni)
            DO i = 1, 3
               fij(i) = SUM(charges(iatom, 1:ni)*gchrg(iatom, 1:ni, i+1))+ &
                        gmcharge(iatom, i+1)*mcharge(iatom)
            END DO
            force(ikind)%rho_elec(1, atom_i) = force(ikind)%rho_elec(1, atom_i)-fij(1)
            force(ikind)%rho_elec(2, atom_i) = force(ikind)%rho_elec(2, atom_i)-fij(2)
            force(ikind)%rho_elec(3, atom_i) = force(ikind)%rho_elec(3, atom_i)-fij(3)
         END DO
      END IF

      IF (.NOT. just_energy) THEN
         CALL get_qs_env(qs_env=qs_env, matrix_s_kp=matrix_s)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)

         nimg = dft_control%nimages
         NULLIFY (cell_to_index)
         IF (nimg > 1) THEN
            NULLIFY (kpoints)
            CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
            CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
         END IF

         IF (calculate_forces .AND. SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimg
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            END DO
         END IF

         NULLIFY (n_list)
         CALL get_qs_env(qs_env=qs_env, sab_orb=n_list)
         CALL neighbor_list_iterator_create(nl_iterator, n_list)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                   iatom=iatom, jatom=jatom, r=rij, cell=cellind)

            icol = MAX(iatom, jatom)
            irow = MIN(iatom, jatom)

            IF (nimg == 1) THEN
               ic = 1
            ELSE
               ic = cell_to_index(cellind(1), cellind(2), cellind(3))
               CPASSERT(ic > 0)
            END IF

            NULLIFY (sblock)
            CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                   row=irow, col=icol, block=sblock, found=found)
            CPASSERT(found)
            ni = SIZE(sblock, 1)
            nj = SIZE(sblock, 2)
            ALLOCATE (gcij(ni, nj))
            DO i = 1, ni
               DO j = 1, nj
                  gcij(i, j) = 0.5_dp*(gchrg(iatom, i, 1)+gchrg(jatom, j, 1))
               END DO
            END DO
            gmij = 0.5_dp*(gmcharge(iatom, 1)+gmcharge(jatom, 1))
            DO is = 1, SIZE(ks_matrix, 1)
               NULLIFY (ksblock)
               CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                      row=irow, col=icol, block=ksblock, found=found)
               CPASSERT(found)
               ksblock = ksblock-gcij*sblock-gmij*sblock
            END DO

            IF (calculate_forces) THEN
               ikind = kind_of(iatom)
               atom_i = atom_of_kind(iatom)
               jkind = kind_of(jatom)
               atom_j = atom_of_kind(jatom)
               IF (irow == jatom) gmij = -gmij
               NULLIFY (pblock)
               CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                      row=irow, col=icol, block=pblock, found=found)
               CPASSERT(found)
               DO i = 1, 3
                  NULLIFY (dsblock)
                  CALL dbcsr_get_block_p(matrix=matrix_s(1+i, ic)%matrix, &
                                         row=irow, col=icol, block=dsblock, found=found)
                  CPASSERT(found)
                  fi = -2.0_dp*(SUM(pblock*dsblock*gcij)+gmij*SUM(pblock*dsblock))
                  force(ikind)%rho_elec(i, atom_i) = force(ikind)%rho_elec(i, atom_i)+fi
                  force(jkind)%rho_elec(i, atom_j) = force(jkind)%rho_elec(i, atom_j)-fi
                  fij(i) = fi
               END DO
               IF (use_virial) THEN
                  CALL virial_pair_force(virial%pv_virial, 1._dp, fij, rij)
                  IF (atprop%stress) THEN
                     CALL virial_pair_force(atprop%atstress(:, :, iatom), 0.5_dp, fij, rij)
                     CALL virial_pair_force(atprop%atstress(:, :, jatom), 0.5_dp, fij, rij)
                  END IF
               END IF
            END IF
            DEALLOCATE (gcij)

         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

         IF (calculate_forces .AND. SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimg
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
            END DO
         END IF
      END IF

      CALL get_qs_env(qs_env, nkind=nkind)
      ALLOCATE (zeffk(nkind), xgamma(nkind))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
         CALL get_xtb_atom_param(xtb_kind, xgamma=xgamma(ikind), zeff=zeffk(ikind))
      END DO
      ! Diagonal 3rd order correction (DFTB3)
!deb  CALL build_dftb3_diagonal(qs_env, ks_matrix, rho, mcharge, energy, xgamma, zeffk, &
!deb                            calculate_forces, just_energy)
      DEALLOCATE (zeffk, xgamma)

      DEALLOCATE (gmcharge, gchrg, atom_of_kind, kind_of)

      CALL timestop(handle)

   END SUBROUTINE build_xtb_coulomb

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
   SUBROUTINE spme_evaluate(ewald_env, ewald_pw, particle_set, box, &
                            gmcharge, mcharge, calculate_forces, virial, use_virial, atprop)

      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      TYPE(cell_type), POINTER                           :: box
      REAL(KIND=dp), DIMENSION(:, :), INTENT(inout)      :: gmcharge
      REAL(KIND=dp), DIMENSION(:), INTENT(inout)         :: mcharge
      LOGICAL, INTENT(in)                                :: calculate_forces
      TYPE(virial_type), POINTER                         :: virial
      LOGICAL, INTENT(in)                                :: use_virial
      TYPE(atprop_type), POINTER                         :: atprop

      CHARACTER(len=*), PARAMETER :: routineN = 'spme_evaluate', routineP = moduleN//':'//routineN

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
            atprop%atstress(1, 1, p1) = atprop%atstress(1, 1, p1)-0.5_dp*mcharge(p1)*fint*dvols
            atprop%atstress(2, 2, p1) = atprop%atstress(2, 2, p1)-0.5_dp*mcharge(p1)*fint*dvols
            atprop%atstress(3, 3, p1) = atprop%atstress(3, 3, p1)-0.5_dp*mcharge(p1)*fint*dvols
         END DO

         CALL pw_pool_create_pw(pw_pool, phib_g, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         ffa = (0.5_dp/alpha)**2
         ffb = 1.0_dp/fourpi
         DO i = 1, 3
            DO ig = grid_spme%first_gne0, grid_spme%ngpts_cut_local
               phib_g%cc(ig) = ffb*dphi_g(i)%pw%cc(ig)*(1.0_dp-ffa*grid_spme%gsq(ig))
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
                  atprop%atstress(i, j, p1) = atprop%atstress(i, j, p1)-fint*dvols*mcharge(p1)
                  IF (i /= j) atprop%atstress(j, i, p1) = atprop%atstress(j, i, p1)-fint*dvols*mcharge(p1)
               END DO

            END DO
         END DO
         CALL pw_pool_give_back_pw(pw_pool, phib_g)

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
         virial%pv_virial = virial%pv_virial-(ffa*f_stress+h_stress)/REAL(para_env%num_pe, dp)
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
         gmcharge(p1, 1) = gmcharge(p1, 1)+fint*dvols

         IF (calculate_forces) THEN
            CALL dg_sum_patch_force_3d(drpot, rhos, center(:, p1), fat)
            gmcharge(p1, 2) = gmcharge(p1, 2)-fat(1)*dvols
            gmcharge(p1, 3) = gmcharge(p1, 3)-fat(2)*dvols
            gmcharge(p1, 4) = gmcharge(p1, 4)-fat(3)*dvols
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

   END SUBROUTINE spme_evaluate

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
   SUBROUTINE ewald_overlap(gmcharge, mcharge, alpha, n_list, virial, use_virial, atprop)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(inout)      :: gmcharge
      REAL(KIND=dp), DIMENSION(:), INTENT(in)            :: mcharge
      REAL(KIND=dp), INTENT(in)                          :: alpha
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list
      TYPE(virial_type), POINTER                         :: virial
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(atprop_type), POINTER                         :: atprop

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ewald_overlap', routineP = moduleN//':'//routineN

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
            gmcharge(iatom, 1) = gmcharge(iatom, 1)+mcharge(jatom)*fr
            gmcharge(jatom, 1) = gmcharge(jatom, 1)+mcharge(iatom)*fr
            IF (nmat > 1) THEN
               dfr = -2._dp*alpha*EXP(-alpha*alpha*dr*dr)*oorootpi/dr-fr/dr
               dfr = -dfr/dr
               DO i = 2, nmat
                  gmcharge(iatom, i) = gmcharge(iatom, i)-rij(i-1)*mcharge(jatom)*dfr
                  gmcharge(jatom, i) = gmcharge(jatom, i)+rij(i-1)*mcharge(iatom)*dfr
               END DO
            END IF
            IF (use_virial) THEN
               pfr = -dfr*mcharge(iatom)*mcharge(jatom)
               CALL virial_pair_force(virial%pv_virial, -pfr, rij, rij)
               IF (atprop%stress) THEN
                  CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp*pfr, rij, rij)
                  CALL virial_pair_force(atprop%atstress(:, :, jatom), -0.5_dp*pfr, rij, rij)
               END IF
            END IF
         END IF

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      CALL timestop(handle)

   END SUBROUTINE ewald_overlap

END MODULE xtb_coulomb

