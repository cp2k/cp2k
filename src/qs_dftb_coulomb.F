!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of Coulomb contributions in DFTB
!> \author JGH
! **************************************************************************************************
MODULE qs_dftb_coulomb
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set,&
                                              is_hydrogen
   USE atprop_types,                    ONLY: atprop_array_init,&
                                              atprop_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type,&
                                              dftb_control_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_get_block_p,&
                                              dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop,&
                                              dbcsr_iterator_type,&
                                              dbcsr_p_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE ewald_environment_types,         ONLY: ewald_env_get,&
                                              ewald_environment_type
   USE ewald_methods_tb,                ONLY: tb_ewald_overlap,&
                                              tb_spme_evaluate
   USE ewald_pw_types,                  ONLY: ewald_pw_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE mathconstants,                   ONLY: oorootpi,&
                                              pi
   USE message_passing,                 ONLY: mp_sum
   USE particle_types,                  ONLY: particle_type
   USE pw_poisson_types,                ONLY: do_ewald_ewald,&
                                              do_ewald_none,&
                                              do_ewald_pme,&
                                              do_ewald_spme
   USE qmmm_tb_coulomb,                 ONLY: build_tb_coulomb_qmqm
   USE qs_dftb3_methods,                ONLY: build_dftb3_diagonal
   USE qs_dftb_types,                   ONLY: qs_dftb_atom_type,&
                                              qs_dftb_pairpot_type
   USE qs_dftb_utils,                   ONLY: compute_block_sk,&
                                              get_dftb_atom_param
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
   USE sap_kind_types,                  ONLY: clist_type,&
                                              release_sap_int,&
                                              sap_int_type
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb_coulomb'

   ! screening for gamma function
   REAL(dp), PARAMETER                    :: tol_gamma = 1.e-4_dp
   ! small real number
   REAL(dp), PARAMETER                    :: rtiny = 1.e-10_dp

   PUBLIC :: build_dftb_coulomb, gamma_rab_sr

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ks_matrix ...
!> \param rho ...
!> \param mcharge ...
!> \param energy ...
!> \param calculate_forces ...
!> \param just_energy ...
! **************************************************************************************************
   SUBROUTINE build_dftb_coulomb(qs_env, ks_matrix, rho, mcharge, energy, &
                                 calculate_forces, just_energy)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(qs_rho_type), POINTER                         :: rho
      REAL(dp), DIMENSION(:)                             :: mcharge
      TYPE(qs_energy_type), POINTER                      :: energy
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy

      CHARACTER(len=*), PARAMETER :: routineN = 'build_dftb_coulomb', &
         routineP = moduleN//':'//routineN

      INTEGER :: atom_i, atom_j, blk, ewald_type, handle, i, ia, iac, iatom, ic, icol, ikind, img, &
         irow, is, jatom, jkind, natom, natorb_a, natorb_b, nimg, nkind, nmat
      INTEGER, DIMENSION(3)                              :: cellind, periodic
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind, kind_of
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: defined, do_ewald, do_gamma_stress, &
                                                            found, hb_sr_damp, use_virial
      REAL(KIND=dp)                                      :: alpha, ddr, deth, dgam, dr, drm, drp, &
                                                            fi, ga, gb, gmat, gmij, hb_para, zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: xgamma, zeffk
      REAL(KIND=dp), DIMENSION(0:3)                      :: eta_a, eta_b
      REAL(KIND=dp), DIMENSION(3)                        :: fij, rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dsblock, gmcharge, ksblock, pblock, &
                                                            sblock
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: dsint
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
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
      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_kind, dftb_kind_a, dftb_kind_b
      TYPE(qs_dftb_pairpot_type), DIMENSION(:, :), &
         POINTER                                         :: dftb_potential
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (matrix_p, matrix_s, virial, atprop, dft_control)

      use_virial = .FALSE.

      IF (calculate_forces) THEN
         nmat = 4
      ELSE
         nmat = 1
      END IF

      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)
      ALLOCATE (gmcharge(natom, nmat))
      gmcharge = 0._dp

      CALL get_qs_env(qs_env, &
                      particle_set=particle_set, &
                      cell=cell, &
                      virial=virial, &
                      atprop=atprop, &
                      dft_control=dft_control, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      force=force, para_env=para_env)

      IF (calculate_forces) THEN
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      END IF

      NULLIFY (dftb_potential)
      CALL get_qs_env(qs_env=qs_env, dftb_potential=dftb_potential)
      NULLIFY (n_list)
      CALL get_qs_env(qs_env=qs_env, sab_orb=n_list)
      CALL neighbor_list_iterator_create(nl_iterator, n_list)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cellind)

         CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind_a)
         CALL get_dftb_atom_param(dftb_kind_a, &
                                  defined=defined, eta=eta_a, natorb=natorb_a)
         IF (.NOT. defined .OR. natorb_a < 1) CYCLE
         CALL get_qs_kind(qs_kind_set(jkind), dftb_parameter=dftb_kind_b)
         CALL get_dftb_atom_param(dftb_kind_b, &
                                  defined=defined, eta=eta_b, natorb=natorb_b)
         IF (.NOT. defined .OR. natorb_b < 1) CYCLE

         ! gamma matrix
         hb_sr_damp = dft_control%qs_control%dftb_control%hb_sr_damp
         IF (hb_sr_damp) THEN
            ! short range correction enabled only when iatom XOR jatom are hydrogens
            hb_sr_damp = is_hydrogen(particle_set(iatom)%atomic_kind) .NEQV. &
                         is_hydrogen(particle_set(jatom)%atomic_kind)
         END IF
         IF (hb_sr_damp) THEN
            hb_para = dft_control%qs_control%dftb_control%hb_sr_para
         ELSE
            hb_para = 0.0_dp
         END IF
         ga = eta_a(0)
         gb = eta_b(0)
         dr = SQRT(SUM(rij(:)**2))
         gmat = gamma_rab_sr(dr, ga, gb, hb_para)
         gmcharge(jatom, 1) = gmcharge(jatom, 1) + gmat*mcharge(iatom)
         IF (iatom /= jatom) THEN
            gmcharge(iatom, 1) = gmcharge(iatom, 1) + gmat*mcharge(jatom)
         END IF
         IF (calculate_forces .AND. (iatom /= jatom .OR. dr > 0.001_dp)) THEN
            ddr = 0.1_dp*dftb_potential(ikind, jkind)%dgrd
            drp = dr + ddr
            drm = dr - ddr
            dgam = 0.5_dp*(gamma_rab_sr(drp, ga, gb, hb_para) - gamma_rab_sr(drm, ga, gb, hb_para))/ddr
            DO i = 1, 3
               gmcharge(jatom, i + 1) = gmcharge(jatom, i + 1) - dgam*mcharge(iatom)*rij(i)/dr
               IF (dr > 0.001_dp) THEN
                  gmcharge(iatom, i + 1) = gmcharge(iatom, i + 1) + dgam*mcharge(jatom)*rij(i)/dr
               END IF
               IF (use_virial) THEN
                  fij(i) = -mcharge(iatom)*mcharge(jatom)*dgam*rij(i)/dr
               END IF
            END DO
            IF (use_virial) THEN
               fi = 1.0_dp
               IF (iatom == jatom) fi = 0.5_dp
               CALL virial_pair_force(virial%pv_virial, fi, fij, rij)
               IF (atprop%stress) THEN
                  CALL virial_pair_force(atprop%atstress(:, :, iatom), fi*0.5_dp, fij, rij)
                  CALL virial_pair_force(atprop%atstress(:, :, jatom), fi*0.5_dp, fij, rij)
               END IF
            END IF
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      IF (atprop%energy) THEN
         CALL get_qs_env(qs_env=qs_env, particle_set=particle_set)
         natom = SIZE(particle_set)
         CALL atprop_array_init(atprop%atecoul, natom)
      END IF

      ! 1/R contribution
      do_ewald = dft_control%qs_control%dftb_control%do_ewald
      IF (do_ewald) THEN
         ! Ewald sum
         NULLIFY (ewald_env, ewald_pw)
         CALL get_qs_env(qs_env=qs_env, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw)
         CALL get_cell(cell=cell, periodic=periodic, deth=deth)
         CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type)
         CALL get_qs_env(qs_env=qs_env, sab_tbe=n_list)
         CALL tb_ewald_overlap(gmcharge, mcharge, alpha, n_list, &
                               virial, use_virial, atprop)
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
            CALL tb_spme_evaluate(ewald_env, ewald_pw, particle_set, cell, &
                                  gmcharge, mcharge, calculate_forces, virial, use_virial, atprop)
         END SELECT
      ELSE
         ! direct sum
         CALL get_qs_env(qs_env=qs_env, &
                         local_particles=local_particles)
         DO ikind = 1, SIZE(local_particles%n_el)
            DO ia = 1, local_particles%n_el(ikind)
               iatom = local_particles%list(ikind)%array(ia)
               DO jatom = 1, iatom - 1
                  rij = particle_set(iatom)%r - particle_set(jatom)%r
                  rij = pbc(rij, cell)
                  dr = SQRT(SUM(rij(:)**2))
                  gmcharge(iatom, 1) = gmcharge(iatom, 1) + mcharge(jatom)/dr
                  gmcharge(jatom, 1) = gmcharge(jatom, 1) + mcharge(iatom)/dr
                  DO i = 2, nmat
                     gmcharge(iatom, i) = gmcharge(iatom, i) + rij(i - 1)*mcharge(jatom)/dr**3
                     gmcharge(jatom, i) = gmcharge(jatom, i) - rij(i - 1)*mcharge(iatom)/dr**3
                  END DO
               END DO
            END DO
         END DO
         CPASSERT(.NOT. use_virial)
      END IF

      CALL mp_sum(gmcharge(:, 1), para_env%group)

      IF (do_ewald) THEN
         ! add self charge interaction and background charge contribution
         gmcharge(:, 1) = gmcharge(:, 1) - 2._dp*alpha*oorootpi*mcharge(:)
         IF (ANY(periodic(:) == 1)) THEN
            gmcharge(:, 1) = gmcharge(:, 1) - pi/alpha**2/deth
         END IF
      END IF

      energy%hartree = energy%hartree + 0.5_dp*SUM(mcharge(:)*gmcharge(:, 1))
      IF (atprop%energy) THEN
         CALL get_qs_env(qs_env=qs_env, &
                         local_particles=local_particles)
         DO ikind = 1, SIZE(local_particles%n_el)
            CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
            CALL get_dftb_atom_param(dftb_kind, zeff=zeff)
            DO ia = 1, local_particles%n_el(ikind)
               iatom = local_particles%list(ikind)%array(ia)
               atprop%atecoul(iatom) = atprop%atecoul(iatom) + &
                                       0.5_dp*zeff*gmcharge(iatom, 1)
            END DO
         END DO
      END IF

      IF (calculate_forces) THEN
         ALLOCATE (atom_of_kind(natom), kind_of(natom))

         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  kind_of=kind_of, &
                                  atom_of_kind=atom_of_kind)

         gmcharge(:, 2) = gmcharge(:, 2)*mcharge(:)
         gmcharge(:, 3) = gmcharge(:, 3)*mcharge(:)
         gmcharge(:, 4) = gmcharge(:, 4)*mcharge(:)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atom_i = atom_of_kind(iatom)
            force(ikind)%rho_elec(1, atom_i) = force(ikind)%rho_elec(1, atom_i) - gmcharge(iatom, 2)
            force(ikind)%rho_elec(2, atom_i) = force(ikind)%rho_elec(2, atom_i) - gmcharge(iatom, 3)
            force(ikind)%rho_elec(3, atom_i) = force(ikind)%rho_elec(3, atom_i) - gmcharge(iatom, 4)
         END DO
      END IF

      do_gamma_stress = .FALSE.
      IF (.NOT. just_energy .AND. use_virial) THEN
         IF (dft_control%nimages == 1) do_gamma_stress = .TRUE.
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

         NULLIFY (sap_int)
         IF (do_gamma_stress) THEN
            ! derivative overlap integral (non collapsed)
            CALL dftb_dsint_list(qs_env, sap_int)
         END IF

         IF (nimg == 1) THEN
            ! no k-points; all matrices have been transformed to periodic bsf
            CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, sblock, blk)
               gmij = 0.5_dp*(gmcharge(irow, 1) + gmcharge(icol, 1))
               DO is = 1, SIZE(ks_matrix, 1)
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock - gmij*sblock
               END DO
               IF (calculate_forces) THEN
                  ikind = kind_of(irow)
                  atom_i = atom_of_kind(irow)
                  jkind = kind_of(icol)
                  atom_j = atom_of_kind(icol)
                  NULLIFY (pblock)
                  CALL dbcsr_get_block_p(matrix=matrix_p(1, 1)%matrix, &
                                         row=irow, col=icol, block=pblock, found=found)
                  CPASSERT(found)
                  DO i = 1, 3
                     NULLIFY (dsblock)
                     CALL dbcsr_get_block_p(matrix=matrix_s(1 + i, 1)%matrix, &
                                            row=irow, col=icol, block=dsblock, found=found)
                     CPASSERT(found)
                     fi = -2.0_dp*gmij*SUM(pblock*dsblock)
                     force(ikind)%rho_elec(i, atom_i) = force(ikind)%rho_elec(i, atom_i) + fi
                     force(jkind)%rho_elec(i, atom_j) = force(jkind)%rho_elec(i, atom_j) - fi
                     fij(i) = fi
                  END DO
               END IF
            ENDDO
            CALL dbcsr_iterator_stop(iter)
            !
            IF (do_gamma_stress) THEN
               DO ikind = 1, nkind
                  DO jkind = 1, nkind
                     iac = ikind + nkind*(jkind - 1)
                     IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) CYCLE
                     DO ia = 1, sap_int(iac)%nalist
                        IF (.NOT. ASSOCIATED(sap_int(iac)%alist(ia)%clist)) CYCLE
                        iatom = sap_int(iac)%alist(ia)%aatom
                        DO ic = 1, sap_int(iac)%alist(ia)%nclist
                           jatom = sap_int(iac)%alist(ia)%clist(ic)%catom
                           rij = sap_int(iac)%alist(ia)%clist(ic)%rac
                           dr = SQRT(SUM(rij(:)**2))
                           IF (dr > 1.e-6_dp) THEN
                              dsint => sap_int(iac)%alist(ia)%clist(ic)%acint
                              gmij = 0.5_dp*(gmcharge(iatom, 1) + gmcharge(jatom, 1))
                              icol = MAX(iatom, jatom)
                              irow = MIN(iatom, jatom)
                              NULLIFY (pblock)
                              CALL dbcsr_get_block_p(matrix=matrix_p(1, 1)%matrix, &
                                                     row=irow, col=icol, block=pblock, found=found)
                              CPASSERT(found)
                              DO i = 1, 3
                                 IF (irow == iatom) THEN
                                    fij(i) = -2.0_dp*gmij*SUM(pblock*dsint(:, :, i))
                                 ELSE
                                    fij(i) = -2.0_dp*gmij*SUM(TRANSPOSE(pblock)*dsint(:, :, i))
                                 END IF
                              END DO
                              fi = 1.0_dp
                              IF (iatom == jatom) fi = 0.5_dp
                              CALL virial_pair_force(virial%pv_virial, fi, fij, rij)
                              IF (atprop%stress) THEN
                                 CALL virial_pair_force(atprop%atstress(:, :, iatom), 0.5_dp*fi, fij, rij)
                                 CALL virial_pair_force(atprop%atstress(:, :, jatom), 0.5_dp*fi, fij, rij)
                              END IF
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
            END IF
         ELSE
            NULLIFY (n_list)
            CALL get_qs_env(qs_env=qs_env, sab_orb=n_list)
            CALL neighbor_list_iterator_create(nl_iterator, n_list)
            DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
               CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                      iatom=iatom, jatom=jatom, r=rij, cell=cellind)

               icol = MAX(iatom, jatom)
               irow = MIN(iatom, jatom)
               ic = cell_to_index(cellind(1), cellind(2), cellind(3))
               CPASSERT(ic > 0)

               gmij = 0.5_dp*(gmcharge(iatom, 1) + gmcharge(jatom, 1))

               NULLIFY (sblock)
               CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                      row=irow, col=icol, block=sblock, found=found)
               CPASSERT(found)
               DO is = 1, SIZE(ks_matrix, 1)
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock - gmij*sblock
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
                     CALL dbcsr_get_block_p(matrix=matrix_s(1 + i, ic)%matrix, &
                                            row=irow, col=icol, block=dsblock, found=found)
                     CPASSERT(found)
                     fi = -2.0_dp*gmij*SUM(pblock*dsblock)
                     force(ikind)%rho_elec(i, atom_i) = force(ikind)%rho_elec(i, atom_i) + fi
                     force(jkind)%rho_elec(i, atom_j) = force(jkind)%rho_elec(i, atom_j) - fi
                     fij(i) = fi
                  END DO
                  IF (use_virial) THEN
                     fi = 1.0_dp
                     IF (iatom == jatom) fi = 0.5_dp
                     CALL virial_pair_force(virial%pv_virial, fi, fij, rij)
                     IF (atprop%stress) THEN
                        CALL virial_pair_force(atprop%atstress(:, :, iatom), fi*0.5_dp, fij, rij)
                        CALL virial_pair_force(atprop%atstress(:, :, jatom), fi*0.5_dp, fij, rij)
                     END IF
                  END IF
               END IF
            END DO
            CALL neighbor_list_iterator_release(nl_iterator)
         END IF

         IF (calculate_forces .AND. SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimg
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
            END DO
         END IF
      END IF

      IF (dft_control%qs_control%dftb_control%dftb3_diagonal) THEN
         CALL get_qs_env(qs_env, nkind=nkind)
         ALLOCATE (zeffk(nkind), xgamma(nkind))
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
            CALL get_dftb_atom_param(dftb_kind, dudq=xgamma(ikind), zeff=zeffk(ikind))
         END DO
         ! Diagonal 3rd order correction (DFTB3)
         CALL build_dftb3_diagonal(qs_env, ks_matrix, rho, mcharge, energy, xgamma, zeffk, &
                                   sap_int, calculate_forces, just_energy)
         DEALLOCATE (zeffk, xgamma)
      END IF

      ! QMMM
      IF (qs_env%qmmm .AND. qs_env%qmmm_periodic) THEN
         CALL build_tb_coulomb_qmqm(qs_env, ks_matrix, rho, mcharge, energy, &
                                    calculate_forces, just_energy)
      END IF

      IF (do_gamma_stress) THEN
         CALL release_sap_int(sap_int)
      END IF

      IF (calculate_forces) THEN
         DEALLOCATE (atom_of_kind, kind_of)
      END IF
      DEALLOCATE (gmcharge)

      CALL timestop(handle)

   END SUBROUTINE build_dftb_coulomb

! **************************************************************************************************
!> \brief  Computes the short-range gamma parameter from exact Coulomb
!>         interaction of normalized exp(-a*r) charge distribution - 1/r
!> \param r ...
!> \param ga ...
!> \param gb ...
!> \param hb_para ...
!> \return ...
!> \par Literature
!>         Elstner et al, PRB 58 (1998) 7260
!> \par History
!>      10.2008 Axel Kohlmeyer - adding sr_damp
!>      08.2014 JGH - adding flexible exponent for damping
!> \version 1.1
! **************************************************************************************************
   FUNCTION gamma_rab_sr(r, ga, gb, hb_para) RESULT(gamma)
      REAL(dp), INTENT(in)                               :: r, ga, gb, hb_para
      REAL(dp)                                           :: gamma

      REAL(dp)                                           :: a, b, fac, g_sum

      gamma = 0.0_dp
      a = 3.2_dp*ga ! 3.2 = 16/5 in Eq. 18 and ff.
      b = 3.2_dp*gb
      g_sum = a + b
      IF (g_sum < tol_gamma) RETURN ! hardness screening
      IF (r < rtiny) THEN ! This is for short distances but non-onsite terms
         ! This gives also correct diagonal elements (a=b, r=0)
         gamma = 0.5_dp*(a*b/g_sum + (a*b)**2/g_sum**3)
         RETURN
      END IF
      !
      ! distinguish two cases: Gamma's are very close, e.g. for the same atom type,
      !                        and Gamma's are different
      !
      IF (ABS(a - b) < rtiny) THEN
         fac = 1.6_dp*r*a*b/g_sum*(1.0_dp + a*b/g_sum**2)
         gamma = -(48.0_dp + 33._dp*fac + (9.0_dp + fac)*fac**2)*EXP(-fac)/(48._dp*r)
      ELSE
         gamma = -EXP(-a*r)*(0.5_dp*a*b**4/(a**2 - b**2)**2 - &
                             (b**6 - 3._dp*a**2*b**4)/(r*(a**2 - b**2)**3)) - & ! a-> b
                 EXP(-b*r)*(0.5_dp*b*a**4/(b**2 - a**2)**2 - &
                            (a**6 - 3._dp*b**2*a**4)/(r*(b**2 - a**2)**3)) ! b-> a
      END IF
      !
      ! damping function for better short range hydrogen bonds.
      ! functional form from Hu H. et al., J. Phys. Chem. A 2007, 111, 5685-5691
      ! according to Elstner M, Theor. Chem. Acc. 2006, 116, 316-325,
      ! this should only be applied to a-b pairs involving hydrogen.
      IF (hb_para > 0.0_dp) THEN
         gamma = gamma*EXP(-(0.5_dp*(ga + gb))**hb_para*r*r)
      END IF
   END FUNCTION gamma_rab_sr

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param sap_int ...
! **************************************************************************************************
   SUBROUTINE dftb_dsint_list(qs_env, sap_int)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dftb_dsint_list', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, i, iac, iatom, ikind, ilist, jatom, jkind, jneighbor, llm, lmaxi, lmaxj, &
         n1, n2, natorb_a, natorb_b, ngrd, ngrdcut, nkind, nlist, nneighbor
      INTEGER, DIMENSION(3)                              :: cell
      LOGICAL                                            :: defined
      REAL(KIND=dp)                                      :: ddr, dgrd, dr
      REAL(KIND=dp), DIMENSION(3)                        :: drij, rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dsblock, dsblockm, smatij, smatji
      TYPE(clist_type), POINTER                          :: clist
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dftb_control_type), POINTER                   :: dftb_control
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_kind_a, dftb_kind_b
      TYPE(qs_dftb_pairpot_type), DIMENSION(:, :), &
         POINTER                                         :: dftb_potential
      TYPE(qs_dftb_pairpot_type), POINTER                :: dftb_param_ij, dftb_param_ji
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, nkind=nkind)
      CPASSERT(.NOT. ASSOCIATED(sap_int))
      ALLOCATE (sap_int(nkind*nkind))
      DO i = 1, nkind*nkind
         NULLIFY (sap_int(i)%alist, sap_int(i)%asort, sap_int(i)%aindex)
         sap_int(i)%nalist = 0
      END DO

      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control, &
                      sab_orb=sab_orb)

      dftb_control => dft_control%qs_control%dftb_control

      NULLIFY (dftb_potential)
      CALL get_qs_env(qs_env=qs_env, &
                      dftb_potential=dftb_potential)

      ! loop over all atom pairs with a non-zero overlap (sab_orb)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, iatom=iatom, &
                                jatom=jatom, nlist=nlist, ilist=ilist, nnode=nneighbor, &
                                inode=jneighbor, cell=cell, r=rij)
         iac = ikind + nkind*(jkind - 1)
         !
         CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind_a)
         CALL get_dftb_atom_param(dftb_kind_a, &
                                  defined=defined, lmax=lmaxi, natorb=natorb_a)
         IF (.NOT. defined .OR. natorb_a < 1) CYCLE
         CALL get_qs_kind(qs_kind_set(jkind), dftb_parameter=dftb_kind_b)
         CALL get_dftb_atom_param(dftb_kind_b, &
                                  defined=defined, lmax=lmaxj, natorb=natorb_b)
         IF (.NOT. defined .OR. natorb_b < 1) CYCLE

         dr = SQRT(SUM(rij(:)**2))

         ! integral list
         IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) THEN
            sap_int(iac)%a_kind = ikind
            sap_int(iac)%p_kind = jkind
            sap_int(iac)%nalist = nlist
            ALLOCATE (sap_int(iac)%alist(nlist))
            DO i = 1, nlist
               NULLIFY (sap_int(iac)%alist(i)%clist)
               sap_int(iac)%alist(i)%aatom = 0
               sap_int(iac)%alist(i)%nclist = 0
            END DO
         END IF
         IF (.NOT. ASSOCIATED(sap_int(iac)%alist(ilist)%clist)) THEN
            sap_int(iac)%alist(ilist)%aatom = iatom
            sap_int(iac)%alist(ilist)%nclist = nneighbor
            ALLOCATE (sap_int(iac)%alist(ilist)%clist(nneighbor))
            DO i = 1, nneighbor
               sap_int(iac)%alist(ilist)%clist(i)%catom = 0
            END DO
         END IF
         clist => sap_int(iac)%alist(ilist)%clist(jneighbor)
         clist%catom = jatom
         clist%cell = cell
         clist%rac = rij
         ALLOCATE (clist%acint(natorb_a, natorb_b, 3))
         NULLIFY (clist%achint)
         clist%acint = 0._dp
         clist%nsgf_cnt = 0
         NULLIFY (clist%sgf_list)

         ! retrieve information on S matrix
         dftb_param_ij => dftb_potential(ikind, jkind)
         dftb_param_ji => dftb_potential(jkind, ikind)
         ! assume table size and type is symmetric
         ngrd = dftb_param_ij%ngrd
         ngrdcut = dftb_param_ij%ngrdcut
         dgrd = dftb_param_ij%dgrd
         ddr = dgrd*0.1_dp
         CPASSERT(dftb_param_ij%llm == dftb_param_ji%llm)
         llm = dftb_param_ij%llm
         smatij => dftb_param_ij%smat
         smatji => dftb_param_ji%smat

         IF (NINT(dr/dgrd) <= ngrdcut) THEN
            IF (iatom == jatom .AND. dr < 0.001_dp) THEN
               ! diagonal block
            ELSE
               ! off-diagonal block
               n1 = natorb_a
               n2 = natorb_b
               ALLOCATE (dsblock(n1, n2), dsblockm(n1, n2))
               DO i = 1, 3
                  dsblock = 0._dp
                  dsblockm = 0.0_dp
                  drij = rij

                  drij(i) = rij(i) - ddr
                  CALL compute_block_sk(dsblockm, smatij, smatji, drij, ngrd, ngrdcut, dgrd, &
                                        llm, lmaxi, lmaxj, iatom, iatom)

                  drij(i) = rij(i) + ddr
                  CALL compute_block_sk(dsblock, smatij, smatji, drij, ngrd, ngrdcut, dgrd, &
                                        llm, lmaxi, lmaxj, iatom, iatom)

                  dsblock = dsblock - dsblockm
                  dsblock = dsblock/(2.0_dp*ddr)

                  clist%acint(1:n1, 1:n2, i) = -dsblock(1:n1, 1:n2)
               ENDDO
               DEALLOCATE (dsblock, dsblockm)
            END IF
         END IF

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      CALL timestop(handle)

   END SUBROUTINE dftb_dsint_list

END MODULE qs_dftb_coulomb

