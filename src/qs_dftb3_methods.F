!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of DFTB3 Terms
!> \author JGH
! **************************************************************************************************
MODULE qs_dftb3_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_type
   USE cell_types,                      ONLY: cell_type
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
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_sum
   USE particle_types,                  ONLY: particle_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE sap_kind_types,                  ONLY: sap_int_type
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb3_methods'

   PUBLIC :: build_dftb3_diagonal

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ks_matrix ...
!> \param rho ...
!> \param mcharge ...
!> \param energy ...
!> \param xgamma ...
!> \param zeff ...
!> \param sap_int ...
!> \param calculate_forces ...
!> \param just_energy ...
! **************************************************************************************************
   SUBROUTINE build_dftb3_diagonal(qs_env, ks_matrix, rho, mcharge, energy, xgamma, zeff, &
                                   sap_int, calculate_forces, just_energy)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(qs_rho_type), POINTER                         :: rho
      REAL(dp), DIMENSION(:), INTENT(in)                 :: mcharge
      TYPE(qs_energy_type), POINTER                      :: energy
      REAL(dp), DIMENSION(:), INTENT(in)                 :: xgamma, zeff
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy

      CHARACTER(len=*), PARAMETER :: routineN = 'build_dftb3_diagonal'

      INTEGER                                            :: atom_i, atom_j, blk, handle, i, ia, iac, &
                                                            iatom, ic, icol, ikind, irow, is, &
                                                            jatom, jkind, natom, nimg, nkind
      INTEGER, DIMENSION(3)                              :: cellind
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind, kind_of
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found, use_virial
      REAL(KIND=dp)                                      :: dr, eb3, eloc, fi, gmij, ua, ui, uj
      REAL(KIND=dp), DIMENSION(3)                        :: fij, rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dsblock, ksblock, pblock, sblock
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: dsint
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)
      NULLIFY (atprop)

      ! Energy
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, atprop=atprop)

      eb3 = 0.0_dp
      CALL get_qs_env(qs_env=qs_env, local_particles=local_particles)
      DO ikind = 1, SIZE(local_particles%n_el)
         ua = xgamma(ikind)
         DO ia = 1, local_particles%n_el(ikind)
            iatom = local_particles%list(ikind)%array(ia)
            eloc = -1.0_dp/6.0_dp*ua*mcharge(iatom)**3
            eb3 = eb3 + eloc
            IF (atprop%energy) THEN
               ! we have to add the part not covered by 0.5*Tr(FP)
               eloc = -0.5_dp*eloc - 0.25_dp*ua*zeff(ikind)*mcharge(iatom)**2
               atprop%atecoul(iatom) = atprop%atecoul(iatom) + eloc
            END IF
         END DO
      END DO
      CALL get_qs_env(qs_env=qs_env, para_env=para_env)
      CALL mp_sum(eb3, para_env%group)
      energy%dftb3 = eb3

      ! Forces and Virial
      IF (calculate_forces) THEN
         CALL get_qs_env(qs_env=qs_env, matrix_s_kp=matrix_s, natom=natom, force=force, &
                         cell=cell, virial=virial, particle_set=particle_set)
         ! virial
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

         ALLOCATE (atom_of_kind(natom), kind_of(natom))
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  kind_of=kind_of, atom_of_kind=atom_of_kind)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
         IF (SIZE(matrix_p, 1) == 2) THEN
            DO ic = 1, SIZE(matrix_p, 2)
               CALL dbcsr_add(matrix_p(1, ic)%matrix, matrix_p(2, ic)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            END DO
         END IF
         !
         nimg = SIZE(matrix_p, 2)
         NULLIFY (cell_to_index)
         IF (nimg > 1) THEN
            NULLIFY (kpoints)
            CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
            CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
         END IF
         IF (nimg == 1) THEN
            ! no k-points; all matrices have been transformed to periodic bsf
            CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, sblock, blk)
               ikind = kind_of(irow)
               atom_i = atom_of_kind(irow)
               ui = xgamma(ikind)
               jkind = kind_of(icol)
               atom_j = atom_of_kind(icol)
               uj = xgamma(jkind)
               !
               gmij = -0.5_dp*(ui*mcharge(irow)**2 + uj*mcharge(icol)**2)
               !
               NULLIFY (pblock)
               CALL dbcsr_get_block_p(matrix=matrix_p(1, 1)%matrix, &
                                      row=irow, col=icol, block=pblock, found=found)
               CPASSERT(found)
               DO i = 1, 3
                  NULLIFY (dsblock)
                  CALL dbcsr_get_block_p(matrix=matrix_s(1 + i, 1)%matrix, &
                                         row=irow, col=icol, block=dsblock, found=found)
                  CPASSERT(found)
                  fi = -gmij*SUM(pblock*dsblock)
                  force(ikind)%rho_elec(i, atom_i) = force(ikind)%rho_elec(i, atom_i) + fi
                  force(jkind)%rho_elec(i, atom_j) = force(jkind)%rho_elec(i, atom_j) - fi
                  fij(i) = fi
               END DO
            END DO
            CALL dbcsr_iterator_stop(iter)
            ! use dsint list
            IF (use_virial) THEN
               CPASSERT(ASSOCIATED(sap_int))
               CALL get_qs_env(qs_env, nkind=nkind)
               DO ikind = 1, nkind
                  DO jkind = 1, nkind
                     iac = ikind + nkind*(jkind - 1)
                     IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) CYCLE
                     ui = xgamma(ikind)
                     uj = xgamma(jkind)
                     DO ia = 1, sap_int(iac)%nalist
                        IF (.NOT. ASSOCIATED(sap_int(iac)%alist(ia)%clist)) CYCLE
                        iatom = sap_int(iac)%alist(ia)%aatom
                        DO ic = 1, sap_int(iac)%alist(ia)%nclist
                           jatom = sap_int(iac)%alist(ia)%clist(ic)%catom
                           rij = sap_int(iac)%alist(ia)%clist(ic)%rac
                           dr = SQRT(SUM(rij(:)**2))
                           IF (dr > 1.e-6_dp) THEN
                              dsint => sap_int(iac)%alist(ia)%clist(ic)%acint
                              gmij = -0.5_dp*(ui*mcharge(iatom)**2 + uj*mcharge(jatom)**2)
                              icol = MAX(iatom, jatom)
                              irow = MIN(iatom, jatom)
                              NULLIFY (pblock)
                              CALL dbcsr_get_block_p(matrix=matrix_p(1, 1)%matrix, &
                                                     row=irow, col=icol, block=pblock, found=found)
                              CPASSERT(found)
                              DO i = 1, 3
                                 IF (irow == iatom) THEN
                                    fij(i) = -gmij*SUM(pblock*dsint(:, :, i))
                                 ELSE
                                    fij(i) = -gmij*SUM(TRANSPOSE(pblock)*dsint(:, :, i))
                                 END IF
                              END DO
                              fi = 1.0_dp
                              IF (iatom == jatom) fi = 0.5_dp
                              CALL virial_pair_force(virial%pv_virial, fi, fij, rij)
                              IF (atprop%stress) THEN
                                 CALL virial_pair_force(atprop%atstress(:, :, irow), fi*0.5_dp, fij, rij)
                                 CALL virial_pair_force(atprop%atstress(:, :, icol), fi*0.5_dp, fij, rij)
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

               dr = SQRT(SUM(rij**2))
               IF (iatom == jatom .AND. dr < 1.0e-6_dp) CYCLE

               icol = MAX(iatom, jatom)
               irow = MIN(iatom, jatom)

               ic = cell_to_index(cellind(1), cellind(2), cellind(3))
               CPASSERT(ic > 0)

               ikind = kind_of(iatom)
               atom_i = atom_of_kind(iatom)
               ui = xgamma(ikind)
               jkind = kind_of(jatom)
               atom_j = atom_of_kind(jatom)
               uj = xgamma(jkind)
               !
               gmij = -0.5_dp*(ui*mcharge(iatom)**2 + uj*mcharge(jatom)**2)
               !
               NULLIFY (pblock)
               CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                      row=irow, col=icol, block=pblock, found=found)
               CPASSERT(found)
               DO i = 1, 3
                  NULLIFY (dsblock)
                  CALL dbcsr_get_block_p(matrix=matrix_s(1 + i, ic)%matrix, &
                                         row=irow, col=icol, block=dsblock, found=found)
                  CPASSERT(found)
                  IF (irow == iatom) THEN
                     fi = -gmij*SUM(pblock*dsblock)
                  ELSE
                     fi = gmij*SUM(pblock*dsblock)
                  END IF
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

            END DO
            CALL neighbor_list_iterator_release(nl_iterator)
            !
         END IF

         DEALLOCATE (atom_of_kind, kind_of)
         IF (SIZE(matrix_p, 1) == 2) THEN
            DO ic = 1, SIZE(matrix_p, 2)
               CALL dbcsr_add(matrix_p(1, ic)%matrix, matrix_p(2, ic)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
            END DO
         END IF
      END IF

      ! KS matrix
      IF (.NOT. just_energy) THEN
         CALL get_qs_env(qs_env=qs_env, matrix_s_kp=matrix_s, natom=natom)
         ALLOCATE (kind_of(natom))
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
         !
         nimg = SIZE(ks_matrix, 2)
         NULLIFY (cell_to_index)
         IF (nimg > 1) THEN
            NULLIFY (kpoints)
            CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
            CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
         END IF

         IF (nimg == 1) THEN
            ! no k-points; all matrices have been transformed to periodic bsf
            CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, sblock, blk)
               ikind = kind_of(irow)
               ui = xgamma(ikind)
               jkind = kind_of(icol)
               uj = xgamma(jkind)
               gmij = -0.5_dp*(ui*mcharge(irow)**2 + uj*mcharge(icol)**2)
               DO is = 1, SIZE(ks_matrix, 1)
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock - 0.5_dp*gmij*sblock
               END DO
            END DO
            CALL dbcsr_iterator_stop(iter)
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

               ikind = kind_of(iatom)
               ui = xgamma(ikind)
               jkind = kind_of(jatom)
               uj = xgamma(jkind)
               gmij = -0.5_dp*(ui*mcharge(iatom)**2 + uj*mcharge(jatom)**2)
               !
               NULLIFY (sblock)
               CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                      row=irow, col=icol, block=sblock, found=found)
               CPASSERT(found)
               DO is = 1, SIZE(ks_matrix, 1)
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock - 0.5_dp*gmij*sblock
               END DO

            END DO
            CALL neighbor_list_iterator_release(nl_iterator)
            !
         END IF
         DEALLOCATE (kind_of)
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_dftb3_diagonal

END MODULE qs_dftb3_methods

