!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of Spin Polarisation contributions in xTB
!> \author JGH
! **************************************************************************************************
MODULE xtb_spinpol
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_type
   USE bibliography,                    ONLY: Neugebauer2023,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_block_p,&
                                              dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop,&
                                              dbcsr_iterator_type,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE mulliken,                        ONLY: ao_charges
   USE particle_types,                  ONLY: particle_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE sap_kind_types,                  ONLY: sap_int_type
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_spinpol'

   PUBLIC :: build_xtb_spinpol, xtb_spinpol_hessian, xtb_spinpol_hforce

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ks_matrix ...
!> \param matrix_p ...
!> \param energy ...
!> \param sap_int ...
!> \param calculate_forces ...
!> \param just_energy ...
! **************************************************************************************************
   SUBROUTINE build_xtb_spinpol(qs_env, ks_matrix, matrix_p, energy, &
                                sap_int, calculate_forces, just_energy)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix, matrix_p
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy

      CHARACTER(len=*), PARAMETER                        :: routineN = 'build_xtb_spinpol'

      INTEGER :: atom_a, atom_i, atom_j, handle, i, ia, iac, iatom, ib, ic, icol, ikind, iknd, &
         irow, jatom, jkind, jknd, la, lb, na, natom, natorb, nb, nimg, nkind, nsgf, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(25)                             :: lao
      INTEGER, DIMENSION(3)                              :: cellind
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: defined, found, use_virial
      REAL(KIND=dp)                                      :: dr, espin, fi, fval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: docg
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aocg, bocg, pam, pbm, wab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: wabk
      REAL(KIND=dp), DIMENSION(3)                        :: fij, rij
      REAL(KIND=dp), DIMENSION(3, 3)                     :: wall
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: aksb, bksb, dsblock, pamat, pbmat, sblock
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: dsint
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_matrix
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p_kp, matrix_s, matrix_s_kp
      TYPE(dbcsr_type), POINTER                          :: s_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
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

      energy%xtb_spinpol = 0.0_dp

      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins
      nimg = dft_control%nimages

      IF (nspins == 2) THEN

         CALL cite_reference(Neugebauer2023)

         CALL get_qs_env(qs_env, &
                         qs_kind_set=qs_kind_set, &
                         particle_set=particle_set, &
                         atomic_kind_set=atomic_kind_set, &
                         cell=cell, &
                         virial=virial, &
                         atprop=atprop)

         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  kind_of=kind_of, &
                                  atom_of_kind=atom_of_kind)

         use_virial = .FALSE.
         IF (calculate_forces) THEN
            use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
         END IF

         CALL get_qs_env(qs_env, nkind=nkind, natom=natom)
         CALL get_qs_kind_set(qs_kind_set, maxsgf=nsgf)
         CALL get_qs_env(qs_env, matrix_s_kp=matrix_s, para_env=para_env)

         ! expand parameters
         ALLOCATE (wabk(nsgf, nsgf, nkind))
         wabk = 0.0_dp
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, natorb=natorb, lao=lao, wall=wall)
            DO ia = 1, natorb
               la = lao(ia) + 1
               DO ib = 1, natorb
                  lb = lao(ib) + 1
                  wabk(ia, ib, ikind) = wall(la, lb)
               END DO
            END DO
         END DO

         ! Calculate charges
         ALLOCATE (aocg(nsgf, natom), bocg(nsgf, natom))
         aocg = 0.0_dp
         bocg = 0.0_dp
         IF (nimg > 1) THEN
            matrix_s_kp => matrix_s(:, :)
            matrix_p_kp => matrix_p(1:1, :)
            CALL ao_charges(matrix_p_kp, matrix_s_kp, aocg, para_env)
            matrix_p_kp => matrix_p(2:2, :)
            CALL ao_charges(matrix_p_kp, matrix_s_kp, bocg, para_env)
         ELSE
            s_matrix => matrix_s(1, 1)%matrix
            p_matrix => matrix_p(1:1, 1)
            CALL ao_charges(p_matrix, s_matrix, aocg, para_env)
            p_matrix => matrix_p(2:2, 1)
            CALL ao_charges(p_matrix, s_matrix, bocg, para_env)
         END IF

         ! calculate energy
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=na)
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, defined=defined, natorb=natorb)
            IF (.NOT. defined .OR. natorb < 1) CYCLE
            ALLOCATE (docg(natorb), wab(natorb, natorb))
            wab(1:natorb, 1:natorb) = wabk(1:natorb, 1:natorb, ikind)
            DO iatom = 1, na
               atom_a = atomic_kind_set(ikind)%atom_list(iatom)
               docg = 0.0_dp
               docg(1:natorb) = aocg(1:natorb, atom_a) - bocg(1:natorb, atom_a)
               espin = 0.5_dp*DOT_PRODUCT(docg, MATMUL(wab, docg))
               energy%xtb_spinpol = energy%xtb_spinpol + espin
               IF (atprop%energy) THEN
                  atprop%atecoul(iatom) = atprop%atecoul(iatom) + espin
               END IF
            END DO
            DEALLOCATE (docg, wab)
         END DO

         ! Forces and Virial
         IF (calculate_forces) THEN
            CALL get_qs_env(qs_env=qs_env, force=force)
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
                  CALL dbcsr_iterator_next_block(iter, irow, icol, sblock)
                  ikind = kind_of(irow)
                  atom_i = atom_of_kind(irow)
                  jkind = kind_of(icol)
                  atom_j = atom_of_kind(icol)

                  CALL dbcsr_get_block_p(matrix=matrix_p(1, 1)%matrix, &
                                         row=irow, col=icol, block=pamat, found=found)
                  CPASSERT(found)
                  CALL dbcsr_get_block_p(matrix=matrix_p(2, 1)%matrix, &
                                         row=irow, col=icol, block=pbmat, found=found)
                  CPASSERT(found)

                  na = SIZE(pamat, 1)
                  nb = SIZE(pamat, 2)

                  DO i = 1, 3
                     CALL dbcsr_get_block_p(matrix=matrix_s(1 + i, 1)%matrix, &
                                            row=irow, col=icol, block=dsblock, found=found)
                     CPASSERT(found)

                     CALL fupdate(fi, pamat, pbmat, dsblock, na, nb, &
                                  wabk(1:na, 1:na, ikind), wabk(1:nb, 1:nb, jkind), &
                                  aocg(1:na, irow), aocg(1:nb, icol), bocg(1:na, irow), bocg(1:nb, icol))

                     force(ikind)%rho_elec(i, atom_i) = force(ikind)%rho_elec(i, atom_i) + fi
                     force(jkind)%rho_elec(i, atom_j) = force(jkind)%rho_elec(i, atom_j) - fi
                  END DO

               END DO
               CALL dbcsr_iterator_stop(iter)
               ! use dsint list
               IF (use_virial .AND. 0 == 0) THEN
                  CPASSERT(ASSOCIATED(sap_int))
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
                                 icol = MAX(iatom, jatom)
                                 irow = MIN(iatom, jatom)
                                 CALL dbcsr_get_block_p(matrix=matrix_p(1, 1)%matrix, &
                                                        row=irow, col=icol, block=pamat, found=found)
                                 CPASSERT(found)
                                 CALL dbcsr_get_block_p(matrix=matrix_p(2, 1)%matrix, &
                                                        row=irow, col=icol, block=pbmat, found=found)
                                 CPASSERT(found)
                                 IF (irow == iatom) THEN
                                    na = SIZE(pamat, 1)
                                    nb = SIZE(pamat, 2)
                                    ALLOCATE (pam(na, nb), pbm(na, nb))
                                    pam(1:na, 1:nb) = pamat(1:na, 1:nb)
                                    pbm(1:na, 1:nb) = pbmat(1:na, 1:nb)
                                 ELSE
                                    na = SIZE(pamat, 2)
                                    nb = SIZE(pamat, 1)
                                    ALLOCATE (pam(na, nb), pbm(na, nb))
                                    pam(1:na, 1:nb) = TRANSPOSE(pamat(1:nb, 1:na))
                                    pbm(1:na, 1:nb) = TRANSPOSE(pbmat(1:nb, 1:na))
                                 END IF

                                 DO i = 1, 3
                                    CALL fupdate(fi, pam, pbm, dsint(:, :, i), na, nb, &
                                                 wabk(1:na, 1:na, ikind), wabk(1:nb, 1:nb, jkind), &
                                                 aocg(1:na, iatom), aocg(1:nb, jatom), &
                                                 bocg(1:na, iatom), bocg(1:nb, jatom))
                                    fij(i) = fi
                                 END DO
                                 fi = 1.0_dp
                                 IF (iatom == jatom) fi = 0.5_dp
                                 CALL virial_pair_force(virial%pv_virial, fi, fij, rij)
                                 DEALLOCATE (pam, pbm)

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

                  IF (irow == iatom) THEN
                     iknd = ikind
                     jknd = jkind
                     atom_i = atom_of_kind(iatom)
                     atom_j = atom_of_kind(jatom)
                     rij = rij
                  ELSE
                     iknd = jkind
                     jknd = ikind
                     atom_i = atom_of_kind(jatom)
                     atom_j = atom_of_kind(iatom)
                     rij = -rij
                  END IF
                  !
                  CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                         row=irow, col=icol, block=pamat, found=found)
                  CPASSERT(found)
                  CALL dbcsr_get_block_p(matrix=matrix_p(2, ic)%matrix, &
                                         row=irow, col=icol, block=pbmat, found=found)
                  CPASSERT(found)

                  na = SIZE(pamat, 1)
                  nb = SIZE(pamat, 2)

                  fij = 0.0_dp
                  DO i = 1, 3
                     CALL dbcsr_get_block_p(matrix=matrix_s(1 + i, ic)%matrix, &
                                            row=irow, col=icol, block=dsblock, found=found)
                     CPASSERT(found)

                     CALL fupdate(fi, pamat, pbmat, dsblock, na, nb, &
                                  wabk(1:na, 1:na, iknd), wabk(1:nb, 1:nb, jknd), &
                                  aocg(1:na, irow), aocg(1:nb, icol), bocg(1:na, irow), bocg(1:nb, icol))

                     force(iknd)%rho_elec(i, atom_i) = force(iknd)%rho_elec(i, atom_i) + fi
                     force(jknd)%rho_elec(i, atom_j) = force(jknd)%rho_elec(i, atom_j) - fi
                     fij(i) = fi
                  END DO
                  IF (use_virial) THEN
                     fi = 1.0_dp
                     IF (iatom == jatom) fi = 0.5_dp
                     CALL virial_pair_force(virial%pv_virial, fi, fij, rij)
                  END IF

               END DO
               CALL neighbor_list_iterator_release(nl_iterator)

            END IF
         END IF

         ! KS matrix
         IF (.NOT. just_energy) THEN
            IF (nimg > 1) THEN
               CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
               CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
            END IF
            IF (nimg == 1) THEN
               ! no k-points; all matrices have been transformed to periodic bsf
               CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
               DO WHILE (dbcsr_iterator_blocks_left(iter))
                  CALL dbcsr_iterator_next_block(iter, irow, icol, sblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(1, 1)%matrix, &
                                         row=irow, col=icol, block=aksb, found=found)
                  CPASSERT(found)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(2, 1)%matrix, &
                                         row=irow, col=icol, block=bksb, found=found)
                  CPASSERT(found)
                  na = SIZE(aksb, 1)
                  nb = SIZE(aksb, 2)
                  ikind = kind_of(irow)
                  jkind = kind_of(icol)
                  fval = 0.5_dp
                  CALL ksupdate(aksb, bksb, sblock, na, nb, fval, &
                                wabk(1:na, 1:na, ikind), wabk(1:nb, 1:nb, jkind), &
                                aocg(1:na, irow), aocg(1:nb, icol), bocg(1:na, irow), bocg(1:nb, icol))
               END DO
               CALL dbcsr_iterator_stop(iter)
            ELSE
               CALL get_qs_env(qs_env=qs_env, sab_orb=n_list)
               CALL neighbor_list_iterator_create(nl_iterator, n_list)
               DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
                  CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                         iatom=iatom, jatom=jatom, r=rij, cell=cellind)

                  icol = MAX(iatom, jatom)
                  irow = MIN(iatom, jatom)

                  ic = cell_to_index(cellind(1), cellind(2), cellind(3))
                  CPASSERT(ic > 0)

                  ikind = kind_of(irow)
                  jkind = kind_of(icol)

                  CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                         row=irow, col=icol, block=sblock, found=found)
                  CPASSERT(found)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(1, ic)%matrix, &
                                         row=irow, col=icol, block=aksb, found=found)
                  CPASSERT(found)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(2, ic)%matrix, &
                                         row=irow, col=icol, block=bksb, found=found)
                  CPASSERT(found)

                  na = SIZE(aksb, 1)
                  nb = SIZE(aksb, 2)
                  fval = 0.5_dp
                  CALL ksupdate(aksb, bksb, sblock, na, nb, fval, &
                                wabk(1:na, 1:na, ikind), wabk(1:nb, 1:nb, jkind), &
                                aocg(1:na, irow), aocg(1:nb, icol), bocg(1:na, irow), bocg(1:nb, icol))
               END DO
               CALL neighbor_list_iterator_release(nl_iterator)
            END IF

         END IF

         DEALLOCATE (wabk)
         DEALLOCATE (aocg, bocg)
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_xtb_spinpol

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ks_matrix ...
!> \param matrix_p1 ...
! **************************************************************************************************
   SUBROUTINE xtb_spinpol_hessian(qs_env, ks_matrix, matrix_p1)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_matrix, matrix_p1

      CHARACTER(len=*), PARAMETER :: routineN = 'xtb_spinpol_hessian'

      INTEGER                                            :: handle, ia, ib, icol, ikind, irow, &
                                                            jkind, la, lb, na, natom, natorb, nb, &
                                                            nimg, nkind, nsgf, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(25)                             :: lao
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: fval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aocg1, bocg1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: wabk
      REAL(KIND=dp), DIMENSION(3, 3)                     :: wall
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: aksb, bksb, sblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_matrix
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: s_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins
      nimg = dft_control%nimages

      IF (nimg /= 1) THEN
         CPABORT("No kpoints allowed in xTB response calculation")
      END IF

      IF (nspins == 2) THEN

         CALL get_qs_env(qs_env, &
                         qs_kind_set=qs_kind_set, &
                         atomic_kind_set=atomic_kind_set)

         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

         CALL get_qs_env(qs_env, nkind=nkind, natom=natom)
         CALL get_qs_kind_set(qs_kind_set, maxsgf=nsgf)
         CALL get_qs_env(qs_env, matrix_s_kp=matrix_s, para_env=para_env)

         ! expand parameters
         ALLOCATE (wabk(nsgf, nsgf, nkind))
         wabk = 0.0_dp
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, natorb=natorb, lao=lao, wall=wall)
            DO ia = 1, natorb
               la = lao(ia) + 1
               DO ib = 1, natorb
                  lb = lao(ib) + 1
                  wabk(ia, ib, ikind) = wall(la, lb)
               END DO
            END DO
         END DO

         ! Calculate response charges
         ALLOCATE (aocg1(nsgf, natom), bocg1(nsgf, natom))
         aocg1 = 0.0_dp
         bocg1 = 0.0_dp
         s_matrix => matrix_s(1, 1)%matrix
         p_matrix => matrix_p1(1:1)
         CALL ao_charges(p_matrix, s_matrix, aocg1, para_env)
         p_matrix => matrix_p1(2:2)
         CALL ao_charges(p_matrix, s_matrix, bocg1, para_env)
         aocg1 = 0.5_dp*aocg1
         bocg1 = 0.5_dp*bocg1

         CALL dbcsr_iterator_start(iter, s_matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, irow, icol, sblock)
            CALL dbcsr_get_block_p(matrix=ks_matrix(1)%matrix, &
                                   row=irow, col=icol, block=aksb, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=ks_matrix(2)%matrix, &
                                   row=irow, col=icol, block=bksb, found=found)
            CPASSERT(found)
            na = SIZE(aksb, 1)
            nb = SIZE(aksb, 2)
            ikind = kind_of(irow)
            jkind = kind_of(icol)
            fval = 1.00_dp
            CALL ksupdate(aksb, bksb, sblock, na, nb, fval, &
                          wabk(1:na, 1:na, ikind), wabk(1:nb, 1:nb, jkind), &
                          aocg1(1:na, irow), aocg1(1:nb, icol), &
                          bocg1(1:na, irow), bocg1(1:nb, icol))
         END DO
         CALL dbcsr_iterator_stop(iter)

         DEALLOCATE (wabk)
         DEALLOCATE (aocg1, bocg1)
      END IF

      CALL timestop(handle)

   END SUBROUTINE xtb_spinpol_hessian

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_p0 ...
!> \param matrix_p1 ...
! **************************************************************************************************
   SUBROUTINE xtb_spinpol_hforce(qs_env, matrix_p0, matrix_p1)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p0, matrix_p1

      CHARACTER(len=*), PARAMETER :: routineN = 'xtb_spinpol_hforce'

      INTEGER                                            :: atom_i, atom_j, handle, i, ia, ib, icol, &
                                                            ikind, irow, jkind, la, lb, na, natom, &
                                                            natorb, nb, nimg, nkind, nsgf, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(25)                             :: lao
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: fi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aocg, aocg1, bocg, bocg1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: wabk
      REAL(KIND=dp), DIMENSION(3, 3)                     :: wall
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dsblock, p0amat, p0bmat, p1amat, p1bmat, &
                                                            sblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, p_matrix
      TYPE(dbcsr_type), POINTER                          :: s_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins
      nimg = dft_control%nimages
      IF (nimg /= 1) THEN
         CPABORT("xTB response forces for spin polarisation Hamiltonian not available")
      END IF

      IF (nspins == 2) THEN

         CALL get_qs_env(qs_env, &
                         qs_kind_set=qs_kind_set, &
                         particle_set=particle_set, &
                         atomic_kind_set=atomic_kind_set)

         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  kind_of=kind_of, &
                                  atom_of_kind=atom_of_kind)

         CALL get_qs_env(qs_env, nkind=nkind, natom=natom)
         CALL get_qs_kind_set(qs_kind_set, maxsgf=nsgf)

         CALL get_qs_env(qs_env, matrix_s=matrix_s, para_env=para_env)

         ! expand parameters
         ALLOCATE (wabk(nsgf, nsgf, nkind))
         wabk = 0.0_dp
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, natorb=natorb, lao=lao, wall=wall)
            DO ia = 1, natorb
               la = lao(ia) + 1
               DO ib = 1, natorb
                  lb = lao(ib) + 1
                  wabk(ia, ib, ikind) = wall(la, lb)
               END DO
            END DO
         END DO

         ! Calculate charges
         s_matrix => matrix_s(1)%matrix
         ALLOCATE (aocg(nsgf, natom), bocg(nsgf, natom))
         aocg = 0.0_dp
         bocg = 0.0_dp
         p_matrix => matrix_p0(1:1)
         CALL ao_charges(p_matrix, s_matrix, aocg, para_env)
         p_matrix => matrix_p0(2:2)
         CALL ao_charges(p_matrix, s_matrix, bocg, para_env)
         ! Calculate response charges
         ALLOCATE (aocg1(nsgf, natom), bocg1(nsgf, natom))
         aocg1 = 0.0_dp
         bocg1 = 0.0_dp
         p_matrix => matrix_p1(1:1)
         CALL ao_charges(p_matrix, s_matrix, aocg1, para_env)
         p_matrix => matrix_p1(2:2)
         CALL ao_charges(p_matrix, s_matrix, bocg1, para_env)
         aocg1 = 0.5_dp*aocg1
         bocg1 = 0.5_dp*bocg1

         ! calculate forces
         CALL get_qs_env(qs_env=qs_env, force=force)
         ! no k-points; all matrices have been transformed to periodic bsf
         CALL dbcsr_iterator_start(iter, s_matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, irow, icol, sblock)
            ikind = kind_of(irow)
            atom_i = atom_of_kind(irow)
            jkind = kind_of(icol)
            atom_j = atom_of_kind(icol)

            CALL dbcsr_get_block_p(matrix=matrix_p0(1)%matrix, &
                                   row=irow, col=icol, block=p0amat, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=matrix_p0(2)%matrix, &
                                   row=irow, col=icol, block=p0bmat, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=matrix_p1(1)%matrix, &
                                   row=irow, col=icol, block=p1amat, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=matrix_p1(2)%matrix, &
                                   row=irow, col=icol, block=p1bmat, found=found)
            CPASSERT(found)

            na = SIZE(p0amat, 1)
            nb = SIZE(p0amat, 2)

            DO i = 1, 3
               CALL dbcsr_get_block_p(matrix=matrix_s(1 + i)%matrix, &
                                      row=irow, col=icol, block=dsblock, found=found)
               CPASSERT(found)

               fi = 0.0_dp
               CALL f2update(fi, p0amat, p0bmat, p1amat, p1bmat, dsblock, na, nb, &
                             wabk(1:na, 1:na, ikind), wabk(1:nb, 1:nb, jkind), &
                             aocg(1:na, irow), aocg(1:nb, icol), &
                             bocg(1:na, irow), bocg(1:nb, icol), &
                             aocg1(1:na, irow), aocg1(1:nb, icol), &
                             bocg1(1:na, irow), bocg1(1:nb, icol))

               force(ikind)%rho_elec(i, atom_i) = force(ikind)%rho_elec(i, atom_i) + fi
               force(jkind)%rho_elec(i, atom_j) = force(jkind)%rho_elec(i, atom_j) - fi
            END DO

         END DO
         CALL dbcsr_iterator_stop(iter)

      END IF

      CALL timestop(handle)

   END SUBROUTINE xtb_spinpol_hforce

! **************************************************************************************************
!> \brief ...
!> \param aksb ...
!> \param bksb ...
!> \param sb ...
!> \param na ...
!> \param nb ...
!> \param fval ...
!> \param wabi ...
!> \param wabj ...
!> \param qai ...
!> \param qaj ...
!> \param qbi ...
!> \param qbj ...
! **************************************************************************************************
   SUBROUTINE ksupdate(aksb, bksb, sb, na, nb, fval, &
                       wabi, wabj, qai, qaj, qbi, qbj)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: aksb, bksb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sb
      INTEGER, INTENT(IN)                                :: na, nb
      REAL(KIND=dp), INTENT(IN)                          :: fval
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: wabi, wabj
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: qai, qaj, qbi, qbj

      INTEGER                                            :: ia, ib
      REAL(KIND=dp), DIMENSION(na)                       :: dqa, wa
      REAL(KIND=dp), DIMENSION(na, nb)                   :: wqab
      REAL(KIND=dp), DIMENSION(nb)                       :: dqb, wb

      dqa = qai - qbi
      dqb = qaj - qbj
      wa = MATMUL(wabi, dqa)
      wb = MATMUL(wabj, dqb)
      DO ib = 1, nb
         DO ia = 1, na
            wqab(ia, ib) = fval*sb(ia, ib)*(wa(ia) + wb(ib))
         END DO
      END DO

      aksb = aksb + wqab
      bksb = bksb - wqab

   END SUBROUTINE ksupdate

! **************************************************************************************************
!> \brief ...
!> \param fij ...
!> \param pa ...
!> \param pb ...
!> \param ds ...
!> \param na ...
!> \param nb ...
!> \param wabi ...
!> \param wabj ...
!> \param qai ...
!> \param qaj ...
!> \param qbi ...
!> \param qbj ...
! **************************************************************************************************
   SUBROUTINE fupdate(fij, pa, pb, ds, na, nb, &
                      wabi, wabj, qai, qaj, qbi, qbj)

      REAL(KIND=dp), INTENT(OUT)                         :: fij
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: pa, pb, ds
      INTEGER, INTENT(IN)                                :: na, nb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: wabi, wabj
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: qai, qaj, qbi, qbj

      INTEGER                                            :: ia, ib
      REAL(KIND=dp), DIMENSION(na)                       :: dpsa, dqa, wa
      REAL(KIND=dp), DIMENSION(na, nb)                   :: dpab
      REAL(KIND=dp), DIMENSION(nb)                       :: dpsb, dqb, wb

      dqa = qai - qbi
      dqb = qaj - qbj
      wa = MATMUL(wabi, dqa)
      wb = MATMUL(wabj, dqb)
      dpab = pa - pb
      dpsa = 0.0_dp
      dpsb = 0.0_dp
      DO ib = 1, nb
         DO ia = 1, na
            dpsa(ia) = dpsa(ia) + dpab(ia, ib)*ds(ia, ib)
            dpsb(ib) = dpsb(ib) + dpab(ia, ib)*ds(ia, ib)
         END DO
      END DO

      fij = SUM(wa*dpsa) + SUM(wb*dpsb)

   END SUBROUTINE fupdate

! **************************************************************************************************
!> \brief ...
!> \param fij ...
!> \param p0a ...
!> \param p0b ...
!> \param p1a ...
!> \param p1b ...
!> \param ds ...
!> \param na ...
!> \param nb ...
!> \param wabi ...
!> \param wabj ...
!> \param qai ...
!> \param qaj ...
!> \param qbi ...
!> \param qbj ...
!> \param rai ...
!> \param raj ...
!> \param rbi ...
!> \param rbj ...
! **************************************************************************************************
   SUBROUTINE f2update(fij, p0a, p0b, p1a, p1b, ds, na, nb, wabi, wabj, &
                       qai, qaj, qbi, qbj, rai, raj, rbi, rbj)

      REAL(KIND=dp), INTENT(OUT)                         :: fij
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: p0a, p0b, p1a, p1b, ds
      INTEGER, INTENT(IN)                                :: na, nb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: wabi, wabj
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: qai, qaj, qbi, qbj, rai, raj, rbi, rbj

      INTEGER                                            :: ia, ib
      REAL(KIND=dp), DIMENSION(na)                       :: dpsa, dqa, wa
      REAL(KIND=dp), DIMENSION(na, nb)                   :: dpab
      REAL(KIND=dp), DIMENSION(nb)                       :: dpsb, dqb, wb

      fij = 0.0_dp

      dqa = qai - qbi
      dqb = qaj - qbj
      wa = MATMUL(wabi, dqa)
      wb = MATMUL(wabj, dqb)
      dpab = p1a - p1b
      dpsa = 0.0_dp
      dpsb = 0.0_dp
      DO ib = 1, nb
         DO ia = 1, na
            dpsa(ia) = dpsa(ia) + dpab(ia, ib)*ds(ia, ib)
            dpsb(ib) = dpsb(ib) + dpab(ia, ib)*ds(ia, ib)
         END DO
      END DO

      fij = fij + SUM(wa*dpsa) + SUM(wb*dpsb)

      dqa = rai - rbi
      dqb = raj - rbj
      wa = MATMUL(wabi, dqa)
      wb = MATMUL(wabj, dqb)
      dpab = p0a - p0b
      dpsa = 0.0_dp
      dpsb = 0.0_dp
      DO ib = 1, nb
         DO ia = 1, na
            dpsa(ia) = dpsa(ia) + dpab(ia, ib)*ds(ia, ib)
            dpsb(ib) = dpsb(ib) + dpab(ia, ib)*ds(ia, ib)
         END DO
      END DO

      fij = fij + SUM(wa*dpsa) + SUM(wb*dpsb)

   END SUBROUTINE f2update

END MODULE xtb_spinpol
