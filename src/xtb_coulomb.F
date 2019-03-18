!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
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

      INTEGER :: atom_i, atom_j, blk, ewald_type, handle, i, ia, iatom, ic, icol, ikind, img, &
         irow, is, j, jatom, jkind, la, lb, lmaxa, lmaxb, natom, natorb_a, natorb_b, ni, nimg, nj, &
         nkind, nmat, za, zb
      INTEGER, DIMENSION(25)                             :: laoa, laob
      INTEGER, DIMENSION(3)                              :: cellind, periodic
      INTEGER, DIMENSION(5)                              :: occ
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind, kind_of
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: defined, do_ewald, found, use_virial
      REAL(KIND=dp)                                      :: alpha, deth, dr, ecsr, etaa, etab, fi, &
                                                            gmij, kg, rcut, rcuta, rcutb, zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: xgamma, zeffk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gammab, gcij, gmcharge
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: gchrg
      REAL(KIND=dp), DIMENSION(25)                       :: gcint
      REAL(KIND=dp), DIMENSION(3)                        :: fij, rij
      REAL(KIND=dp), DIMENSION(5)                        :: kappaa, kappab
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dsblock, ksblock, pblock, sblock
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
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b, xtb_kind

      CALL timeset(routineN, handle)

      NULLIFY (matrix_p, matrix_s, virial, atprop, dft_control)

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
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
      ! loop over all atom pairs with a non-zero overlap (sab_orb)
      kg = dft_control%qs_control%xtb_control%kg
      NULLIFY (n_list)
      CALL get_qs_env(qs_env=qs_env, sab_orb=n_list)
      CALL neighbor_list_iterator_create(nl_iterator, n_list)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cellind)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
         CALL get_xtb_atom_param(xtb_atom_a, defined=defined, natorb=natorb_a)
         IF (.NOT. defined .OR. natorb_a < 1) CYCLE
         CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
         CALL get_xtb_atom_param(xtb_atom_b, defined=defined, natorb=natorb_b)
         IF (.NOT. defined .OR. natorb_b < 1) CYCLE
         ! atomic parameters
         CALL get_xtb_atom_param(xtb_atom_a, eta=etaa, lmax=lmaxa, kappa=kappaa, rcut=rcuta)
         CALL get_xtb_atom_param(xtb_atom_b, eta=etab, lmax=lmaxb, kappa=kappab, rcut=rcutb)
         ! gamma matrix
         ni = lmaxa+1
         nj = lmaxb+1
         ALLOCATE (gammab(ni, nj))
         rcut = rcuta+rcutb
         dr = SQRT(SUM(rij(:)**2))
         CALL gamma_rab_sr(gammab, dr, ni, kappaa, etaa, nj, kappab, etab, kg, rcut)
         gchrg(iatom, 1:ni, 1) = gchrg(iatom, 1:ni, 1)+MATMUL(gammab, charges(jatom, 1:nj))
         IF (iatom /= jatom) THEN
            gchrg(jatom, 1:nj, 1) = gchrg(jatom, 1:nj, 1)+MATMUL(charges(iatom, 1:ni), gammab)
         END IF
         IF (calculate_forces) THEN
            IF (dr > 1.e-6_dp) THEN
               CALL dgamma_rab_sr(gammab, dr, ni, kappaa, etaa, nj, kappab, etab, kg, rcut)
               DO i = 1, 3
                  gchrg(iatom, 1:ni, i+1) = gchrg(iatom, 1:ni, i+1) &
                                            +MATMUL(gammab, charges(jatom, 1:nj))*rij(i)/dr
                  IF (iatom /= jatom) THEN
                     gchrg(jatom, 1:nj, i+1) = gchrg(jatom, 1:nj, i+1) &
                                               -MATMUL(charges(iatom, 1:ni), gammab)*rij(i)/dr
                  END IF
               END DO
               IF (use_virial) THEN
                  gcint(1:ni) = MATMUL(gammab, charges(jatom, 1:nj))
                  DO i = 1, 3
                     fij(i) = -SUM(charges(iatom, 1:ni)*gcint(1:ni))*rij(i)/dr
                  END DO
                  fi = 1.0_dp
                  IF (iatom == jatom) fi = 0.5_dp
                  CALL virial_pair_force(virial%pv_virial, fi, fij, rij)
                  IF (atprop%stress) THEN
                     CALL virial_pair_force(atprop%atstress(:, :, iatom), fi*0.5_dp, fij, rij)
                     CALL virial_pair_force(atprop%atstress(:, :, jatom), fi*0.5_dp, fij, rij)
                  END IF
               END IF
            END IF
         END IF
         DEALLOCATE (gammab)
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

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
         CALL tb_ewald_overlap(gmcharge, mcharge, alpha, n_list, virial, use_virial, atprop)
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
               DO jatom = 1, iatom-1
                  rij = particle_set(iatom)%r-particle_set(jatom)%r
                  rij = pbc(rij, cell)
                  dr = SQRT(SUM(rij(:)**2))
                  IF (dr > 1.e-6_dp) THEN
                     gmcharge(iatom, 1) = gmcharge(iatom, 1)+mcharge(jatom)/dr
                     gmcharge(jatom, 1) = gmcharge(jatom, 1)+mcharge(iatom)/dr
                     DO i = 2, nmat
                        gmcharge(iatom, i) = gmcharge(iatom, i)+rij(i-1)*mcharge(jatom)/dr**3
                        gmcharge(jatom, i) = gmcharge(jatom, i)-rij(i-1)*mcharge(iatom)/dr**3
                     END DO
                  END IF
               END DO
            END DO
         END DO
         CPASSERT(.NOT. use_virial)
      END IF

      ! global sum of gamma*p arrays
      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
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
         CALL get_xtb_atom_param(xtb_kind, lmax=ni)
         ni = ni+1
         ecsr = ecsr+SUM(charges(iatom, 1:ni)*gchrg(iatom, 1:ni, 1))
      END DO

      energy%hartree = energy%hartree+0.5_dp*ecsr
      energy%hartree = energy%hartree+0.5_dp*SUM(mcharge(:)*gmcharge(:, 1))

      IF (atprop%energy) THEN
         CALL get_qs_env(qs_env=qs_env, local_particles=local_particles)
         DO ikind = 1, SIZE(local_particles%n_el)
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, lmax=ni, occupation=occ)
            ni = ni+1
            zeff = SUM(REAL(occ, KIND=dp))
            DO ia = 1, local_particles%n_el(ikind)
               iatom = local_particles%list(ikind)%array(ia)
               atprop%atecoul(iatom) = atprop%atecoul(iatom)+ &
                                       0.5_dp*SUM(REAL(occ(1:ni), KIND=dp)*gchrg(iatom, 1:ni, 1))
               atprop%atecoul(iatom) = atprop%atecoul(iatom)+ &
                                       0.5_dp*zeff*gmcharge(iatom, 1)
            END DO
         END DO
      END IF

      IF (calculate_forces) THEN
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atom_i = atom_of_kind(iatom)
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, lmax=ni)
            ! short range
            ni = ni+1
            DO i = 1, 3
               fij(i) = SUM(charges(iatom, 1:ni)*gchrg(iatom, 1:ni, i+1))
            END DO
            force(ikind)%rho_elec(1, atom_i) = force(ikind)%rho_elec(1, atom_i)-fij(1)
            force(ikind)%rho_elec(2, atom_i) = force(ikind)%rho_elec(2, atom_i)-fij(2)
            force(ikind)%rho_elec(3, atom_i) = force(ikind)%rho_elec(3, atom_i)-fij(3)
            ! long range
            DO i = 1, 3
               fij(i) = gmcharge(iatom, i+1)*mcharge(iatom)
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

         IF (nimg == 1) THEN
            ! no k-points; all matrices have been transformed to periodic bsf
            CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, sblock, blk)
               ikind = kind_of(irow)
               jkind = kind_of(icol)

               ! atomic parameters
               CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
               CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
               CALL get_xtb_atom_param(xtb_atom_a, z=za, lao=laoa)
               CALL get_xtb_atom_param(xtb_atom_b, z=zb, lao=laob)

               ni = SIZE(sblock, 1)
               nj = SIZE(sblock, 2)
               ALLOCATE (gcij(ni, nj))
               DO i = 1, ni
                  DO j = 1, nj
                     la = laoa(i)+1
                     lb = laob(j)+1
                     gcij(i, j) = 0.5_dp*(gchrg(irow, la, 1)+gchrg(icol, lb, 1))
                  END DO
               END DO
               gmij = 0.5_dp*(gmcharge(irow, 1)+gmcharge(icol, 1))
               DO is = 1, SIZE(ks_matrix, 1)
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock-gcij*sblock
                  ksblock = ksblock-gmij*sblock
               END DO
               IF (calculate_forces) THEN
                  atom_i = atom_of_kind(irow)
                  atom_j = atom_of_kind(icol)
                  NULLIFY (pblock)
                  CALL dbcsr_get_block_p(matrix=matrix_p(1, 1)%matrix, &
                                         row=irow, col=icol, block=pblock, found=found)
                  CPASSERT(found)
                  DO i = 1, 3
                     NULLIFY (dsblock)
                     CALL dbcsr_get_block_p(matrix=matrix_s(1+i, 1)%matrix, &
                                            row=irow, col=icol, block=dsblock, found=found)
                     CPASSERT(found)
                     fij(i) = 0.0_dp
                     ! short range
                     fi = -2.0_dp*SUM(pblock*dsblock*gcij)
                     force(ikind)%rho_elec(i, atom_i) = force(ikind)%rho_elec(i, atom_i)+fi
                     force(jkind)%rho_elec(i, atom_j) = force(jkind)%rho_elec(i, atom_j)-fi
                     fij(i) = fij(i)+fi
                     ! long range
                     fi = -2.0_dp*gmij*SUM(pblock*dsblock)
                     force(ikind)%rho_elec(i, atom_i) = force(ikind)%rho_elec(i, atom_i)+fi
                     force(jkind)%rho_elec(i, atom_j) = force(jkind)%rho_elec(i, atom_j)-fi
                     fij(i) = fij(i)+fi
                  END DO
                  ! this assumes a MIC case
                  IF (use_virial) THEN
                     rij = particle_set(irow)%r-particle_set(icol)%r
                     rij = pbc(rij, cell)
                     dr = SQRT(SUM(rij(:)**2))
                     IF (dr > 1.e-6_dp) THEN
                        fi = -1.0_dp
                        CALL virial_pair_force(virial%pv_virial, fi, fij, rij)
                        IF (atprop%stress) THEN
                           CALL virial_pair_force(atprop%atstress(:, :, irow), fi*0.5_dp, fij, rij)
                           CALL virial_pair_force(atprop%atstress(:, :, icol), fi*0.5_dp, fij, rij)
                        END IF
                     END IF
                  END IF
               END IF
               DEALLOCATE (gcij)
            ENDDO
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

               NULLIFY (sblock)
               CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                      row=irow, col=icol, block=sblock, found=found)
               CPASSERT(found)

               ! atomic parameters
               CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
               CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
               CALL get_xtb_atom_param(xtb_atom_a, z=za, lao=laoa)
               CALL get_xtb_atom_param(xtb_atom_b, z=zb, lao=laob)

               ni = SIZE(sblock, 1)
               nj = SIZE(sblock, 2)
               ALLOCATE (gcij(ni, nj))
               DO i = 1, ni
                  DO j = 1, nj
                     IF (irow == iatom) THEN
                        la = laoa(i)+1
                        lb = laob(j)+1
                        gcij(i, j) = 0.5_dp*(gchrg(iatom, la, 1)+gchrg(jatom, lb, 1))
                     ELSE
                        la = laoa(j)+1
                        lb = laob(i)+1
                        gcij(i, j) = 0.5_dp*(gchrg(iatom, la, 1)+gchrg(jatom, lb, 1))
                     END IF
                  END DO
               END DO
               gmij = 0.5_dp*(gmcharge(iatom, 1)+gmcharge(jatom, 1))
               DO is = 1, SIZE(ks_matrix, 1)
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock-gcij*sblock
                  ksblock = ksblock-gmij*sblock
               END DO

               IF (calculate_forces) THEN
                  atom_i = atom_of_kind(iatom)
                  atom_j = atom_of_kind(jatom)
                  IF (irow == jatom) THEN
                     gmij = -gmij
                     gcij = -gcij
                  END IF
                  NULLIFY (pblock)
                  CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                         row=irow, col=icol, block=pblock, found=found)
                  CPASSERT(found)
                  DO i = 1, 3
                     NULLIFY (dsblock)
                     CALL dbcsr_get_block_p(matrix=matrix_s(1+i, ic)%matrix, &
                                            row=irow, col=icol, block=dsblock, found=found)
                     CPASSERT(found)
                     fij(i) = 0.0_dp
                     ! short range
                     fi = -2.0_dp*SUM(pblock*dsblock*gcij)
                     force(ikind)%rho_elec(i, atom_i) = force(ikind)%rho_elec(i, atom_i)+fi
                     force(jkind)%rho_elec(i, atom_j) = force(jkind)%rho_elec(i, atom_j)-fi
                     fij(i) = fij(i)+fi
                     ! long range
                     fi = -2.0_dp*gmij*SUM(pblock*dsblock)
                     force(ikind)%rho_elec(i, atom_i) = force(ikind)%rho_elec(i, atom_i)+fi
                     force(jkind)%rho_elec(i, atom_j) = force(jkind)%rho_elec(i, atom_j)-fi
                     fij(i) = fij(i)+fi
                  END DO
                  IF (use_virial) THEN
                     dr = SQRT(SUM(rij(:)**2))
                     IF (dr > 1.e-6_dp) THEN
                        fi = 1.0_dp
                        IF (iatom == jatom) fi = 0.5_dp
                        CALL virial_pair_force(virial%pv_virial, fi, fij, rij)
                        IF (atprop%stress) THEN
                           CALL virial_pair_force(atprop%atstress(:, :, iatom), fi*0.5_dp, fij, rij)
                           CALL virial_pair_force(atprop%atstress(:, :, jatom), fi*0.5_dp, fij, rij)
                        END IF
                     END IF
                  END IF
               END IF
               DEALLOCATE (gcij)

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

      IF (dft_control%qs_control%xtb_control%tb3_interaction) THEN
         CALL get_qs_env(qs_env, nkind=nkind)
         ALLOCATE (zeffk(nkind), xgamma(nkind))
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, xgamma=xgamma(ikind), zeff=zeffk(ikind))
         END DO
         ! Diagonal 3rd order correction (DFTB3)
         CALL build_dftb3_diagonal(qs_env, ks_matrix, rho, mcharge, energy, xgamma, zeffk, &
                                   calculate_forces, just_energy)
         DEALLOCATE (zeffk, xgamma)
      END IF

      ! QMMM
      IF (qs_env%qmmm .AND. qs_env%qmmm_periodic) THEN
         CALL build_tb_coulomb_qmqm(qs_env, ks_matrix, rho, mcharge, energy, &
                                    calculate_forces, just_energy)
      END IF

      DEALLOCATE (gmcharge, gchrg, atom_of_kind, kind_of)

      CALL timestop(handle)

   END SUBROUTINE build_xtb_coulomb

! **************************************************************************************************
!> \brief  Computes the short-range gamma parameter from
!>         Nataga-Mishimoto-Ohno-Klopman formula for xTB
!>         WARNING: The xTB function (gamma - 1/r) has still an l-dependent longrange (1/r^3)
!>                  behaviour. We use a cutoff function to smoothly remove this part.
!>                  However, this will change energies and effect final results.
!>
!> \param gmat ...
!> \param rab ...
!> \param nla ...
!> \param kappaa ...
!> \param etaa ...
!> \param nlb ...
!> \param kappab ...
!> \param etab ...
!> \param kg ...
!> \param rcut ...
!> \par History
!>      10.2018 JGH
!> \version 1.1
! **************************************************************************************************
   SUBROUTINE gamma_rab_sr(gmat, rab, nla, kappaa, etaa, nlb, kappab, etab, kg, rcut)
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: gmat
      REAL(dp), INTENT(IN)                               :: rab
      INTEGER, INTENT(IN)                                :: nla
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: kappaa
      REAL(dp), INTENT(IN)                               :: etaa
      INTEGER, INTENT(IN)                                :: nlb
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: kappab
      REAL(dp), INTENT(IN)                               :: etab, kg, rcut

      CHARACTER(len=*), PARAMETER :: routineN = 'gamma_rab_sr', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: rsmooth = 1.0_dp

      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: fcut, r, rk, x
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eta

      ALLOCATE (eta(nla, nlb))
      eta = 0.0_dp

      DO j = 1, nlb
         DO i = 1, nla
            eta(i, j) = 1._dp/(etaa*(1._dp+kappaa(i)))+1._dp/(etab*(1._dp+kappab(j)))
            eta(i, j) = 2._dp/eta(i, j)
         END DO
      END DO

      gmat = 0.0_dp
      IF (rab < 1.e-6_dp) THEN
         ! on site terms
         gmat(:, :) = eta(:, :)
      ELSEIF (rab > rcut) THEN
         ! do nothing
      ELSE
         rk = rab**kg
         eta = eta**(-kg)
         IF (rab < rcut-rsmooth) THEN
            fcut = 1.0_dp
         ELSE
            r = rab-(rcut-rsmooth)
            x = r/rsmooth
            fcut = -10._dp*x**5+15._dp*x**4-10._dp*x**3+1._dp
         END IF
         gmat(:, :) = fcut*(1._dp/(rk+eta(:, :)))**(1._dp/kg)-fcut/rab
      END IF

      DEALLOCATE (eta)

   END SUBROUTINE gamma_rab_sr

! **************************************************************************************************
!> \brief  Computes the derivative of the short-range gamma parameter from
!>         Nataga-Mishimoto-Ohno-Klopman formula for xTB
!>         WARNING: The xTB function (gamma - 1/r) has still an l-dependent longrange (1/r^3)
!>                  behaviour. We use a cutoff function to smoothly remove this part.
!>                  However, this will change energies and effect final results.
!>
!> \param dgmat ...
!> \param rab ...
!> \param nla ...
!> \param kappaa ...
!> \param etaa ...
!> \param nlb ...
!> \param kappab ...
!> \param etab ...
!> \param kg ...
!> \param rcut ...
!> \par History
!>      10.2018 JGH
!> \version 1.1
! **************************************************************************************************
   SUBROUTINE dgamma_rab_sr(dgmat, rab, nla, kappaa, etaa, nlb, kappab, etab, kg, rcut)
      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: dgmat
      REAL(dp), INTENT(IN)                               :: rab
      INTEGER, INTENT(IN)                                :: nla
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: kappaa
      REAL(dp), INTENT(IN)                               :: etaa
      INTEGER, INTENT(IN)                                :: nlb
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: kappab
      REAL(dp), INTENT(IN)                               :: etab, kg, rcut

      CHARACTER(len=*), PARAMETER :: routineN = 'dgamma_rab_sr', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: rsmooth = 1.0_dp

      INTEGER                                            :: i, j
      REAL(KIND=dp)                                      :: dfcut, fcut, r, rk, x
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eta

      ALLOCATE (eta(nla, nlb))

      DO j = 1, nlb
         DO i = 1, nla
            eta(i, j) = 1._dp/(etaa*(1._dp+kappaa(i)))+1._dp/(etab*(1._dp+kappab(j)))
            eta(i, j) = 2._dp/eta(i, j)
         END DO
      END DO

      IF (rab < 1.e-6) THEN
         ! on site terms
         dgmat(:, :) = 0.0_dp
      ELSEIF (rab > rcut) THEN
         dgmat(:, :) = 0.0_dp
      ELSE
         eta = eta**(-kg)
         rk = rab**kg
         IF (rab < rcut-rsmooth) THEN
            fcut = 1.0_dp
            dfcut = 0.0_dp
         ELSE
            r = rab-(rcut-rsmooth)
            x = r/rsmooth
            fcut = -10._dp*x**5+15._dp*x**4-10._dp*x**3+1._dp
            dfcut = -50._dp*x**4+60._dp*x**3-30._dp*x**2
            dfcut = dfcut/rsmooth
         END IF
         dgmat(:, :) = dfcut*(1._dp/(rk+eta(:, :)))**(1._dp/kg)
         dgmat(:, :) = dgmat(:, :)-dfcut/rab+fcut/rab**2
         dgmat(:, :) = dgmat(:, :)-fcut/(rk+eta(:, :))*(1._dp/(rk+eta(:, :)))**(1._dp/kg)*rk/rab
      END IF

      DEALLOCATE (eta)

   END SUBROUTINE dgamma_rab_sr

END MODULE xtb_coulomb

