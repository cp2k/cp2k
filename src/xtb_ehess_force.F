!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of forces for Coulomb contributions in response xTB
!> \author JGH
! **************************************************************************************************
MODULE xtb_ehess_force
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_get_block_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, dbcsr_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE ewald_environment_types,         ONLY: ewald_env_get,&
                                              ewald_environment_type
   USE ewald_methods_tb,                ONLY: tb_ewald_overlap,&
                                              tb_spme_zforce
   USE ewald_pw_types,                  ONLY: ewald_pw_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: oorootpi,&
                                              pi
   USE message_passing,                 ONLY: mp_sum
   USE particle_types,                  ONLY: particle_type
   USE pw_poisson_types,                ONLY: do_ewald_ewald,&
                                              do_ewald_none,&
                                              do_ewald_pme,&
                                              do_ewald_spme
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
   USE qs_rho_types,                    ONLY: qs_rho_type
   USE virial_types,                    ONLY: virial_type
   USE xtb_coulomb,                     ONLY: dgamma_rab_sr,&
                                              gamma_rab_sr
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_ehess_force'

   PUBLIC :: calc_xtb_ehess_force

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_p0 ...
!> \param matrix_p1 ...
!> \param charges0 ...
!> \param mcharge0 ...
!> \param charges1 ...
!> \param mcharge1 ...
!> \param debug_forces ...
! **************************************************************************************************
   SUBROUTINE calc_xtb_ehess_force(qs_env, matrix_p0, matrix_p1, charges0, mcharge0, &
                                   charges1, mcharge1, debug_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p0, matrix_p1
      REAL(KIND=dp), DIMENSION(:, :), INTENT(in)         :: charges0
      REAL(KIND=dp), DIMENSION(:), INTENT(in)            :: mcharge0
      REAL(KIND=dp), DIMENSION(:, :), INTENT(in)         :: charges1
      REAL(KIND=dp), DIMENSION(:), INTENT(in)            :: mcharge1
      LOGICAL, INTENT(IN)                                :: debug_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'calc_xtb_ehess_force'

      INTEGER :: atom_i, atom_j, blk, ewald_type, handle, i, ia, iatom, icol, ikind, iounit, irow, &
         j, jatom, jkind, la, lb, lmaxa, lmaxb, natom, natorb_a, natorb_b, ni, nimg, nj, nkind, &
         nmat, za, zb
      INTEGER, DIMENSION(25)                             :: laoa, laob
      INTEGER, DIMENSION(3)                              :: cellind, periodic
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind, kind_of
      LOGICAL                                            :: calculate_forces, defined, do_ewald, &
                                                            found, just_energy, use_virial
      REAL(KIND=dp)                                      :: alpha, deth, dr, etaa, etab, fi, gmij0, &
                                                            gmij1, kg, rcut, rcuta, rcutb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: xgamma
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gammab, gcij0, gcij1, gmcharge0, &
                                                            gmcharge1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: gchrg0, gchrg1
      REAL(KIND=dp), DIMENSION(3)                        :: fij, fodeb, rij
      REAL(KIND=dp), DIMENSION(5)                        :: kappaa, kappab
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dsblock, pblock0, pblock1, sblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b, xtb_kind

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_unit_nr(logger)

      CPASSERT(ASSOCIATED(matrix_p1))

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      cell=cell, &
                      rho=rho, &
                      energy=energy, &
                      virial=virial, &
                      atprop=atprop, &
                      dft_control=dft_control)

      calculate_forces = .TRUE.
      just_energy = .FALSE.
      use_virial = .FALSE.
      nmat = 4
      nimg = dft_control%nimages
      IF (nimg > 1) THEN
         CPABORT('xTB-sTDA forces for k-points not available')
      END IF

      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)
      ALLOCATE (gchrg0(natom, 5, nmat))
      gchrg0 = 0._dp
      ALLOCATE (gmcharge0(natom, nmat))
      gmcharge0 = 0._dp
      ALLOCATE (gchrg1(natom, 5, nmat))
      gchrg1 = 0._dp
      ALLOCATE (gmcharge1(natom, nmat))
      gmcharge1 = 0._dp

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
         ni = lmaxa + 1
         nj = lmaxb + 1
         ALLOCATE (gammab(ni, nj))
         rcut = rcuta + rcutb
         dr = SQRT(SUM(rij(:)**2))
         CALL gamma_rab_sr(gammab, dr, ni, kappaa, etaa, nj, kappab, etab, kg, rcut)
         gchrg0(iatom, 1:ni, 1) = gchrg0(iatom, 1:ni, 1) + MATMUL(gammab, charges0(jatom, 1:nj))
         gchrg1(iatom, 1:ni, 1) = gchrg1(iatom, 1:ni, 1) + MATMUL(gammab, charges1(jatom, 1:nj))
         IF (iatom /= jatom) THEN
            gchrg0(jatom, 1:nj, 1) = gchrg0(jatom, 1:nj, 1) + MATMUL(charges0(iatom, 1:ni), gammab)
            gchrg1(jatom, 1:nj, 1) = gchrg1(jatom, 1:nj, 1) + MATMUL(charges1(iatom, 1:ni), gammab)
         END IF
         IF (dr > 1.e-6_dp) THEN
            CALL dgamma_rab_sr(gammab, dr, ni, kappaa, etaa, nj, kappab, etab, kg, rcut)
            DO i = 1, 3
               gchrg0(iatom, 1:ni, i + 1) = gchrg0(iatom, 1:ni, i + 1) &
                                            + MATMUL(gammab, charges0(jatom, 1:nj))*rij(i)/dr
               gchrg1(iatom, 1:ni, i + 1) = gchrg1(iatom, 1:ni, i + 1) &
                                            + MATMUL(gammab, charges1(jatom, 1:nj))*rij(i)/dr
               IF (iatom /= jatom) THEN
                  gchrg0(jatom, 1:nj, i + 1) = gchrg0(jatom, 1:nj, i + 1) &
                                               - MATMUL(charges0(iatom, 1:ni), gammab)*rij(i)/dr
                  gchrg1(jatom, 1:nj, i + 1) = gchrg1(jatom, 1:nj, i + 1) &
                                               - MATMUL(charges1(iatom, 1:ni), gammab)*rij(i)/dr
               END IF
            END DO
         END IF
         DEALLOCATE (gammab)
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ! 1/R contribution

      IF (dft_control%qs_control%xtb_control%coulomb_lr) THEN
         do_ewald = dft_control%qs_control%xtb_control%do_ewald
         IF (do_ewald) THEN
            ! Ewald sum
            NULLIFY (ewald_env, ewald_pw)
            CALL get_qs_env(qs_env=qs_env, &
                            ewald_env=ewald_env, ewald_pw=ewald_pw)
            CALL get_cell(cell=cell, periodic=periodic, deth=deth)
            CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type)
            CALL get_qs_env(qs_env=qs_env, sab_tbe=n_list)
            CALL tb_ewald_overlap(gmcharge0, mcharge0, alpha, n_list, virial, use_virial, atprop)
            CALL tb_ewald_overlap(gmcharge1, mcharge1, alpha, n_list, virial, use_virial, atprop)
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
               CALL tb_spme_zforce(ewald_env, ewald_pw, particle_set, cell, gmcharge0, mcharge0)
               CALL tb_spme_zforce(ewald_env, ewald_pw, particle_set, cell, gmcharge1, mcharge1)
            END SELECT
         ELSE
            ! direct sum
            CALL get_qs_env(qs_env=qs_env, local_particles=local_particles)
            DO ikind = 1, SIZE(local_particles%n_el)
               DO ia = 1, local_particles%n_el(ikind)
                  iatom = local_particles%list(ikind)%array(ia)
                  DO jatom = 1, iatom - 1
                     rij = particle_set(iatom)%r - particle_set(jatom)%r
                     rij = pbc(rij, cell)
                     dr = SQRT(SUM(rij(:)**2))
                     IF (dr > 1.e-6_dp) THEN
                        gmcharge0(iatom, 1) = gmcharge0(iatom, 1) + mcharge0(jatom)/dr
                        gmcharge0(jatom, 1) = gmcharge0(jatom, 1) + mcharge0(iatom)/dr
                        gmcharge1(iatom, 1) = gmcharge1(iatom, 1) + mcharge1(jatom)/dr
                        gmcharge1(jatom, 1) = gmcharge1(jatom, 1) + mcharge1(iatom)/dr
                        DO i = 2, nmat
                           gmcharge0(iatom, i) = gmcharge0(iatom, i) + rij(i - 1)*mcharge0(jatom)/dr**3
                           gmcharge0(jatom, i) = gmcharge0(jatom, i) - rij(i - 1)*mcharge0(iatom)/dr**3
                           gmcharge1(iatom, i) = gmcharge1(iatom, i) + rij(i - 1)*mcharge1(jatom)/dr**3
                           gmcharge1(jatom, i) = gmcharge1(jatom, i) - rij(i - 1)*mcharge1(iatom)/dr**3
                        END DO
                     END IF
                  END DO
               END DO
            END DO
            CPASSERT(.NOT. use_virial)
         END IF
      END IF

      ! global sum of gamma*p arrays
      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      force=force, para_env=para_env)
      CALL mp_sum(gmcharge0(:, 1), para_env%group)
      CALL mp_sum(gchrg0(:, :, 1), para_env%group)
      CALL mp_sum(gmcharge1(:, 1), para_env%group)
      CALL mp_sum(gchrg1(:, :, 1), para_env%group)

      IF (dft_control%qs_control%xtb_control%coulomb_lr) THEN
         IF (do_ewald) THEN
            ! add self charge interaction and background charge contribution
            gmcharge0(:, 1) = gmcharge0(:, 1) - 2._dp*alpha*oorootpi*mcharge0(:)
            IF (ANY(periodic(:) == 1)) THEN
               gmcharge0(:, 1) = gmcharge0(:, 1) - pi/alpha**2/deth
            END IF
            gmcharge1(:, 1) = gmcharge1(:, 1) - 2._dp*alpha*oorootpi*mcharge1(:)
            IF (ANY(periodic(:) == 1)) THEN
               gmcharge1(:, 1) = gmcharge1(:, 1) - pi/alpha**2/deth
            END IF
         END IF
      END IF

      ALLOCATE (atom_of_kind(natom), kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               kind_of=kind_of, &
                               atom_of_kind=atom_of_kind)

      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         atom_i = atom_of_kind(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
         CALL get_xtb_atom_param(xtb_kind, lmax=ni)
         ni = ni + 1
         ! short range
         fij = 0.0_dp
         DO i = 1, 3
            fij(i) = SUM(charges0(iatom, 1:ni)*gchrg1(iatom, 1:ni, i + 1)) + &
                     SUM(charges1(iatom, 1:ni)*gchrg0(iatom, 1:ni, i + 1))
         END DO
         force(ikind)%rho_elec(1, atom_i) = force(ikind)%rho_elec(1, atom_i) - fij(1)
         force(ikind)%rho_elec(2, atom_i) = force(ikind)%rho_elec(2, atom_i) - fij(2)
         force(ikind)%rho_elec(3, atom_i) = force(ikind)%rho_elec(3, atom_i) - fij(3)
         ! long range
         fij = 0.0_dp
         DO i = 1, 3
            fij(i) = gmcharge1(iatom, i + 1)*mcharge0(iatom) + &
                     gmcharge0(iatom, i + 1)*mcharge1(iatom)
         END DO
         force(ikind)%rho_elec(1, atom_i) = force(ikind)%rho_elec(1, atom_i) - fij(1)
         force(ikind)%rho_elec(2, atom_i) = force(ikind)%rho_elec(2, atom_i) - fij(2)
         force(ikind)%rho_elec(3, atom_i) = force(ikind)%rho_elec(3, atom_i) - fij(3)
      END DO
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*dH[Pz]    ", fodeb
      END IF

      CALL get_qs_env(qs_env=qs_env, matrix_s_kp=matrix_s)

      IF (SIZE(matrix_p0) == 2) THEN
         CALL dbcsr_add(matrix_p0(1)%matrix, matrix_p0(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         CALL dbcsr_add(matrix_p1(1)%matrix, matrix_p1(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF

      ! no k-points; all matrices have been transformed to periodic bsf
      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
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
         ALLOCATE (gcij0(ni, nj))
         ALLOCATE (gcij1(ni, nj))
         DO i = 1, ni
            DO j = 1, nj
               la = laoa(i) + 1
               lb = laob(j) + 1
               gcij0(i, j) = 0.5_dp*(gchrg0(irow, la, 1) + gchrg0(icol, lb, 1))
               gcij1(i, j) = 0.5_dp*(gchrg1(irow, la, 1) + gchrg1(icol, lb, 1))
            END DO
         END DO
         gmij0 = 0.5_dp*(gmcharge0(irow, 1) + gmcharge0(icol, 1))
         gmij1 = 0.5_dp*(gmcharge1(irow, 1) + gmcharge1(icol, 1))
         atom_i = atom_of_kind(irow)
         atom_j = atom_of_kind(icol)
         NULLIFY (pblock0)
         CALL dbcsr_get_block_p(matrix=matrix_p0(1)%matrix, &
                                row=irow, col=icol, block=pblock0, found=found)
         CPASSERT(found)
         NULLIFY (pblock1)
         CALL dbcsr_get_block_p(matrix=matrix_p1(1)%matrix, &
                                row=irow, col=icol, block=pblock1, found=found)
         CPASSERT(found)
         DO i = 1, 3
            NULLIFY (dsblock)
            CALL dbcsr_get_block_p(matrix=matrix_s(1 + i, 1)%matrix, &
                                   row=irow, col=icol, block=dsblock, found=found)
            CPASSERT(found)
            ! short range
            fi = -2.0_dp*SUM(pblock0*dsblock*gcij1) - 2.0_dp*SUM(pblock1*dsblock*gcij0)
            force(ikind)%rho_elec(i, atom_i) = force(ikind)%rho_elec(i, atom_i) + fi
            force(jkind)%rho_elec(i, atom_j) = force(jkind)%rho_elec(i, atom_j) - fi
            ! long range
            fi = -2.0_dp*gmij1*SUM(pblock0*dsblock) - 2.0_dp*gmij0*SUM(pblock1*dsblock)
            force(ikind)%rho_elec(i, atom_i) = force(ikind)%rho_elec(i, atom_i) + fi
            force(jkind)%rho_elec(i, atom_j) = force(jkind)%rho_elec(i, atom_j) - fi
         END DO
         DEALLOCATE (gcij0, gcij1)
      ENDDO
      CALL dbcsr_iterator_stop(iter)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*H[P]*dS  ", fodeb
      END IF

      IF (dft_control%qs_control%xtb_control%tb3_interaction) THEN
         CALL get_qs_env(qs_env, nkind=nkind)
         ALLOCATE (xgamma(nkind))
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, xgamma=xgamma(ikind))
         END DO
         ! Diagonal 3rd order correction (DFTB3)
         IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
         CALL dftb3_diagonal_hessian_force(qs_env, mcharge0, mcharge1, &
                                           matrix_p0(1)%matrix, matrix_p1(1)%matrix, xgamma)
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*H3[P]    ", fodeb
         END IF
         DEALLOCATE (xgamma)
      END IF

      IF (SIZE(matrix_p0) == 2) THEN
         CALL dbcsr_add(matrix_p0(1)%matrix, matrix_p0(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
         CALL dbcsr_add(matrix_p1(1)%matrix, matrix_p1(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
      END IF

      ! QMMM
      IF (qs_env%qmmm .AND. qs_env%qmmm_periodic) THEN
         CPABORT("Not Available")
      END IF

      DEALLOCATE (gmcharge0, gchrg0, gmcharge1, gchrg1)
      DEALLOCATE (atom_of_kind, kind_of)

      CALL timestop(handle)

   END SUBROUTINE calc_xtb_ehess_force

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mcharge0 ...
!> \param mcharge1 ...
!> \param matrixp0 ...
!> \param matrixp1 ...
!> \param xgamma ...
! **************************************************************************************************
   SUBROUTINE dftb3_diagonal_hessian_force(qs_env, mcharge0, mcharge1, &
                                           matrixp0, matrixp1, xgamma)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), DIMENSION(:)                             :: mcharge0, mcharge1
      TYPE(dbcsr_type), POINTER                          :: matrixp0, matrixp1
      REAL(dp), DIMENSION(:)                             :: xgamma

      CHARACTER(len=*), PARAMETER :: routineN = 'dftb3_diagonal_hessian_force'

      INTEGER                                            :: atom_i, atom_j, blk, handle, i, icol, &
                                                            ikind, irow, jkind, natom
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind, kind_of
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: fi, gmijp, gmijq, ui, uj
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dsblock, p0block, p1block, sblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, natom=natom)
      ALLOCATE (atom_of_kind(natom), kind_of(natom))
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               kind_of=kind_of, atom_of_kind=atom_of_kind)
      CALL get_qs_env(qs_env=qs_env, force=force)
      ! no k-points; all matrices have been transformed to periodic bsf
      CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, sblock, blk)
         ikind = kind_of(irow)
         atom_i = atom_of_kind(irow)
         ui = xgamma(ikind)
         jkind = kind_of(icol)
         atom_j = atom_of_kind(icol)
         uj = xgamma(jkind)
         !
         gmijp = ui*mcharge0(irow)*mcharge1(irow) + uj*mcharge0(icol)*mcharge1(icol)
         gmijq = 0.5_dp*ui*mcharge0(irow)**2 + 0.5_dp*uj*mcharge0(icol)**2
         !
         NULLIFY (p0block)
         CALL dbcsr_get_block_p(matrix=matrixp0, &
                                row=irow, col=icol, block=p0block, found=found)
         CPASSERT(found)
         NULLIFY (p1block)
         CALL dbcsr_get_block_p(matrix=matrixp1, &
                                row=irow, col=icol, block=p1block, found=found)
         CPASSERT(found)
         DO i = 1, 3
            NULLIFY (dsblock)
            CALL dbcsr_get_block_p(matrix=matrix_s(1 + i)%matrix, &
                                   row=irow, col=icol, block=dsblock, found=found)
            CPASSERT(found)
            fi = gmijp*SUM(p0block*dsblock) + gmijq*SUM(p1block*dsblock)
            force(ikind)%rho_elec(i, atom_i) = force(ikind)%rho_elec(i, atom_i) + fi
            force(jkind)%rho_elec(i, atom_j) = force(jkind)%rho_elec(i, atom_j) - fi
         END DO
      ENDDO
      CALL dbcsr_iterator_stop(iter)
      DEALLOCATE (atom_of_kind, kind_of)

      CALL timestop(handle)

   END SUBROUTINE dftb3_diagonal_hessian_force

END MODULE xtb_ehess_force

