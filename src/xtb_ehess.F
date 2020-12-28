!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of Coulomb Hessian contributions in xTB
!> \author JGH
! **************************************************************************************************
MODULE xtb_ehess
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_get_block_p,&
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
   USE mathconstants,                   ONLY: oorootpi,&
                                              pi
   USE message_passing,                 ONLY: mp_sum
   USE particle_types,                  ONLY: particle_type
   USE pw_poisson_types,                ONLY: do_ewald_ewald,&
                                              do_ewald_none,&
                                              do_ewald_pme,&
                                              do_ewald_spme
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE virial_types,                    ONLY: virial_type
   USE xtb_coulomb,                     ONLY: gamma_rab_sr
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xtb_ehess'

   PUBLIC :: xtb_coulomb_hessian

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ks_matrix ...
!> \param charges1 ...
!> \param mcharge1 ...
!> \param mcharge ...
! **************************************************************************************************
   SUBROUTINE xtb_coulomb_hessian(qs_env, ks_matrix, charges1, mcharge1, mcharge)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_matrix
      REAL(dp), DIMENSION(:, :)                          :: charges1
      REAL(dp), DIMENSION(:)                             :: mcharge1, mcharge

      CHARACTER(len=*), PARAMETER :: routineN = 'xtb_coulomb_hessian'

      INTEGER :: blk, ewald_type, handle, i, ia, iatom, icol, ikind, irow, is, j, jatom, jkind, &
         la, lb, lmaxa, lmaxb, natom, natorb_a, natorb_b, ni, nj, nkind, nmat, za, zb
      INTEGER, DIMENSION(25)                             :: laoa, laob
      INTEGER, DIMENSION(3)                              :: cellind, periodic
      INTEGER, DIMENSION(:), POINTER                     :: kind_of
      LOGICAL                                            :: defined, do_ewald, found
      REAL(KIND=dp)                                      :: alpha, deth, dr, etaa, etab, gmij, kg, &
                                                            rcut, rcuta, rcutb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: xgamma
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gammab, gcij, gmcharge
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: gchrg
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      REAL(KIND=dp), DIMENSION(5)                        :: kappaa, kappab
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: ksblock, sblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
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
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b, xtb_kind

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      matrix_s_kp=matrix_s, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      cell=cell, &
                      dft_control=dft_control)

      IF (dft_control%nimages /= 1) THEN
         CPABORT("No kpoints allowed in xTB response calculation")
      END IF

      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)
      nmat = 1
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
         ni = lmaxa + 1
         nj = lmaxb + 1
         ALLOCATE (gammab(ni, nj))
         rcut = rcuta + rcutb
         dr = SQRT(SUM(rij(:)**2))
         CALL gamma_rab_sr(gammab, dr, ni, kappaa, etaa, nj, kappab, etab, kg, rcut)
         gchrg(iatom, 1:ni, 1) = gchrg(iatom, 1:ni, 1) + MATMUL(gammab, charges1(jatom, 1:nj))
         IF (iatom /= jatom) THEN
            gchrg(jatom, 1:nj, 1) = gchrg(jatom, 1:nj, 1) + MATMUL(charges1(iatom, 1:ni), gammab)
         END IF
         DEALLOCATE (gammab)
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ! 1/R contribution

      do_ewald = dft_control%qs_control%xtb_control%do_ewald
      IF (do_ewald) THEN
         ! Ewald sum
         NULLIFY (ewald_env, ewald_pw)
         NULLIFY (virial, atprop)
         CALL get_qs_env(qs_env=qs_env, &
                         ewald_env=ewald_env, ewald_pw=ewald_pw)
         CALL get_cell(cell=cell, periodic=periodic, deth=deth)
         CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type)
         CALL get_qs_env(qs_env=qs_env, sab_tbe=n_list)
         CALL tb_ewald_overlap(gmcharge, mcharge1, alpha, n_list, virial, .FALSE., atprop)
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
                                  gmcharge, mcharge1, .FALSE., virial, .FALSE., atprop)
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
                  IF (dr > 1.e-6_dp) THEN
                     gmcharge(iatom, 1) = gmcharge(iatom, 1) + mcharge1(jatom)/dr
                     gmcharge(jatom, 1) = gmcharge(jatom, 1) + mcharge1(iatom)/dr
                  END IF
               END DO
            END DO
         END DO
      END IF

      ! global sum of gamma*p arrays
      CALL get_qs_env(qs_env=qs_env, para_env=para_env)
      CALL mp_sum(gmcharge(:, 1), para_env%group)
      CALL mp_sum(gchrg(:, :, 1), para_env%group)

      IF (do_ewald) THEN
         ! add self charge interaction and background charge contribution
         gmcharge(:, 1) = gmcharge(:, 1) - 2._dp*alpha*oorootpi*mcharge1(:)
         IF (ANY(periodic(:) == 1)) THEN
            gmcharge(:, 1) = gmcharge(:, 1) - pi/alpha**2/deth
         END IF
      END IF

      ALLOCATE (kind_of(natom))
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

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
               la = laoa(i) + 1
               lb = laob(j) + 1
               gcij(i, j) = gchrg(irow, la, 1) + gchrg(icol, lb, 1)
            END DO
         END DO
         gmij = gmcharge(irow, 1) + gmcharge(icol, 1)
         DO is = 1, SIZE(ks_matrix)
            NULLIFY (ksblock)
            CALL dbcsr_get_block_p(matrix=ks_matrix(is)%matrix, &
                                   row=irow, col=icol, block=ksblock, found=found)
            CPASSERT(found)
            ksblock = ksblock - gcij*sblock
            ksblock = ksblock - gmij*sblock
         END DO
         DEALLOCATE (gcij)
      ENDDO
      CALL dbcsr_iterator_stop(iter)

      IF (dft_control%qs_control%xtb_control%tb3_interaction) THEN
         CALL get_qs_env(qs_env, nkind=nkind)
         ALLOCATE (xgamma(nkind))
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, xgamma=xgamma(ikind))
         END DO
         ! Diagonal 3rd order correction (DFTB3)
         CALL dftb3_diagonal_hessian(qs_env, ks_matrix, mcharge, mcharge1, xgamma)
         DEALLOCATE (xgamma)
      END IF

      IF (qs_env%qmmm .AND. qs_env%qmmm_periodic) THEN
         CPABORT("QMMM not available in xTB response calculations")
      END IF

      DEALLOCATE (gmcharge, gchrg, kind_of)

      CALL timestop(handle)

   END SUBROUTINE xtb_coulomb_hessian

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ks_matrix ...
!> \param mcharge ...
!> \param mcharge1 ...
!> \param xgamma ...
! **************************************************************************************************
   SUBROUTINE dftb3_diagonal_hessian(qs_env, ks_matrix, mcharge, mcharge1, xgamma)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_matrix
      REAL(dp), DIMENSION(:)                             :: mcharge, mcharge1, xgamma

      CHARACTER(len=*), PARAMETER :: routineN = 'dftb3_diagonal_hessian'

      INTEGER                                            :: blk, handle, icol, ikind, irow, is, &
                                                            jkind, natom
      INTEGER, DIMENSION(:), POINTER                     :: kind_of
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: gmij, ui, uj
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: ksblock, sblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, matrix_s_kp=matrix_s, natom=natom)
      ALLOCATE (kind_of(natom))
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      ! no k-points; all matrices have been transformed to periodic bsf
      CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, sblock, blk)
         ikind = kind_of(irow)
         ui = xgamma(ikind)
         jkind = kind_of(icol)
         uj = xgamma(jkind)
         gmij = ui*mcharge(irow)*mcharge1(irow) + uj*mcharge(icol)*mcharge1(icol)
         DO is = 1, SIZE(ks_matrix)
            NULLIFY (ksblock)
            CALL dbcsr_get_block_p(matrix=ks_matrix(is)%matrix, &
                                   row=irow, col=icol, block=ksblock, found=found)
            CPASSERT(found)
            ksblock = ksblock + gmij*sblock
         END DO
      ENDDO
      CALL dbcsr_iterator_stop(iter)
      DEALLOCATE (kind_of)

      CALL timestop(handle)

   END SUBROUTINE dftb3_diagonal_hessian

END MODULE xtb_ehess

