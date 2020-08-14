!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of electric field contributions in TB
!> \author JGH
! **************************************************************************************************
MODULE efield_tb_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type,&
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
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE message_passing,                 ONLY: mp_sum
   USE particle_types,                  ONLY: particle_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_period_efield_types,          ONLY: efield_berry_type,&
                                              init_efield_matrices
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'efield_tb_methods'

   PUBLIC :: efield_tb_matrix

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
   SUBROUTINE efield_tb_matrix(qs_env, ks_matrix, rho, mcharge, energy, calculate_forces, just_energy)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(qs_rho_type), POINTER                         :: rho
      REAL(dp), DIMENSION(:), INTENT(in)                 :: mcharge
      TYPE(qs_energy_type), POINTER                      :: energy
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy

      CHARACTER(len=*), PARAMETER                        :: routineN = 'efield_tb_matrix'

      INTEGER                                            :: handle
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      energy%efield = 0.0_dp
      CALL get_qs_env(qs_env, dft_control=dft_control)
      IF (dft_control%qs_control%dftb .OR. dft_control%qs_control%xtb) THEN
         IF (dft_control%apply_period_efield) THEN
            IF (dft_control%period_efield%displacement_field) THEN
               CALL dfield_tb_berry(qs_env, ks_matrix, rho, mcharge, energy, calculate_forces, just_energy)
            ELSE
               CALL efield_tb_berry(qs_env, ks_matrix, rho, mcharge, energy, calculate_forces, just_energy)
            END IF
         ELSE IF (dft_control%apply_efield) THEN
            CALL efield_tb_local(qs_env, ks_matrix, rho, mcharge, energy, calculate_forces, just_energy)
         ELSE IF (dft_control%apply_efield_field) THEN
            CPABORT("efield_filed")
         END IF
      ELSE
         CPABORT("This routine should only be called from TB")
      END IF

      CALL timestop(handle)

   END SUBROUTINE efield_tb_matrix

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
   SUBROUTINE efield_tb_local(qs_env, ks_matrix, rho, mcharge, energy, calculate_forces, just_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(qs_rho_type), POINTER                         :: rho
      REAL(dp), DIMENSION(:), INTENT(in)                 :: mcharge
      TYPE(qs_energy_type), POINTER                      :: energy
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'efield_tb_local'

      INTEGER                                            :: atom_a, atom_b, blk, handle, ia, icol, &
                                                            idir, ikind, irow, ispin, jkind, &
                                                            natom, nspin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      LOGICAL                                            :: do_kpoints, found, use_virial
      REAL(dp)                                           :: charge, fdir
      REAL(dp), DIMENSION(3)                             :: ci, fieldpol, fij, ria, rib
      REAL(dp), DIMENSION(:, :), POINTER                 :: ds_block, ks_block, p_block, s_block
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control, cell=cell, particle_set=particle_set)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, energy=energy, para_env=para_env)
      CALL get_qs_env(qs_env=qs_env, do_kpoints=do_kpoints, virial=virial)
      IF (do_kpoints) THEN
         CPABORT("Local electric field with kpoints not possible. Use Berry phase periodic version")
      END IF
      ! disable stress calculation
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      IF (use_virial) THEN
         CPABORT("Stress tensor for non-periodic E-field not possible")
      END IF

      fieldpol = dft_control%efield_fields(1)%efield%polarisation* &
                 dft_control%efield_fields(1)%efield%strength

      natom = SIZE(particle_set)
      ci = 0.0_dp
      DO ia = 1, natom
         charge = mcharge(ia)
         ria = particle_set(ia)%r
         ria = pbc(ria, cell)
         ci(:) = ci(:) + charge*ria(:)
      END DO
      energy%efield = -SUM(ci(:)*fieldpol(:))

      IF (.NOT. just_energy) THEN

         IF (calculate_forces) THEN
            CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, force=force)
            ALLOCATE (atom_of_kind(natom), kind_of(natom))
            CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
            IF (para_env%mepos == 0) THEN
               DO ia = 1, natom
                  charge = mcharge(ia)
                  ikind = kind_of(ia)
                  atom_a = atom_of_kind(ia)
                  force(ikind)%efield(1:3, atom_a) = -charge*fieldpol(:)
               END DO
            ELSE
               DO ia = 1, natom
                  ikind = kind_of(ia)
                  atom_a = atom_of_kind(ia)
                  force(ikind)%efield(1:3, atom_a) = 0.0_dp
               END DO
            END IF
            CALL qs_rho_get(rho, rho_ao=matrix_p)
         END IF

         ! Update KS matrix
         nspin = SIZE(ks_matrix, 1)
         NULLIFY (matrix_s)
         CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
         CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            NULLIFY (ks_block, s_block, p_block)
            CALL dbcsr_iterator_next_block(iter, irow, icol, s_block, blk)
            ria = particle_set(irow)%r
            ria = pbc(ria, cell)
            rib = particle_set(icol)%r
            rib = pbc(rib, cell)
            fdir = 0.5_dp*SUM(fieldpol(1:3)*(ria(1:3) + rib(1:3)))
            DO ispin = 1, nspin
               CALL dbcsr_get_block_p(matrix=ks_matrix(ispin, 1)%matrix, &
                                      row=irow, col=icol, BLOCK=ks_block, found=found)
               ks_block = ks_block + fdir*s_block
               CPASSERT(found)
            END DO
            IF (calculate_forces) THEN
               ikind = kind_of(irow)
               jkind = kind_of(icol)
               atom_a = atom_of_kind(irow)
               atom_b = atom_of_kind(icol)
               fij = 0.0_dp
               DO ispin = 1, nspin
                  CALL dbcsr_get_block_p(matrix=matrix_p(ispin)%matrix, &
                                         row=irow, col=icol, BLOCK=p_block, found=found)
                  CPASSERT(found)
                  DO idir = 1, 3
                     CALL dbcsr_get_block_p(matrix=matrix_s(idir + 1)%matrix, &
                                            row=irow, col=icol, BLOCK=ds_block, found=found)
                     CPASSERT(found)
                     fij(idir) = fij(idir) + SUM(p_block*ds_block)
                  END DO
               END DO
               fdir = SUM(ria(1:3)*fieldpol(1:3))
               force(ikind)%efield(1:3, atom_a) = force(ikind)%efield(1:3, atom_a) + fdir*fij(1:3)
               force(jkind)%efield(1:3, atom_b) = force(jkind)%efield(1:3, atom_b) - fdir*fij(1:3)
               fdir = SUM(rib(1:3)*fieldpol(1:3))
               force(ikind)%efield(1:3, atom_a) = force(ikind)%efield(1:3, atom_a) + fdir*fij(1:3)
               force(jkind)%efield(1:3, atom_b) = force(jkind)%efield(1:3, atom_b) - fdir*fij(1:3)
            END IF
         ENDDO
         CALL dbcsr_iterator_stop(iter)

         IF (calculate_forces) THEN
            DO ikind = 1, SIZE(atomic_kind_set)
               CALL mp_sum(force(ikind)%efield, para_env%group)
            END DO
            DEALLOCATE (atom_of_kind, kind_of)
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE efield_tb_local

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
   SUBROUTINE efield_tb_berry(qs_env, ks_matrix, rho, mcharge, energy, calculate_forces, just_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(qs_rho_type), POINTER                         :: rho
      REAL(dp), DIMENSION(:), INTENT(in)                 :: mcharge
      TYPE(qs_energy_type), POINTER                      :: energy
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'efield_tb_berry'

      COMPLEX(KIND=dp)                                   :: zdeta
      COMPLEX(KIND=dp), DIMENSION(3)                     :: zi(3)
      INTEGER                                            :: atom_a, atom_b, blk, handle, ia, iatom, &
                                                            ic, icol, idir, ikind, irow, is, &
                                                            ispin, jatom, jkind, natom, nimg, nspin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(3)                              :: cellind
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found, use_virial
      REAL(KIND=dp)                                      :: charge, dd, fdir
      REAL(KIND=dp), DIMENSION(3)                        :: fieldpol, fij, forcea, fpolvec, kvec, &
                                                            qi, rab, ria, rib
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: ds_block, ks_block, p_block, s_block
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, cell, particle_set)
      CALL get_qs_env(qs_env, dft_control=dft_control, cell=cell, &
                      particle_set=particle_set, virial=virial)
      NULLIFY (qs_kind_set, para_env, sab_orb)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, &
                      energy=energy, para_env=para_env, sab_orb=sab_orb)

      ! calculate stress only if forces requested also
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      use_virial = use_virial .AND. calculate_forces
      ! disable stress calculation
      IF (use_virial) THEN
         CPABORT("Stress tensor for periodic E-field not implemented")
      END IF

      fieldpol = dft_control%period_efield%polarisation
      fieldpol = fieldpol/SQRT(DOT_PRODUCT(fieldpol, fieldpol))
      fieldpol = -fieldpol*dft_control%period_efield%strength
      hmat = cell%hmat(:, :)/twopi
      DO idir = 1, 3
         fpolvec(idir) = fieldpol(1)*hmat(1, idir) + fieldpol(2)*hmat(2, idir) + fieldpol(3)*hmat(3, idir)
      END DO

      natom = SIZE(particle_set)
      nspin = SIZE(ks_matrix, 1)

      zi(:) = CMPLX(1._dp, 0._dp, dp)
      DO ia = 1, natom
         charge = mcharge(ia)
         ria = particle_set(ia)%r
         DO idir = 1, 3
            kvec(:) = twopi*cell%h_inv(idir, :)
            dd = SUM(kvec(:)*ria(:))
            zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)**charge
            zi(idir) = zi(idir)*zdeta
         END DO
      END DO
      qi = AIMAG(LOG(zi))
      energy%efield = -SUM(fpolvec(:)*qi(:))

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

         IF (calculate_forces) THEN
            CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, force=force)
            ALLOCATE (atom_of_kind(natom), kind_of(natom))
            CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
            IF (para_env%mepos == 0) THEN
               DO ia = 1, natom
                  charge = -mcharge(ia)
                  iatom = atom_of_kind(ia)
                  ikind = kind_of(ia)
                  force(ikind)%efield(:, iatom) = fieldpol(:)*charge
                  IF (use_virial) THEN
                     ria = particle_set(ia)%r
                     ria = pbc(ria, cell)
                     forcea(1:3) = fieldpol(1:3)*charge
                     CALL virial_pair_force(virial%pv_virial, 1.0_dp, forcea, ria)
                  END IF
               END DO
            ELSE
               DO ia = 1, natom
                  iatom = atom_of_kind(ia)
                  ikind = kind_of(ia)
                  force(ikind)%efield(:, iatom) = 0.0_dp
               END DO
            END IF
         END IF

         IF (nimg == 1) THEN
            ! no k-points; all matrices have been transformed to periodic bsf
            CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, s_block, blk)

               fdir = 0.0_dp
               ria = particle_set(irow)%r
               rib = particle_set(icol)%r
               DO idir = 1, 3
                  kvec(:) = twopi*cell%h_inv(idir, :)
                  dd = SUM(kvec(:)*ria(:))
                  zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)
                  fdir = fdir + fpolvec(idir)*AIMAG(LOG(zdeta))
                  dd = SUM(kvec(:)*rib(:))
                  zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)
                  fdir = fdir + fpolvec(idir)*AIMAG(LOG(zdeta))
               END DO

               DO is = 1, nspin
                  NULLIFY (ks_block)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ks_block, found=found)
                  CPASSERT(found)
                  ks_block = ks_block + 0.5_dp*fdir*s_block
               END DO
               IF (calculate_forces) THEN
                  ikind = kind_of(irow)
                  jkind = kind_of(icol)
                  atom_a = atom_of_kind(irow)
                  atom_b = atom_of_kind(icol)
                  fij = 0.0_dp
                  DO ispin = 1, nspin
                     CALL dbcsr_get_block_p(matrix=matrix_p(ispin, 1)%matrix, &
                                            row=irow, col=icol, BLOCK=p_block, found=found)
                     CPASSERT(found)
                     DO idir = 1, 3
                        CALL dbcsr_get_block_p(matrix=matrix_s(idir + 1, 1)%matrix, &
                                               row=irow, col=icol, BLOCK=ds_block, found=found)
                        CPASSERT(found)
                        fij(idir) = fij(idir) + SUM(p_block*ds_block)
                     END DO
                  END DO
                  force(ikind)%efield(1:3, atom_a) = force(ikind)%efield(1:3, atom_a) + fdir*fij(1:3)
                  force(jkind)%efield(1:3, atom_b) = force(jkind)%efield(1:3, atom_b) - fdir*fij(1:3)
               END IF

            ENDDO
            CALL dbcsr_iterator_stop(iter)
         ELSE
            CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
            DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
               CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                      iatom=iatom, jatom=jatom, r=rab, cell=cellind)

               icol = MAX(iatom, jatom)
               irow = MIN(iatom, jatom)

               ic = cell_to_index(cellind(1), cellind(2), cellind(3))
               CPASSERT(ic > 0)

               fdir = 0.0_dp
               ria = particle_set(irow)%r
               rib = particle_set(icol)%r
               DO idir = 1, 3
                  kvec(:) = twopi*cell%h_inv(idir, :)
                  dd = SUM(kvec(:)*ria(:))
                  zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)
                  fdir = fdir + fpolvec(idir)*AIMAG(LOG(zdeta))
                  dd = SUM(kvec(:)*rib(:))
                  zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)
                  fdir = fdir + fpolvec(idir)*AIMAG(LOG(zdeta))
               END DO

               NULLIFY (s_block)
               CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                      row=irow, col=icol, block=s_block, found=found)
               CPASSERT(found)
               DO is = 1, nspin
                  NULLIFY (ks_block)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                         row=irow, col=icol, block=ks_block, found=found)
                  CPASSERT(found)
                  ks_block = ks_block + 0.5_dp*fdir*s_block
               END DO
               IF (calculate_forces) THEN
                  atom_a = atom_of_kind(iatom)
                  atom_b = atom_of_kind(jatom)
                  fij = 0.0_dp
                  DO ispin = 1, nspin
                     CALL dbcsr_get_block_p(matrix=matrix_p(ispin, ic)%matrix, &
                                            row=irow, col=icol, BLOCK=p_block, found=found)
                     CPASSERT(found)
                     DO idir = 1, 3
                        CALL dbcsr_get_block_p(matrix=matrix_s(idir + 1, ic)%matrix, &
                                               row=irow, col=icol, BLOCK=ds_block, found=found)
                        CPASSERT(found)
                        fij(idir) = fij(idir) + SUM(p_block*ds_block)
                     END DO
                  END DO
                  IF (irow == iatom) fij = -fij
                  force(ikind)%efield(1:3, atom_a) = force(ikind)%efield(1:3, atom_a) - fdir*fij(1:3)
                  force(jkind)%efield(1:3, atom_b) = force(jkind)%efield(1:3, atom_b) + fdir*fij(1:3)
               END IF

            END DO
            CALL neighbor_list_iterator_release(nl_iterator)
         END IF

         IF (calculate_forces) THEN
            DO ikind = 1, SIZE(atomic_kind_set)
               CALL mp_sum(force(ikind)%efield, para_env%group)
            END DO
            DEALLOCATE (atom_of_kind, kind_of)
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE efield_tb_berry

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
   SUBROUTINE dfield_tb_berry(qs_env, ks_matrix, rho, mcharge, energy, calculate_forces, just_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(qs_rho_type), POINTER                         :: rho
      REAL(dp), DIMENSION(:), INTENT(in)                 :: mcharge
      TYPE(qs_energy_type), POINTER                      :: energy
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dfield_tb_berry'

      COMPLEX(KIND=dp)                                   :: zdeta
      COMPLEX(KIND=dp), DIMENSION(3)                     :: zi(3)
      INTEGER                                            :: atom_a, atom_b, blk, handle, i, ia, &
                                                            iatom, ic, icol, idir, ikind, irow, &
                                                            is, ispin, jatom, jkind, natom, nimg, &
                                                            nspin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(3)                              :: cellind
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found, use_virial
      REAL(KIND=dp)                                      :: charge, dd, ener_field, fdir, omega
      REAL(KIND=dp), DIMENSION(3)                        :: ci, cqi, dfilter, di, fieldpol, fij, &
                                                            hdi, kvec, qi, rab, ria, rib
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: ds_block, ks_block, p_block, s_block
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(efield_berry_type), POINTER                   :: efield
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, cell, particle_set)
      CALL get_qs_env(qs_env, dft_control=dft_control, cell=cell, &
                      particle_set=particle_set, virial=virial)
      NULLIFY (qs_kind_set, para_env, sab_orb)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, &
                      efield=efield, energy=energy, para_env=para_env, sab_orb=sab_orb)

      ! efield history
      CALL init_efield_matrices(efield)
      CALL set_qs_env(qs_env, efield=efield)

      ! calculate stress only if forces requested also
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      use_virial = use_virial .AND. calculate_forces
      ! disable stress calculation
      IF (use_virial) THEN
         CPABORT("Stress tensor for periodic E-field not implemented")
      END IF

      dfilter(1:3) = dft_control%period_efield%d_filter(1:3)

      fieldpol = dft_control%period_efield%polarisation
      fieldpol = fieldpol/SQRT(DOT_PRODUCT(fieldpol, fieldpol))
      fieldpol = fieldpol*dft_control%period_efield%strength

      omega = cell%deth
      hmat = cell%hmat(:, :)/twopi

      natom = SIZE(particle_set)
      nspin = SIZE(ks_matrix, 1)

      zi(:) = CMPLX(1._dp, 0._dp, dp)
      DO ia = 1, natom
         charge = mcharge(ia)
         ria = particle_set(ia)%r
         DO idir = 1, 3
            kvec(:) = twopi*cell%h_inv(idir, :)
            dd = SUM(kvec(:)*ria(:))
            zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)**charge
            zi(idir) = zi(idir)*zdeta
         END DO
      END DO
      qi = AIMAG(LOG(zi))

      ! make sure the total normalized polarization is within [-1:1]
      DO idir = 1, 3
         cqi(idir) = qi(idir)/omega
         IF (cqi(idir) > pi) cqi(idir) = cqi(idir) - twopi
         IF (cqi(idir) < -pi) cqi(idir) = cqi(idir) + twopi
         ! now check for log branch
         IF (calculate_forces) THEN
            IF (ABS(efield%polarisation(idir) - cqi(idir)) > pi) THEN
               di(idir) = (efield%polarisation(idir) - cqi(idir))/pi
               DO i = 1, 10
                  cqi(idir) = cqi(idir) + SIGN(1.0_dp, di(idir))*twopi
                  IF (ABS(efield%polarisation(idir) - cqi(idir)) < pi) EXIT
               END DO
            END IF
         END IF
         cqi(idir) = cqi(idir)*omega
      END DO
      DO idir = 1, 3
         ci(idir) = 0.0_dp
         DO i = 1, 3
            ci(idir) = ci(idir) + hmat(idir, i)*cqi(i)
         END DO
      END DO
      ! update the references
      IF (calculate_forces) THEN
         ener_field = SUM(ci)
         ! check for smoothness of energy surface
         IF (ABS(efield%field_energy - ener_field) > pi*ABS(SUM(hmat))) THEN
            CPWARN("Large change of e-field energy detected. Correct for non-smooth energy surface")
         END IF
         efield%field_energy = ener_field
         efield%polarisation(:) = cqi(:)/omega
      END IF

      ! Energy
      ener_field = 0.0_dp
      DO idir = 1, 3
         ener_field = ener_field + dfilter(idir)*(fieldpol(idir) - 2._dp*twopi/omega*ci(idir))**2
      END DO
      energy%efield = 0.25_dp/twopi*ener_field

      IF (.NOT. just_energy) THEN
         di(:) = -(fieldpol(:) - 2._dp*twopi/omega*ci(:))*dfilter(:)/omega

         CALL get_qs_env(qs_env=qs_env, matrix_s_kp=matrix_s)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)

         nimg = dft_control%nimages
         NULLIFY (cell_to_index)
         IF (nimg > 1) THEN
            NULLIFY (kpoints)
            CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
            CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
         END IF

         IF (calculate_forces) THEN
            CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, force=force)
            ALLOCATE (atom_of_kind(natom), kind_of(natom))
            CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
            IF (para_env%mepos == 0) THEN
               DO ia = 1, natom
                  charge = mcharge(ia)
                  iatom = atom_of_kind(ia)
                  ikind = kind_of(ia)
                  force(ikind)%efield(:, iatom) = force(ikind)%efield(:, iatom) + di(:)*charge
                  IF (use_virial) THEN
                     ria = particle_set(ia)%r
                     ria = pbc(ria, cell)
                     CALL virial_pair_force(virial%pv_virial, 1.0_dp, di*charge, ria)
                  END IF
               END DO
            END IF
         END IF

         IF (nimg == 1) THEN
            ! no k-points; all matrices have been transformed to periodic bsf
            CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, s_block, blk)

               DO idir = 1, 3
                  hdi(idir) = -SUM(di(1:3)*hmat(1:3, idir))
               END DO
               fdir = 0.0_dp
               ria = particle_set(irow)%r
               rib = particle_set(icol)%r
               DO idir = 1, 3
                  kvec(:) = twopi*cell%h_inv(idir, :)
                  dd = SUM(kvec(:)*ria(:))
                  zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)
                  fdir = fdir + hdi(idir)*AIMAG(LOG(zdeta))
                  dd = SUM(kvec(:)*rib(:))
                  zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)
                  fdir = fdir + hdi(idir)*AIMAG(LOG(zdeta))
               END DO

               DO is = 1, nspin
                  NULLIFY (ks_block)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ks_block, found=found)
                  CPASSERT(found)
                  ks_block = ks_block + 0.5_dp*fdir*s_block
               END DO
               IF (calculate_forces) THEN
                  ikind = kind_of(irow)
                  jkind = kind_of(icol)
                  atom_a = atom_of_kind(irow)
                  atom_b = atom_of_kind(icol)
                  fij = 0.0_dp
                  DO ispin = 1, nspin
                     CALL dbcsr_get_block_p(matrix=matrix_p(ispin, 1)%matrix, &
                                            row=irow, col=icol, BLOCK=p_block, found=found)
                     CPASSERT(found)
                     DO idir = 1, 3
                        CALL dbcsr_get_block_p(matrix=matrix_s(idir + 1, 1)%matrix, &
                                               row=irow, col=icol, BLOCK=ds_block, found=found)
                        CPASSERT(found)
                        fij(idir) = fij(idir) + SUM(p_block*ds_block)
                     END DO
                  END DO
                  force(ikind)%efield(1:3, atom_a) = force(ikind)%efield(1:3, atom_a) + fdir*fij(1:3)
                  force(jkind)%efield(1:3, atom_b) = force(jkind)%efield(1:3, atom_b) - fdir*fij(1:3)
               END IF

            ENDDO
            CALL dbcsr_iterator_stop(iter)
         ELSE
            CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
            DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
               CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                      iatom=iatom, jatom=jatom, r=rab, cell=cellind)

               icol = MAX(iatom, jatom)
               irow = MIN(iatom, jatom)

               ic = cell_to_index(cellind(1), cellind(2), cellind(3))
               CPASSERT(ic > 0)

               DO idir = 1, 3
                  hdi(idir) = -SUM(di(1:3)*hmat(1:3, idir))
               END DO
               fdir = 0.0_dp
               ria = particle_set(irow)%r
               rib = particle_set(icol)%r
               DO idir = 1, 3
                  kvec(:) = twopi*cell%h_inv(idir, :)
                  dd = SUM(kvec(:)*ria(:))
                  zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)
                  fdir = fdir + hdi(idir)*AIMAG(LOG(zdeta))
                  dd = SUM(kvec(:)*rib(:))
                  zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)
                  fdir = fdir + hdi(idir)*AIMAG(LOG(zdeta))
               END DO

               NULLIFY (s_block)
               CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                      row=irow, col=icol, block=s_block, found=found)
               CPASSERT(found)
               DO is = 1, nspin
                  NULLIFY (ks_block)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                         row=irow, col=icol, block=ks_block, found=found)
                  CPASSERT(found)
                  ks_block = ks_block + 0.5_dp*fdir*s_block
               END DO
               IF (calculate_forces) THEN
                  atom_a = atom_of_kind(iatom)
                  atom_b = atom_of_kind(jatom)
                  fij = 0.0_dp
                  DO ispin = 1, nspin
                     CALL dbcsr_get_block_p(matrix=matrix_p(ispin, ic)%matrix, &
                                            row=irow, col=icol, BLOCK=p_block, found=found)
                     CPASSERT(found)
                     DO idir = 1, 3
                        CALL dbcsr_get_block_p(matrix=matrix_s(idir + 1, ic)%matrix, &
                                               row=irow, col=icol, BLOCK=ds_block, found=found)
                        CPASSERT(found)
                        fij(idir) = fij(idir) + SUM(p_block*ds_block)
                     END DO
                  END DO
                  IF (irow == iatom) fij = -fij
                  force(ikind)%efield(1:3, atom_a) = force(ikind)%efield(1:3, atom_a) - fdir*fij(1:3)
                  force(jkind)%efield(1:3, atom_b) = force(jkind)%efield(1:3, atom_b) + fdir*fij(1:3)
               END IF

            END DO
            CALL neighbor_list_iterator_release(nl_iterator)
         END IF

         IF (calculate_forces) THEN
            DEALLOCATE (atom_of_kind, kind_of)
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE dfield_tb_berry

END MODULE efield_tb_methods
