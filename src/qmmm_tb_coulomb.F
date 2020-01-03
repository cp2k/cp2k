!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of QMMM Coulomb contributions in TB
!> \author JGH
! **************************************************************************************************
MODULE qmmm_tb_coulomb
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
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
   USE ewald_methods_tb,                ONLY: tb_spme_evaluate
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
   USE qmmm_types_low,                  ONLY: qmmm_env_qm_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_tb_coulomb'

   PUBLIC :: build_tb_coulomb_qmqm

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
   SUBROUTINE build_tb_coulomb_qmqm(qs_env, ks_matrix, rho, mcharge, energy, &
                                    calculate_forces, just_energy)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(qs_rho_type), POINTER                         :: rho
      REAL(dp), DIMENSION(:)                             :: mcharge
      TYPE(qs_energy_type), POINTER                      :: energy
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy

      CHARACTER(len=*), PARAMETER :: routineN = 'build_tb_coulomb_qmqm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_i, atom_j, blk, ewald_type, handle, &
                                                            i, ia, iatom, ikind, jatom, jkind, &
                                                            natom, nmat
      INTEGER, DIMENSION(3)                              :: periodic
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind, kind_of
      LOGICAL                                            :: found, use_virial
      REAL(KIND=dp)                                      :: alpha, deth, dfr, dr, fi, fr, gmij
      REAL(KIND=dp), DIMENSION(3)                        :: fij, rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dsblock, gmcharge, ksblock, ksblock_2, &
                                                            pblock, sblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell, mm_cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qmmm_env_qm_type), POINTER                    :: qmmm_env_qm
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (matrix_p, matrix_s, virial, atprop, dft_control)

      use_virial = .FALSE.

      IF (calculate_forces) THEN
         nmat = 4
      ELSE
         nmat = 1
      END IF

      natom = SIZE(mcharge)
      ALLOCATE (gmcharge(natom, nmat))
      gmcharge = 0._dp

      CALL get_qs_env(qs_env, &
                      particle_set=particle_set, &
                      cell=cell, &
                      virial=virial, &
                      atprop=atprop, &
                      dft_control=dft_control)

      IF (calculate_forces) THEN
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      END IF

      ! Qm-QM long range correction for QMMM calculations
      ! no atomic energy evaluation
      CPASSERT(.NOT. atprop%energy)
      ! no stress tensor possible for QMMM
      CPASSERT(.NOT. use_virial)
      qmmm_env_qm => qs_env%qmmm_env_qm
      ewald_env => qmmm_env_qm%ewald_env
      ewald_pw => qmmm_env_qm%ewald_pw
      CALL get_qs_env(qs_env=qs_env, super_cell=mm_cell)
      CALL get_cell(cell=mm_cell, periodic=periodic, deth=deth)
      CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type)
      gmcharge = 0.0_dp
      ! direct sum for overlap and local correction
      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      local_particles=local_particles, &
                      force=force, para_env=para_env)
      DO ikind = 1, SIZE(local_particles%n_el)
         DO ia = 1, local_particles%n_el(ikind)
            iatom = local_particles%list(ikind)%array(ia)
            DO jatom = 1, iatom - 1
               rij = particle_set(iatom)%r - particle_set(jatom)%r
               ! no pbc(rij,mm_cell) at this point, we assume that QM particles are
               ! inside QM box and QM box << MM box
               dr = SQRT(SUM(rij(:)**2))
               ! local (unit cell) correction 1/R
               gmcharge(iatom, 1) = gmcharge(iatom, 1) - mcharge(jatom)/dr
               gmcharge(jatom, 1) = gmcharge(jatom, 1) - mcharge(iatom)/dr
               DO i = 2, nmat
                  gmcharge(iatom, i) = gmcharge(iatom, i) - rij(i - 1)*mcharge(jatom)/dr**3
                  gmcharge(jatom, i) = gmcharge(jatom, i) + rij(i - 1)*mcharge(iatom)/dr**3
               END DO
               ! overlap correction
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
            END DO
         END DO
      END DO

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
         CALL tb_spme_evaluate(ewald_env, ewald_pw, particle_set, mm_cell, &
                               gmcharge, mcharge, calculate_forces, virial, use_virial, atprop)
      END SELECT
      !
      CALL mp_sum(gmcharge(:, 1), para_env%group)
      !
      ! add self charge interaction and background charge contribution
      gmcharge(:, 1) = gmcharge(:, 1) - 2._dp*alpha*oorootpi*mcharge(:)
      IF (ANY(periodic(:) == 1)) THEN
         gmcharge(:, 1) = gmcharge(:, 1) - pi/alpha**2/deth
      END IF
      !
      energy%qmmm_el = energy%qmmm_el + 0.5_dp*SUM(mcharge(:)*gmcharge(:, 1))
      !
      IF (calculate_forces) THEN
         ALLOCATE (atom_of_kind(natom), kind_of(natom))
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  kind_of=kind_of, &
                                  atom_of_kind=atom_of_kind)
      END IF
      !
      IF (.NOT. just_energy) THEN
         CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
         CALL qs_rho_get(rho, rho_ao=matrix_p)

         IF (calculate_forces .AND. SIZE(matrix_p) == 2) THEN
            CALL dbcsr_add(matrix_p(1)%matrix, matrix_p(2)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         END IF

         CALL dbcsr_iterator_start(iter, ks_matrix(1, 1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, iatom, jatom, ksblock, blk)
            NULLIFY (sblock, ksblock_2)
            IF (SIZE(ks_matrix, 1) > 1) THEN
               CALL dbcsr_get_block_p(matrix=ks_matrix(2, 1)%matrix, &
                                      row=iatom, col=jatom, block=ksblock_2, found=found)
            END IF
            CALL dbcsr_get_block_p(matrix=matrix_s(1)%matrix, &
                                   row=iatom, col=jatom, block=sblock, found=found)
            gmij = 0.5_dp*(gmcharge(iatom, 1) + gmcharge(jatom, 1))
            ksblock = ksblock - gmij*sblock
            IF (SIZE(ks_matrix, 1) > 1) ksblock_2 = ksblock_2 - gmij*sblock
            IF (calculate_forces) THEN
               ikind = kind_of(iatom)
               atom_i = atom_of_kind(iatom)
               jkind = kind_of(jatom)
               atom_j = atom_of_kind(jatom)
               NULLIFY (pblock)
               CALL dbcsr_get_block_p(matrix=matrix_p(1)%matrix, &
                                      row=iatom, col=jatom, block=pblock, found=found)
               DO i = 1, 3
                  NULLIFY (dsblock)
                  CALL dbcsr_get_block_p(matrix=matrix_s(1 + i)%matrix, &
                                         row=iatom, col=jatom, block=dsblock, found=found)
                  fi = -2.0_dp*gmij*SUM(pblock*dsblock)
                  force(ikind)%rho_elec(i, atom_i) = force(ikind)%rho_elec(i, atom_i) + fi
                  force(jkind)%rho_elec(i, atom_j) = force(jkind)%rho_elec(i, atom_j) - fi
                  fij(i) = fi
               END DO
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)
         IF (calculate_forces .AND. SIZE(matrix_p) == 2) THEN
            CALL dbcsr_add(matrix_p(1)%matrix, matrix_p(2)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
         END IF
      END IF

      IF (calculate_forces) THEN
         DEALLOCATE (atom_of_kind, kind_of)
      END IF
      DEALLOCATE (gmcharge)

      CALL timestop(handle)

   END SUBROUTINE build_tb_coulomb_qmqm

END MODULE qmmm_tb_coulomb

