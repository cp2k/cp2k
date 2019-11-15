!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE qs_mixing_utils

   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_symmetric
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_max,&
                                              mp_min
   USE pw_types,                        ONLY: pw_p_type
   USE qs_density_mixing_types,         ONLY: broyden_mixing_new_nr,&
                                              broyden_mixing_nr,&
                                              gspace_mixing_nr,&
                                              mixing_storage_type,&
                                              multisecant_mixing_nr,&
                                              pulay_mixing_nr
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_atom_types,               ONLY: rho_atom_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_methods,                  ONLY: cp_sm_mix
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_mixing_utils'

   PUBLIC :: mixing_allocate, mixing_init, charge_mixing_init, &
             self_consistency_check

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param rho_ao ...
!> \param p_delta ...
!> \param para_env ...
!> \param p_out ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE self_consistency_check(rho_ao, p_delta, para_env, p_out, delta)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao, p_delta
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: p_out
      REAL(KIND=dp), INTENT(INOUT)                       :: delta

      CHARACTER(len=*), PARAMETER :: routineN = 'self_consistency_check', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ic, ispin, nimg, nspins
      REAL(KIND=dp)                                      :: tmp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_q, p_in

      CALL timeset(routineN, handle)

      NULLIFY (matrix_q, p_in)

      CPASSERT(ASSOCIATED(p_out))
      NULLIFY (matrix_q, p_in)
      p_in => rho_ao
      matrix_q => p_delta
      nspins = SIZE(p_in, 1)
      nimg = SIZE(p_in, 2)

      ! Compute the difference (p_out - p_in)and check convergence
      delta = 0.0_dp
      DO ispin = 1, nspins
         DO ic = 1, nimg
            CALL dbcsr_set(matrix_q(ispin, ic)%matrix, 0.0_dp)
            CALL cp_sm_mix(m1=p_out(ispin, ic)%matrix, m2=p_in(ispin, ic)%matrix, &
                           p_mix=1.0_dp, delta=tmp, para_env=para_env, &
                           m3=matrix_q(ispin, ic)%matrix)
            delta = MAX(tmp, delta)
         END DO
      END DO
      CALL timestop(handle)

   END SUBROUTINE self_consistency_check

! **************************************************************************************************
!> \brief  allocation needed when density mixing is used
!> \param qs_env ...
!> \param mixing_method ...
!> \param p_mix_new ...
!> \param p_delta ...
!> \param nspins ...
!> \param mixing_store ...
!> \par History
!>      05.2009 created [MI]
!>      08.2014 kpoints [JGH]
!>      02.2015 changed input to make g-space mixing available in linear scaling SCF [Patrick Seewald]
!>      02.2019 charge mixing [JGH]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE mixing_allocate(qs_env, mixing_method, p_mix_new, p_delta, nspins, mixing_store)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: mixing_method
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: p_mix_new, p_delta
      INTEGER, INTENT(IN)                                :: nspins
      TYPE(mixing_storage_type), POINTER                 :: mixing_store

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mixing_allocate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ia, iat, ic, ikind, ispin, &
                                                            max_shell, na, natom, nbuffer, nel, &
                                                            nimg, nkind
      LOGICAL                                            :: charge_mixing
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: refmatrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom

      CALL timeset(routineN, handle)

      NULLIFY (matrix_s, dft_control, sab_orb, refmatrix, rho_atom)
      CALL get_qs_env(qs_env, &
                      sab_orb=sab_orb, &
                      matrix_s_kp=matrix_s, &
                      dft_control=dft_control)

      charge_mixing = dft_control%qs_control%dftb .OR. dft_control%qs_control%xtb &
                      .OR. dft_control%qs_control%semi_empirical
      refmatrix => matrix_s(1, 1)%matrix
      nimg = dft_control%nimages

      !   *** allocate p_mix_new ***
      IF (PRESENT(p_mix_new)) THEN
         IF (.NOT. ASSOCIATED(p_mix_new)) THEN
            CALL dbcsr_allocate_matrix_set(p_mix_new, nspins, nimg)
            DO i = 1, nspins
               DO ic = 1, nimg
                  ALLOCATE (p_mix_new(i, ic)%matrix)
                  CALL dbcsr_create(matrix=p_mix_new(i, ic)%matrix, template=refmatrix, &
                                    name="SCF DENSITY", matrix_type=dbcsr_type_symmetric, nze=0)
                  CALL cp_dbcsr_alloc_block_from_nbl(p_mix_new(i, ic)%matrix, sab_orb)
                  CALL dbcsr_set(p_mix_new(i, ic)%matrix, 0.0_dp)
               ENDDO
            ENDDO
         END IF
      END IF

      !   *** allocate p_delta ***
      IF (PRESENT(p_delta)) THEN
         IF (mixing_method >= gspace_mixing_nr) THEN
            IF (.NOT. ASSOCIATED(p_delta)) THEN
               CALL dbcsr_allocate_matrix_set(p_delta, nspins, nimg)
               DO i = 1, nspins
                  DO ic = 1, nimg
                     ALLOCATE (p_delta(i, ic)%matrix)
                     CALL dbcsr_create(matrix=p_delta(i, ic)%matrix, template=refmatrix, &
                                       name="SCF DENSITY", matrix_type=dbcsr_type_symmetric, nze=0)
                     CALL cp_dbcsr_alloc_block_from_nbl(p_delta(i, ic)%matrix, sab_orb)
                     CALL dbcsr_set(p_delta(i, ic)%matrix, 0.0_dp)
                  ENDDO
               ENDDO
            END IF
            CPASSERT(ASSOCIATED(mixing_store))
         END IF
      END IF

      IF (charge_mixing) THEN
         !   *** allocate buffer for charge mixing ***
         IF (mixing_method >= gspace_mixing_nr) THEN
            CPASSERT(.NOT. mixing_store%gmix_p)
            IF (dft_control%qs_control%dftb) THEN
               max_shell = 1
            ELSEIF (dft_control%qs_control%xtb) THEN
               max_shell = 5
            ELSE
               CPABORT('UNKNOWN METHOD')
            END IF
            nbuffer = mixing_store%nbuffer
            mixing_store%ncall = 0
            CALL get_qs_env(qs_env, local_particles=distribution_1d)
            nkind = SIZE(distribution_1d%n_el)
            na = SUM(distribution_1d%n_el(:))
            IF (ASSOCIATED(mixing_store%atlist)) DEALLOCATE (mixing_store%atlist)
            ALLOCATE (mixing_store%atlist(na))
            mixing_store%nat_local = na
            mixing_store%max_shell = max_shell
            ia = 0
            DO ikind = 1, nkind
               nel = distribution_1d%n_el(ikind)
               DO iat = 1, nel
                  ia = ia + 1
                  mixing_store%atlist(ia) = distribution_1d%list(ikind)%array(iat)
               END DO
            END DO
            IF (ASSOCIATED(mixing_store%acharge)) DEALLOCATE (mixing_store%acharge)
            ALLOCATE (mixing_store%acharge(na, max_shell, nbuffer))
            IF (ASSOCIATED(mixing_store%dacharge)) DEALLOCATE (mixing_store%dacharge)
            ALLOCATE (mixing_store%dacharge(na, max_shell, nbuffer))
         END IF
         IF (mixing_method == pulay_mixing_nr) THEN
            IF (ASSOCIATED(mixing_store%pulay_matrix)) DEALLOCATE (mixing_store%pulay_matrix)
            ALLOCATE (mixing_store%pulay_matrix(nbuffer, nbuffer))
            mixing_store%pulay_matrix = 0.0_dp
         ELSEIF (mixing_method == broyden_mixing_nr) THEN
            IF (ASSOCIATED(mixing_store%abroy)) DEALLOCATE (mixing_store%abroy)
            ALLOCATE (mixing_store%abroy(nbuffer, nbuffer))
            mixing_store%abroy = 0.0_dp
            IF (ASSOCIATED(mixing_store%wbroy)) DEALLOCATE (mixing_store%wbroy)
            ALLOCATE (mixing_store%wbroy(nbuffer))
            mixing_store%wbroy = 0.0_dp
            IF (ASSOCIATED(mixing_store%dfbroy)) DEALLOCATE (mixing_store%dfbroy)
            ALLOCATE (mixing_store%dfbroy(na, max_shell, nbuffer))
            mixing_store%dfbroy = 0.0_dp
            IF (ASSOCIATED(mixing_store%ubroy)) DEALLOCATE (mixing_store%ubroy)
            ALLOCATE (mixing_store%ubroy(na, max_shell, nbuffer))
            mixing_store%ubroy = 0.0_dp
         ELSEIF (mixing_method == broyden_mixing_new_nr) THEN
            CPABORT("broyden_mixing_new not available")
         ELSEIF (mixing_method == multisecant_mixing_nr) THEN
            CPABORT("multisecant_mixing not available")
         END IF
      ELSE
         !   *** allocate buffer for gspace mixing ***
         IF (mixing_method >= gspace_mixing_nr) THEN
            nbuffer = mixing_store%nbuffer
            mixing_store%ncall = 0
            IF (.NOT. ASSOCIATED(mixing_store%rhoin)) THEN
               ALLOCATE (mixing_store%rhoin(nspins))
               DO ispin = 1, nspins
                  NULLIFY (mixing_store%rhoin(ispin)%cc)
               END DO
            END IF
            IF (mixing_store%gmix_p) THEN
               IF (.NOT. ASSOCIATED(mixing_store%rho_ao_in)) THEN
                  CALL dbcsr_allocate_matrix_set(mixing_store%rho_ao_in, nspins, nimg)
                  DO i = 1, nspins
                     DO ic = 1, nimg
                        ALLOCATE (mixing_store%rho_ao_in(i, ic)%matrix)
                        CALL dbcsr_create(matrix=mixing_store%rho_ao_in(i, ic)%matrix, template=refmatrix, &
                                          name="SCF DENSITY", matrix_type=dbcsr_type_symmetric, nze=0)
                        CALL cp_dbcsr_alloc_block_from_nbl(mixing_store%rho_ao_in(i, ic)%matrix, sab_orb)
                        CALL dbcsr_set(mixing_store%rho_ao_in(i, ic)%matrix, 0.0_dp)
                     ENDDO
                  END DO
               END IF
               IF (.NOT. ASSOCIATED(mixing_store%rho_ao_mix)) THEN
                  ALLOCATE (mixing_store%rho_ao_mix)
                  CALL dbcsr_create(matrix=mixing_store%rho_ao_mix, template=refmatrix, &
                                    name="DENSITY MIX", matrix_type=dbcsr_type_symmetric, nze=0)
                  CALL cp_dbcsr_alloc_block_from_nbl(mixing_store%rho_ao_mix, sab_orb)
                  CALL dbcsr_set(mixing_store%rho_ao_mix, 0.0_dp)
               END IF
            END IF
            IF (mixing_store%gmix_p .AND. dft_control%qs_control%gapw) THEN
               CALL get_qs_env(qs_env=qs_env, rho_atom_set=rho_atom)
               natom = SIZE(rho_atom)
               IF (.NOT. ASSOCIATED(mixing_store%paw)) THEN
                  ALLOCATE (mixing_store%paw(natom))
                  mixing_store%paw = .FALSE.
                  ALLOCATE (mixing_store%cpc_h_in(natom, nspins))
                  ALLOCATE (mixing_store%cpc_s_in(natom, nspins))
                  DO ispin = 1, nspins
                     DO iat = 1, natom
                        NULLIFY (mixing_store%cpc_h_in(iat, ispin)%r_coef)
                        NULLIFY (mixing_store%cpc_s_in(iat, ispin)%r_coef)
                     END DO
                  END DO
               END IF
            END IF
         END IF

         !   *** allocare res_buffer if needed
         IF (mixing_method >= pulay_mixing_nr) THEN
            IF (.NOT. ASSOCIATED(mixing_store%res_buffer)) THEN
               ALLOCATE (mixing_store%res_buffer(nbuffer, nspins))
               DO ispin = 1, nspins
                  DO i = 1, nbuffer
                     NULLIFY (mixing_store%res_buffer(i, ispin)%cc)
                  END DO
               END DO
            END IF
         END IF

         !   *** allocate pulay
         IF (mixing_method == pulay_mixing_nr) THEN
            IF (.NOT. ASSOCIATED(mixing_store%pulay_matrix)) THEN
               ALLOCATE (mixing_store%pulay_matrix(nbuffer, nbuffer))
            END IF

            IF (.NOT. ASSOCIATED(mixing_store%rhoin_buffer)) THEN
               ALLOCATE (mixing_store%rhoin_buffer(nbuffer, nspins))
               DO ispin = 1, nspins
                  DO i = 1, nbuffer
                     NULLIFY (mixing_store%rhoin_buffer(i, ispin)%cc)
                  END DO
               END DO
            END IF
            IF (mixing_store%gmix_p) THEN
               IF (.NOT. ASSOCIATED(mixing_store%rho_ao_in_buffer)) THEN
                  CALL dbcsr_allocate_matrix_set(mixing_store%rho_ao_in_buffer, nspins, nimg, nbuffer)
                  DO ispin = 1, nspins
                     DO ic = 1, nimg
                        DO i = 1, nbuffer
                           ALLOCATE (mixing_store%rho_ao_in_buffer(ispin, ic, i)%matrix)
                           CALL dbcsr_create(matrix=mixing_store%rho_ao_in_buffer(ispin, ic, i)%matrix, &
                                             template=refmatrix, name="SCF DENSITY", &
                                             matrix_type=dbcsr_type_symmetric, nze=0)
                           CALL cp_dbcsr_alloc_block_from_nbl(mixing_store%rho_ao_in_buffer(ispin, ic, i)%matrix, sab_orb)
                           CALL dbcsr_set(mixing_store%rho_ao_in_buffer(ispin, ic, i)%matrix, 0.0_dp)
                        END DO
                     ENDDO
                  END DO
               END IF
               IF (dft_control%qs_control%gapw) THEN
                  IF (.NOT. ASSOCIATED(mixing_store%cpc_h_in_buffer)) THEN
                     ALLOCATE (mixing_store%cpc_h_in_buffer(nbuffer, natom, nspins))
                     ALLOCATE (mixing_store%cpc_s_in_buffer(nbuffer, natom, nspins))
                     ALLOCATE (mixing_store%cpc_h_res_buffer(nbuffer, natom, nspins))
                     ALLOCATE (mixing_store%cpc_s_res_buffer(nbuffer, natom, nspins))
                     DO ispin = 1, nspins
                        DO iat = 1, natom
                           DO i = 1, nbuffer
                              NULLIFY (mixing_store%cpc_h_in_buffer(i, iat, ispin)%r_coef)
                              NULLIFY (mixing_store%cpc_s_in_buffer(i, iat, ispin)%r_coef)
                              NULLIFY (mixing_store%cpc_h_res_buffer(i, iat, ispin)%r_coef)
                              NULLIFY (mixing_store%cpc_s_res_buffer(i, iat, ispin)%r_coef)
                           END DO
                        END DO
                     END DO
                  END IF
               END IF
            END IF

         END IF
         !   *** allocate broyden buffer ***
         IF (mixing_method == broyden_mixing_nr) THEN
            IF (.NOT. ASSOCIATED(mixing_store%rhoin_old)) THEN
               ALLOCATE (mixing_store%rhoin_old(nspins))
               DO ispin = 1, nspins
                  NULLIFY (mixing_store%rhoin_old(ispin)%cc)
               END DO
            END IF
            IF (.NOT. ASSOCIATED(mixing_store%drho_buffer)) THEN
               ALLOCATE (mixing_store%drho_buffer(nbuffer, nspins))
               ALLOCATE (mixing_store%last_res(nspins))
               DO ispin = 1, nspins
                  DO i = 1, nbuffer
                     NULLIFY (mixing_store%drho_buffer(i, ispin)%cc)
                  END DO
                  NULLIFY (mixing_store%last_res(ispin)%cc)
               END DO
            END IF
            IF (mixing_store%gmix_p) THEN
               IF (.NOT. ASSOCIATED(mixing_store%rho_ao_res)) THEN
                  ALLOCATE (mixing_store%rho_ao_res)
                  CALL dbcsr_create(matrix=mixing_store%rho_ao_res, template=refmatrix, &
                                    name="DENSITY MIX", matrix_type=dbcsr_type_symmetric, nze=0)
                  CALL cp_dbcsr_alloc_block_from_nbl(mixing_store%rho_ao_res, sab_orb)
                  CALL dbcsr_set(mixing_store%rho_ao_res, 0.0_dp)
               END IF

               IF (.NOT. ASSOCIATED(mixing_store%rho_ao_in_old)) THEN
                  CALL dbcsr_allocate_matrix_set(mixing_store%rho_ao_in_old, nspins, nimg)
                  DO i = 1, nspins
                     DO ic = 1, nimg
                        ALLOCATE (mixing_store%rho_ao_in_old(i, ic)%matrix)
                        CALL dbcsr_create(matrix=mixing_store%rho_ao_in_old(i, ic)%matrix, template=refmatrix, &
                                          name="SCF DENSITY", matrix_type=dbcsr_type_symmetric, nze=0)
                        CALL cp_dbcsr_alloc_block_from_nbl(mixing_store%rho_ao_in_old(i, ic)%matrix, sab_orb)
                        CALL dbcsr_set(mixing_store%rho_ao_in_old(i, ic)%matrix, 0.0_dp)
                     ENDDO
                  END DO
               END IF

               IF (.NOT. ASSOCIATED(mixing_store%rho_ao_lastres)) THEN
                  CALL dbcsr_allocate_matrix_set(mixing_store%rho_ao_lastres, nspins, nimg)
                  DO i = 1, nspins
                     DO ic = 1, nimg
                        ALLOCATE (mixing_store%rho_ao_lastres(i, ic)%matrix)
                        CALL dbcsr_create(matrix=mixing_store%rho_ao_lastres(i, ic)%matrix, template=refmatrix, &
                                          name="SCF DENSITY", matrix_type=dbcsr_type_symmetric, nze=0)
                        CALL cp_dbcsr_alloc_block_from_nbl(mixing_store%rho_ao_lastres(i, ic)%matrix, sab_orb)
                        CALL dbcsr_set(mixing_store%rho_ao_lastres(i, ic)%matrix, 0.0_dp)
                     ENDDO
                  END DO
               END IF

               IF (dft_control%qs_control%gapw) THEN
                  IF (.NOT. ASSOCIATED(mixing_store%cpc_h_old)) THEN
                     ALLOCATE (mixing_store%cpc_h_old(natom, nspins))
                     ALLOCATE (mixing_store%cpc_s_old(natom, nspins))
                     DO ispin = 1, nspins
                        DO iat = 1, natom
                           NULLIFY (mixing_store%cpc_h_old(iat, ispin)%r_coef)
                           NULLIFY (mixing_store%cpc_s_old(iat, ispin)%r_coef)
                        END DO
                     END DO
                  END IF
                  IF (.NOT. ASSOCIATED(mixing_store%dcpc_h_in)) THEN
                     ALLOCATE (mixing_store%dcpc_h_in(nbuffer, natom, nspins))
                     ALLOCATE (mixing_store%dcpc_s_in(nbuffer, natom, nspins))
                     ALLOCATE (mixing_store%cpc_h_lastres(natom, nspins))
                     ALLOCATE (mixing_store%cpc_s_lastres(natom, nspins))
                     DO ispin = 1, nspins
                        DO iat = 1, natom
                           DO i = 1, nbuffer
                              NULLIFY (mixing_store%dcpc_h_in(i, iat, ispin)%r_coef)
                              NULLIFY (mixing_store%dcpc_s_in(i, iat, ispin)%r_coef)
                           END DO
                           NULLIFY (mixing_store%cpc_h_lastres(iat, ispin)%r_coef)
                           NULLIFY (mixing_store%cpc_s_lastres(iat, ispin)%r_coef)
                        END DO
                     END DO
                  END IF
               END IF
            END IF
         END IF

         IF (mixing_method == pulay_mixing_nr .OR. mixing_method == broyden_mixing_nr) THEN
            IF (mixing_store%gmix_p) THEN
               IF (.NOT. ASSOCIATED(mixing_store%rho_ao_res_buffer)) THEN
                  CALL dbcsr_allocate_matrix_set(mixing_store%rho_ao_res_buffer, nspins, nimg, nbuffer)
                  DO ispin = 1, nspins
                     DO ic = 1, nimg
                        DO i = 1, nbuffer
                           ALLOCATE (mixing_store%rho_ao_res_buffer(ispin, ic, i)%matrix)
                           CALL dbcsr_create(matrix=mixing_store%rho_ao_res_buffer(ispin, ic, i)%matrix, &
                                             template=refmatrix, name="SCF DENSITY", &
                                             matrix_type=dbcsr_type_symmetric, nze=0)
                           CALL cp_dbcsr_alloc_block_from_nbl(mixing_store%rho_ao_res_buffer(ispin, ic, i)%matrix, sab_orb)
                           CALL dbcsr_set(mixing_store%rho_ao_res_buffer(ispin, ic, i)%matrix, 0.0_dp)
                        END DO
                     ENDDO
                  END DO
               END IF
            END IF
         END IF
         !   *** allocate broyden buffer ***
         IF (mixing_method == broyden_mixing_new_nr) THEN
            IF (.NOT. ASSOCIATED(mixing_store%u_vec)) THEN
               ALLOCATE (mixing_store%last_res(nspins))
               ALLOCATE (mixing_store%delta_res(nbuffer - 1, nspins))
               ALLOCATE (mixing_store%u_vec(nbuffer - 1, nspins))
               ALLOCATE (mixing_store%z_vec(nbuffer - 1, nspins))
               ALLOCATE (mixing_store%weight(nbuffer, nspins))
               ALLOCATE (mixing_store%fmat(nbuffer - 1, nbuffer - 1, nspins))
               DO ispin = 1, nspins
                  DO i = 1, nbuffer - 1
                     NULLIFY (mixing_store%delta_res(i, ispin)%cc)
                     NULLIFY (mixing_store%u_vec(i, ispin)%cc)
                     NULLIFY (mixing_store%z_vec(i, ispin)%cc)
                  END DO
                  NULLIFY (mixing_store%last_res(ispin)%cc)
               END DO
            END IF
         END IF

         !   *** allocate multisecant buffer ***
         IF (mixing_method == multisecant_mixing_nr) THEN
            IF (.NOT. ASSOCIATED(mixing_store%norm_res_buffer)) THEN
               ALLOCATE (mixing_store%norm_res_buffer(nbuffer, nspins))
            END IF
         END IF

         IF (mixing_method == broyden_mixing_new_nr .OR. mixing_method == multisecant_mixing_nr) THEN
            IF (.NOT. ASSOCIATED(mixing_store%rhoin_buffer)) THEN
               ALLOCATE (mixing_store%rhoin_buffer(nbuffer, nspins))
               DO ispin = 1, nspins
                  DO i = 1, nbuffer
                     NULLIFY (mixing_store%rhoin_buffer(i, ispin)%cc)
                  END DO
               END DO
            END IF
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE mixing_allocate

! **************************************************************************************************
!> \brief  initialiation needed when gspace mixing is used
!> \param mixing_method ...
!> \param rho ...
!> \param mixing_store ...
!> \param para_env ...
!> \param rho_atom ...
!> \par History
!>      05.2009 created [MI]
!> \author MI
! **************************************************************************************************
   SUBROUTINE mixing_init(mixing_method, rho, mixing_store, para_env, rho_atom)
      INTEGER, INTENT(IN)                                :: mixing_method
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(mixing_storage_type), POINTER                 :: mixing_store
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(rho_atom_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: rho_atom

      CHARACTER(len=*), PARAMETER :: routineN = 'mixing_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iat, ib, ig, ig1, ig_count, &
                                                            iproc, ispin, n1, n2, natom, nbuffer, &
                                                            ng, nimg, nspin
      REAL(dp)                                           :: bconst, beta, fdamp, g2max, g2min, kmin
      REAL(dp), DIMENSION(:), POINTER                    :: g2
      REAL(dp), DIMENSION(:, :), POINTER                 :: g_vec
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g

      CALL timeset(routineN, handle)

      NULLIFY (g2, g_vec, rho_ao_kp, rho_g)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp, rho_g=rho_g)

      nspin = SIZE(rho_g)
      ng = SIZE(rho_g(1)%pw%pw_grid%gsq, 1)
      nimg = SIZE(rho_ao_kp, 2)
      mixing_store%ig_max = ng
      g2 => rho_g(1)%pw%pw_grid%gsq
      g_vec => rho_g(1)%pw%pw_grid%g

      IF (mixing_store%max_gvec_exp > 0._dp) THEN
         DO ig = 1, ng
            IF (g2(ig) > mixing_store%max_g2) THEN
               mixing_store%ig_max = ig
               EXIT
            END IF
         END DO
      END IF

      IF (.NOT. ASSOCIATED(mixing_store%kerker_factor)) THEN
         ALLOCATE (mixing_store%kerker_factor(ng))
      END IF
      IF (.NOT. ASSOCIATED(mixing_store%special_metric)) THEN
         ALLOCATE (mixing_store%special_metric(ng))
      END IF
      beta = mixing_store%beta
      kmin = 0.1_dp
      mixing_store%kerker_factor = 1.0_dp
      mixing_store%special_metric = 1.0_dp
      ig1 = 1
      IF (rho_g(1)%pw%pw_grid%have_g0) ig1 = 2
      DO ig = ig1, mixing_store%ig_max
         mixing_store%kerker_factor(ig) = MAX(g2(ig)/(g2(ig) + beta*beta), kmin)
         mixing_store%special_metric(ig) = &
            1.0_dp + 50.0_dp/8.0_dp*( &
            1.0_dp + COS(g_vec(1, ig)) + COS(g_vec(2, ig)) + COS(g_vec(3, ig)) + &
            COS(g_vec(1, ig))*COS(g_vec(2, ig)) + &
            COS(g_vec(2, ig))*COS(g_vec(3, ig)) + &
            COS(g_vec(1, ig))*COS(g_vec(3, ig)) + &
            COS(g_vec(1, ig))*COS(g_vec(2, ig))*COS(g_vec(3, ig)))
      END DO

      nbuffer = mixing_store%nbuffer
      DO ispin = 1, nspin
         IF (.NOT. ASSOCIATED(mixing_store%rhoin(ispin)%cc)) THEN
            ALLOCATE (mixing_store%rhoin(ispin)%cc(ng))
         END IF
         mixing_store%rhoin(ispin)%cc = rho_g(ispin)%pw%cc

         IF (ASSOCIATED(mixing_store%rhoin_buffer)) THEN
            IF (.NOT. ASSOCIATED(mixing_store%rhoin_buffer(1, ispin)%cc)) THEN
               DO ib = 1, nbuffer
                  ALLOCATE (mixing_store%rhoin_buffer(ib, ispin)%cc(ng))
               END DO
            END IF
            mixing_store%rhoin_buffer(1, ispin)%cc(1:ng) = &
               rho_g(ispin)%pw%cc(1:ng)
         END IF
         IF (ASSOCIATED(mixing_store%res_buffer)) THEN
            IF (.NOT. ASSOCIATED(mixing_store%res_buffer(1, ispin)%cc)) THEN
               DO ib = 1, nbuffer
                  ALLOCATE (mixing_store%res_buffer(ib, ispin)%cc(ng))
               END DO
            END IF
         END IF
      END DO

      IF (nspin == 2) THEN
         mixing_store%rhoin(1)%cc = rho_g(1)%pw%cc + rho_g(2)%pw%cc
         mixing_store%rhoin(2)%cc = rho_g(1)%pw%cc - rho_g(2)%pw%cc
         IF (ASSOCIATED(mixing_store%rhoin_buffer)) THEN
            mixing_store%rhoin_buffer(1, 1)%cc = rho_g(1)%pw%cc + rho_g(2)%pw%cc
            mixing_store%rhoin_buffer(1, 2)%cc = rho_g(1)%pw%cc - rho_g(2)%pw%cc
         END IF
      END IF

      IF (mixing_store%gmix_p) THEN
         IF (ASSOCIATED(mixing_store%rho_ao_in)) THEN
            DO ispin = 1, nspin
               DO ib = 1, nimg
                  CALL dbcsr_copy(mixing_store%rho_ao_in(ispin, ib)%matrix, rho_ao_kp(ispin, ib)%matrix)
               END DO
            END DO
         END IF
         IF (PRESENT(rho_atom)) THEN
            natom = SIZE(rho_atom)
            DO ispin = 1, nspin
               DO iat = 1, natom
                  IF (ASSOCIATED(rho_atom(iat)%cpc_s(ispin)%r_coef)) THEN
                     mixing_store%paw(iat) = .TRUE.
                     n1 = SIZE(rho_atom(iat)%cpc_s(ispin)%r_coef, 1)
                     n2 = SIZE(rho_atom(iat)%cpc_s(ispin)%r_coef, 2)
                     IF (ASSOCIATED(mixing_store%cpc_s_in)) THEN
                        IF (.NOT. ASSOCIATED(mixing_store%cpc_s_in(iat, ispin)%r_coef)) THEN
                           ALLOCATE (mixing_store%cpc_s_in(iat, ispin)%r_coef(n1, n2))
                           ALLOCATE (mixing_store%cpc_h_in(iat, ispin)%r_coef(n1, n2))
                        END IF
                        mixing_store%cpc_h_in(iat, ispin)%r_coef = rho_atom(iat)%cpc_h(ispin)%r_coef
                        mixing_store%cpc_s_in(iat, ispin)%r_coef = rho_atom(iat)%cpc_s(ispin)%r_coef
                     END IF
                  END IF
               END DO
            END DO
         END IF
      END IF

      IF (mixing_method == gspace_mixing_nr) THEN
      ELSEIF (mixing_method == pulay_mixing_nr) THEN
         IF (mixing_store%gmix_p .AND. PRESENT(rho_atom)) THEN
            DO ispin = 1, nspin
               DO iat = 1, natom
                  IF (mixing_store%paw(iat)) THEN
                     n1 = SIZE(rho_atom(iat)%cpc_s(ispin)%r_coef, 1)
                     n2 = SIZE(rho_atom(iat)%cpc_s(ispin)%r_coef, 2)
                     IF (.NOT. ASSOCIATED(mixing_store%cpc_h_in_buffer(1, iat, ispin)%r_coef)) THEN
                        DO ib = 1, nbuffer
                           ALLOCATE (mixing_store%cpc_s_in_buffer(ib, iat, ispin)%r_coef(n1, n2))
                           ALLOCATE (mixing_store%cpc_h_in_buffer(ib, iat, ispin)%r_coef(n1, n2))
                           ALLOCATE (mixing_store%cpc_s_res_buffer(ib, iat, ispin)%r_coef(n1, n2))
                           ALLOCATE (mixing_store%cpc_h_res_buffer(ib, iat, ispin)%r_coef(n1, n2))
                        END DO
                     END IF
                     DO ib = 1, nbuffer
                        mixing_store%cpc_h_in_buffer(ib, iat, ispin)%r_coef = 0.0_dp
                        mixing_store%cpc_s_in_buffer(ib, iat, ispin)%r_coef = 0.0_dp
                        mixing_store%cpc_h_res_buffer(ib, iat, ispin)%r_coef = 0.0_dp
                        mixing_store%cpc_s_res_buffer(ib, iat, ispin)%r_coef = 0.0_dp
                     END DO
                  END IF
               END DO
            END DO
         END IF
      ELSEIF (mixing_method == broyden_mixing_nr) THEN
         DO ispin = 1, nspin
            IF (.NOT. ASSOCIATED(mixing_store%rhoin_old(ispin)%cc)) THEN
               ALLOCATE (mixing_store%rhoin_old(ispin)%cc(ng))
            END IF
            IF (.NOT. ASSOCIATED(mixing_store%drho_buffer(1, ispin)%cc)) THEN
               DO ib = 1, nbuffer
                  ALLOCATE (mixing_store%drho_buffer(ib, ispin)%cc(ng))
               END DO
               ALLOCATE (mixing_store%last_res(ispin)%cc(ng))
            END IF
            DO ib = 1, nbuffer
               mixing_store%drho_buffer(ib, ispin)%cc = CMPLX(0.0_dp, 0.0_dp, kind=dp)
            END DO
            mixing_store%last_res(ispin)%cc = CMPLX(0.0_dp, 0.0_dp, kind=dp)
            mixing_store%rhoin_old(ispin)%cc = CMPLX(0.0_dp, 0.0_dp, kind=dp)
         END DO
         IF (mixing_store%gmix_p) THEN
            IF (PRESENT(rho_atom)) THEN
               DO ispin = 1, nspin
                  DO iat = 1, natom
                     IF (mixing_store%paw(iat)) THEN
                        n1 = SIZE(rho_atom(iat)%cpc_s(ispin)%r_coef, 1)
                        n2 = SIZE(rho_atom(iat)%cpc_s(ispin)%r_coef, 2)
                        IF (.NOT. ASSOCIATED(mixing_store%cpc_s_old(iat, ispin)%r_coef)) THEN
                           ALLOCATE (mixing_store%cpc_s_old(iat, ispin)%r_coef(n1, n2))
                           ALLOCATE (mixing_store%cpc_h_old(iat, ispin)%r_coef(n1, n2))
                        END IF
                        mixing_store%cpc_h_old(iat, ispin)%r_coef = 0.0_dp
                        mixing_store%cpc_s_old(iat, ispin)%r_coef = 0.0_dp
                        IF (.NOT. ASSOCIATED(mixing_store%dcpc_s_in(1, iat, ispin)%r_coef)) THEN
                           DO ib = 1, nbuffer
                              ALLOCATE (mixing_store%dcpc_h_in(ib, iat, ispin)%r_coef(n1, n2))
                              ALLOCATE (mixing_store%dcpc_s_in(ib, iat, ispin)%r_coef(n1, n2))
                           END DO
                           ALLOCATE (mixing_store%cpc_h_lastres(iat, ispin)%r_coef(n1, n2))
                           ALLOCATE (mixing_store%cpc_s_lastres(iat, ispin)%r_coef(n1, n2))
                        END IF
                        DO ib = 1, nbuffer
                           mixing_store%dcpc_h_in(ib, iat, ispin)%r_coef = 0.0_dp
                           mixing_store%dcpc_s_in(ib, iat, ispin)%r_coef = 0.0_dp
                        END DO
                        mixing_store%cpc_h_lastres(iat, ispin)%r_coef = 0.0_dp
                        mixing_store%cpc_s_lastres(iat, ispin)%r_coef = 0.0_dp
                     END IF
                  END DO
               END DO
            END IF
         END IF

         IF (.NOT. ASSOCIATED(mixing_store%p_metric)) THEN
            ALLOCATE (mixing_store%p_metric(ng))
            bconst = mixing_store%bconst
            g2min = 1.E30_dp
            DO ig = 1, ng
               IF (g2(ig) > 1.E-10_dp) g2min = MIN(g2min, g2(ig))
            END DO
            g2max = -1.E30_dp
            DO ig = 1, ng
               g2max = MAX(g2max, g2(ig))
            END DO
            CALL mp_min(g2min, para_env%group)
            CALL mp_max(g2max, para_env%group)
            ! fdamp/g2 varies between (bconst-1) and 0
            ! i.e. p_metric varies between bconst and 1
            ! fdamp = (bconst-1.0_dp)*g2min
            fdamp = (bconst - 1.0_dp)*g2min*g2max/(g2max - g2min*bconst)
            DO ig = 1, ng
               mixing_store%p_metric(ig) = (g2(ig) + fdamp)/MAX(g2(ig), 1.E-10_dp)
            END DO
            IF (rho_g(1)%pw%pw_grid%have_g0) mixing_store%p_metric(1) = bconst
         END IF

      ELSEIF (mixing_method == broyden_mixing_new_nr) THEN
         DO ispin = 1, nspin
            IF (.NOT. ASSOCIATED(mixing_store%u_vec(1, ispin)%cc)) THEN
               DO ib = 1, nbuffer - 1
                  ALLOCATE (mixing_store%delta_res(ib, ispin)%cc(ng))
                  ALLOCATE (mixing_store%u_vec(ib, ispin)%cc(ng))
                  ALLOCATE (mixing_store%z_vec(ib, ispin)%cc(ng))
               END DO
               ALLOCATE (mixing_store%last_res(ispin)%cc(ng))
            END IF
         END DO
         IF (.NOT. ASSOCIATED(mixing_store%p_metric)) THEN
            ALLOCATE (mixing_store%p_metric(ng))
            bconst = mixing_store%bconst
            g2min = 1.E30_dp
            DO ig = 1, ng
               IF (g2(ig) > 1.E-10_dp) g2min = MIN(g2min, g2(ig))
            END DO
            CALL mp_min(g2min, para_env%group)
            fdamp = (bconst - 1.0_dp)*g2min
            DO ig = 1, ng
               mixing_store%p_metric(ig) = (g2(ig) + fdamp)/MAX(g2(ig), 1.E-10_dp)
            END DO
            IF (rho_g(1)%pw%pw_grid%have_g0) mixing_store%p_metric(1) = bconst
         END IF
      ELSEIF (mixing_method == multisecant_mixing_nr) THEN
         IF (.NOT. ASSOCIATED(mixing_store%ig_global_index)) THEN
            ALLOCATE (mixing_store%ig_global_index(ng))
         END IF
         mixing_store%ig_global_index = 0
         ig_count = 0
         DO iproc = 0, para_env%num_pe - 1
            IF (para_env%mepos == iproc) THEN
               DO ig = 1, ng
                  ig_count = ig_count + 1
                  mixing_store%ig_global_index(ig) = ig_count
               END DO
            END IF
            CALL mp_bcast(ig_count, iproc, para_env%group)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE mixing_init

! **************************************************************************************************
!> \brief initialiation needed when charge mixing is used
!> \param mixing_store ...
!> \par History
!>      02.2019 created [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE charge_mixing_init(mixing_store)
      TYPE(mixing_storage_type), POINTER                 :: mixing_store

      CHARACTER(len=*), PARAMETER :: routineN = 'charge_mixing_init', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      mixing_store%acharge = 0.0_dp

      CALL timestop(handle)

   END SUBROUTINE charge_mixing_init

END MODULE qs_mixing_utils
