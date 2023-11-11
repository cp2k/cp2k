!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief methods of the rho structure (defined in qs_rho_types)
!> \par History
!>      08.2002 created [fawzi]
!>      08.2014 kpoints [JGH]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE qs_rho_methods
   USE admm_types,                      ONLY: get_admm_env
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_to_string
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_p_type, dbcsr_scale, dbcsr_set, dbcsr_type, &
        dbcsr_type_antisymmetric, dbcsr_type_symmetric
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE lri_environment_methods,         ONLY: calculate_lri_densities
   USE lri_environment_types,           ONLY: lri_density_type,&
                                              lri_environment_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_scale,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
   USE qs_collocate_density,            ONLY: calculate_drho_elec,&
                                              calculate_rho_elec
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_local_rho_types,              ONLY: local_rho_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_oce_types,                    ONLY: oce_matrix_type
   USE qs_rho_atom_methods,             ONLY: calculate_rho_atom_coeff
   USE qs_rho_atom_types,               ONLY: rho_atom_type
   USE qs_rho_types,                    ONLY: qs_rho_clear,&
                                              qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE ri_environment_methods,          ONLY: calculate_ri_densities
   USE task_list_types,                 ONLY: task_list_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_rho_methods'

   PUBLIC :: qs_rho_update_rho, qs_rho_update_tddfpt, &
             qs_rho_rebuild, qs_rho_copy, qs_rho_scale_and_add
   PUBLIC :: duplicate_rho_type, allocate_rho_ao_imag_from_real

CONTAINS

! **************************************************************************************************
!> \brief rebuilds rho (if necessary allocating and initializing it)
!> \param rho the rho type to rebuild (defaults to qs_env%rho)
!> \param qs_env the environment to which rho belongs
!> \param rebuild_ao if it is necessary to rebuild rho_ao. Defaults to true.
!> \param rebuild_grids if it in necessary to rebuild rho_r and rho_g.
!>        Defaults to false.
!> \param admm (use aux_fit basis)
!> \param pw_env_external external plane wave environment
!> \par History
!>      11.2002 created replacing qs_rho_create and qs_env_rebuild_rho[fawzi]
!> \author Fawzi Mohamed
!> \note
!>      needs updated  pw pools, s, s_mstruct and h in qs_env.
!>      The use of p to keep the structure of h (needed for the forces)
!>      is ugly and should be removed.
!>      Change so that it does not allocate a subcomponent if it is not
!>      associated and not requested?
! **************************************************************************************************
   SUBROUTINE qs_rho_rebuild(rho, qs_env, rebuild_ao, rebuild_grids, admm, pw_env_external)
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(in), OPTIONAL                      :: rebuild_ao, rebuild_grids, admm
      TYPE(pw_env_type), OPTIONAL, POINTER               :: pw_env_external

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'qs_rho_rebuild'

      CHARACTER(LEN=default_string_length)               :: headline
      INTEGER                                            :: handle, i, ic, j, nimg, nspins
      LOGICAL                                            :: do_kpoints, my_admm, my_rebuild_ao, &
                                                            my_rebuild_grids, rho_ao_is_complex
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, rho_ao_im_kp, rho_ao_kp
      TYPE(dbcsr_type), POINTER                          :: refmatrix, tmatrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type), DIMENSION(:), POINTER               :: rho_g, rho_r, tau_g, tau_r
      TYPE(pw_type), DIMENSION(:, :), POINTER            :: drho_g, drho_r
      TYPE(pw_type), POINTER                             :: rho_r_sccs

      CALL timeset(routineN, handle)

      NULLIFY (pw_env, auxbas_pw_pool, matrix_s_kp, dft_control)
      NULLIFY (tot_rho_r, rho_ao_kp, rho_r, rho_g, drho_r, drho_g, tau_r, tau_g, rho_ao_im_kp)
      NULLIFY (rho_r_sccs)
      NULLIFY (sab_orb)
      my_rebuild_ao = .TRUE.
      my_rebuild_grids = .TRUE.
      my_admm = .FALSE.
      IF (PRESENT(rebuild_ao)) my_rebuild_ao = rebuild_ao
      IF (PRESENT(rebuild_grids)) my_rebuild_grids = rebuild_grids
      IF (PRESENT(admm)) my_admm = admm

      CALL get_qs_env(qs_env, &
                      kpoints=kpoints, &
                      do_kpoints=do_kpoints, &
                      pw_env=pw_env, &
                      dft_control=dft_control)
      IF (PRESENT(pw_env_external)) &
         pw_env => pw_env_external

      nimg = dft_control%nimages

      IF (my_admm) THEN
         CALL get_admm_env(qs_env%admm_env, sab_aux_fit=sab_orb, matrix_s_aux_fit_kp=matrix_s_kp)
      ELSE
         CALL get_qs_env(qs_env, matrix_s_kp=matrix_s_kp)

         IF (do_kpoints) THEN
            CALL get_kpoint_info(kpoints, sab_nl=sab_orb)
         ELSE
            CALL get_qs_env(qs_env, sab_orb=sab_orb)
         END IF
      END IF
      refmatrix => matrix_s_kp(1, 1)%matrix

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      nspins = dft_control%nspins

      CALL qs_rho_get(rho, &
                      tot_rho_r=tot_rho_r, &
                      rho_ao_kp=rho_ao_kp, &
                      rho_ao_im_kp=rho_ao_im_kp, &
                      rho_r=rho_r, &
                      rho_g=rho_g, &
                      drho_r=drho_r, &
                      drho_g=drho_g, &
                      tau_r=tau_r, &
                      tau_g=tau_g, &
                      rho_r_sccs=rho_r_sccs, &
                      complex_rho_ao=rho_ao_is_complex)

      IF (.NOT. ASSOCIATED(tot_rho_r)) THEN
         ALLOCATE (tot_rho_r(nspins))
         tot_rho_r = 0.0_dp
         CALL qs_rho_set(rho, tot_rho_r=tot_rho_r)
      END IF

      ! rho_ao
      IF (my_rebuild_ao .OR. (.NOT. ASSOCIATED(rho_ao_kp))) THEN
         IF (ASSOCIATED(rho_ao_kp)) &
            CALL dbcsr_deallocate_matrix_set(rho_ao_kp)
         ! Create a new density matrix set
         CALL dbcsr_allocate_matrix_set(rho_ao_kp, nspins, nimg)
         CALL qs_rho_set(rho, rho_ao_kp=rho_ao_kp)
         DO i = 1, nspins
            DO ic = 1, nimg
               IF (nspins > 1) THEN
                  IF (i == 1) THEN
                     headline = "DENSITY MATRIX FOR ALPHA SPIN"
                  ELSE
                     headline = "DENSITY MATRIX FOR BETA SPIN"
                  END IF
               ELSE
                  headline = "DENSITY MATRIX"
               END IF
               ALLOCATE (rho_ao_kp(i, ic)%matrix)
               tmatrix => rho_ao_kp(i, ic)%matrix
               CALL dbcsr_create(matrix=tmatrix, template=refmatrix, name=TRIM(headline), &
                                 matrix_type=dbcsr_type_symmetric, nze=0)
               CALL cp_dbcsr_alloc_block_from_nbl(tmatrix, sab_orb)
               CALL dbcsr_set(tmatrix, 0.0_dp)
            END DO
         END DO
         IF (rho_ao_is_complex) THEN
            IF (ASSOCIATED(rho_ao_im_kp)) THEN
               CALL dbcsr_deallocate_matrix_set(rho_ao_im_kp)
            END IF
            CALL dbcsr_allocate_matrix_set(rho_ao_im_kp, nspins, nimg)
            CALL qs_rho_set(rho, rho_ao_im_kp=rho_ao_im_kp)
            DO i = 1, nspins
               DO ic = 1, nimg
                  IF (nspins > 1) THEN
                     IF (i == 1) THEN
                        headline = "IMAGINARY PART OF DENSITY MATRIX FOR ALPHA SPIN"
                     ELSE
                        headline = "IMAGINARY PART OF DENSITY MATRIX FOR BETA SPIN"
                     END IF
                  ELSE
                     headline = "IMAGINARY PART OF DENSITY MATRIX"
                  END IF
                  ALLOCATE (rho_ao_im_kp(i, ic)%matrix)
                  tmatrix => rho_ao_im_kp(i, ic)%matrix
                  CALL dbcsr_create(matrix=tmatrix, template=refmatrix, name=TRIM(headline), &
                                    matrix_type=dbcsr_type_antisymmetric, nze=0)
                  CALL cp_dbcsr_alloc_block_from_nbl(tmatrix, sab_orb)
                  CALL dbcsr_set(tmatrix, 0.0_dp)
               END DO
            END DO
         END IF
      END IF

      ! rho_r
      IF (my_rebuild_grids .OR. .NOT. ASSOCIATED(rho_r)) THEN
         IF (ASSOCIATED(rho_r)) THEN
            DO i = 1, SIZE(rho_r)
               CALL rho_r(i)%release()
            END DO
            DEALLOCATE (rho_r)
         END IF
         ALLOCATE (rho_r(nspins))
         CALL qs_rho_set(rho, rho_r=rho_r)
         DO i = 1, nspins
            CALL auxbas_pw_pool%create_pw(rho_r(i), &
                                          use_data=REALDATA3D, in_space=REALSPACE)
         END DO
      END IF

      ! rho_g
      IF (my_rebuild_grids .OR. .NOT. ASSOCIATED(rho_g)) THEN
         IF (ASSOCIATED(rho_g)) THEN
            DO i = 1, SIZE(rho_g)
               CALL rho_g(i)%release()
            END DO
            DEALLOCATE (rho_g)
         END IF
         ALLOCATE (rho_g(nspins))
         CALL qs_rho_set(rho, rho_g=rho_g)
         DO i = 1, nspins
            CALL auxbas_pw_pool%create_pw(rho_g(i), &
                                          use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         END DO
      END IF

      ! SCCS
      IF (dft_control%do_sccs) THEN
         IF (my_rebuild_grids .OR. (.NOT. ASSOCIATED(rho_r_sccs))) THEN
            IF (ASSOCIATED(rho_r_sccs)) THEN
               CALL rho_r_sccs%release()
               DEALLOCATE (rho_r_sccs)
            END IF
            ALLOCATE (rho_r_sccs)
            CALL qs_rho_set(rho, rho_r_sccs=rho_r_sccs)
            CALL auxbas_pw_pool%create_pw(rho_r_sccs, &
                                          use_data=REALDATA3D, &
                                          in_space=REALSPACE)
            CALL pw_zero(rho_r_sccs)
         END IF
      END IF

      ! allocate drho_r and drho_g if xc_deriv_collocate
      IF (dft_control%drho_by_collocation) THEN
         ! drho_r
         IF (my_rebuild_grids .OR. .NOT. ASSOCIATED(drho_r)) THEN
            IF (ASSOCIATED(drho_r)) THEN
               DO j = 1, SIZE(drho_r, 2)
                  DO i = 1, SIZE(drho_r, 1)
                     CALL drho_r(i, j)%release()
                  END DO
               END DO
               DEALLOCATE (drho_r)
            END IF
            ALLOCATE (drho_r(3, nspins))
            CALL qs_rho_set(rho, drho_r=drho_r)
            DO j = 1, nspins
               DO i = 1, 3
                  CALL auxbas_pw_pool%create_pw(drho_r(i, j), &
                                                use_data=REALDATA3D, in_space=REALSPACE)
               END DO
            END DO
         END IF
         ! drho_g
         IF (my_rebuild_grids .OR. .NOT. ASSOCIATED(drho_g)) THEN
            IF (ASSOCIATED(drho_g)) THEN
               DO j = 1, SIZE(drho_g, 2)
                  DO i = 1, SIZE(drho_r, 1)
                     CALL drho_g(i, j)%release()
                  END DO
               END DO
               DEALLOCATE (drho_g)
            END IF
            ALLOCATE (drho_g(3, nspins))
            CALL qs_rho_set(rho, drho_g=drho_g)
            DO j = 1, nspins
               DO i = 1, 3
                  CALL auxbas_pw_pool%create_pw(drho_g(i, j), &
                                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
               END DO
            END DO
         END IF
      END IF

      ! allocate tau_r and tau_g if use_kinetic_energy_density
      IF (dft_control%use_kinetic_energy_density) THEN
         ! tau_r
         IF (my_rebuild_grids .OR. .NOT. ASSOCIATED(tau_r)) THEN
            IF (ASSOCIATED(tau_r)) THEN
               DO i = 1, SIZE(tau_r)
                  CALL tau_r(i)%release()
               END DO
               DEALLOCATE (tau_r)
            END IF
            ALLOCATE (tau_r(nspins))
            CALL qs_rho_set(rho, tau_r=tau_r)
            DO i = 1, nspins
               CALL auxbas_pw_pool%create_pw(tau_r(i), &
                                             use_data=REALDATA3D, in_space=REALSPACE)
            END DO
         END IF

         ! tau_g
         IF (my_rebuild_grids .OR. .NOT. ASSOCIATED(tau_g)) THEN
            IF (ASSOCIATED(tau_g)) THEN
               DO i = 1, SIZE(tau_g)
                  CALL tau_g(i)%release()
               END DO
               DEALLOCATE (tau_g)
            END IF
            ALLOCATE (tau_g(nspins))
            CALL qs_rho_set(rho, tau_g=tau_g)
            DO i = 1, nspins
               CALL auxbas_pw_pool%create_pw(tau_g(i), &
                                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            END DO
         END IF
      END IF ! use_kinetic_energy_density

      CALL timestop(handle)

   END SUBROUTINE qs_rho_rebuild

! **************************************************************************************************
!> \brief updates rho_r and rho_g to the rho%rho_ao.
!>      if use_kinetic_energy_density also computes tau_r and tau_g
!>      this works for all ground state and ground state response methods
!> \param rho_struct the rho structure that should be updated
!> \param qs_env the qs_env rho_struct refers to
!>        the integrated charge in r space
!> \param rho_xc_external ...
!> \param local_rho_set ...
!> \param task_list_external external task list
!> \param task_list_external_soft external task list (soft_version)
!> \param pw_env_external    external plane wave environment
!> \param para_env_external  external MPI environment
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE qs_rho_update_rho(rho_struct, qs_env, rho_xc_external, local_rho_set, &
                                task_list_external, task_list_external_soft, &
                                pw_env_external, para_env_external)
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_struct
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), OPTIONAL, POINTER               :: rho_xc_external
      TYPE(local_rho_type), OPTIONAL, POINTER            :: local_rho_set
      TYPE(task_list_type), OPTIONAL, POINTER            :: task_list_external, &
                                                            task_list_external_soft
      TYPE(pw_env_type), OPTIONAL, POINTER               :: pw_env_external
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_external

      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL get_qs_env(qs_env, dft_control=dft_control, &
                      atomic_kind_set=atomic_kind_set, &
                      para_env=para_env)
      IF (PRESENT(para_env_external)) para_env => para_env_external

      IF (dft_control%qs_control%semi_empirical .OR. &
          dft_control%qs_control%dftb .OR. dft_control%qs_control%xtb) THEN

         CALL qs_rho_set(rho_struct, rho_r_valid=.FALSE., rho_g_valid=.FALSE.)

      ELSEIF (dft_control%qs_control%lrigpw) THEN
         CPASSERT(.NOT. dft_control%use_kinetic_energy_density)
         CPASSERT(.NOT. dft_control%drho_by_collocation)
         CALL qs_rho_get(rho_struct, rho_ao_kp=rho_ao_kp)
         CALL get_qs_env(qs_env, ks_env=ks_env)
         CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
         CALL get_qs_env(qs_env, lri_env=lri_env, lri_density=lri_density)
         CALL calculate_lri_densities(lri_env, lri_density, qs_env, rho_ao_kp, cell_to_index, &
                                      lri_rho_struct=rho_struct, &
                                      atomic_kind_set=atomic_kind_set, &
                                      para_env=para_env, &
                                      response_density=.FALSE.)
         CALL set_qs_env(qs_env, lri_density=lri_density)
         CALL qs_rho_set(rho_struct, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)

      ELSEIF (dft_control%qs_control%rigpw) THEN
         CPASSERT(.NOT. dft_control%use_kinetic_energy_density)
         CPASSERT(.NOT. dft_control%drho_by_collocation)
         CALL get_qs_env(qs_env, lri_env=lri_env)
         CALL qs_rho_get(rho_struct, rho_ao=rho_ao)
         CALL calculate_ri_densities(lri_env, qs_env, rho_ao, &
                                     lri_rho_struct=rho_struct, &
                                     atomic_kind_set=atomic_kind_set, &
                                     para_env=para_env)
         CALL qs_rho_set(rho_struct, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)

      ELSE
         CALL qs_rho_update_rho_low(rho_struct=rho_struct, qs_env=qs_env, &
                                    rho_xc_external=rho_xc_external, &
                                    local_rho_set=local_rho_set, &
                                    task_list_external=task_list_external, &
                                    task_list_external_soft=task_list_external_soft, &
                                    pw_env_external=pw_env_external, &
                                    para_env_external=para_env_external)

      END IF

   END SUBROUTINE qs_rho_update_rho

! **************************************************************************************************
!> \brief updates rho_r and rho_g to the rho%rho_ao.
!>      if use_kinetic_energy_density also computes tau_r and tau_g
!> \param rho_struct the rho structure that should be updated
!> \param qs_env the qs_env rho_struct refers to
!>        the integrated charge in r space
!> \param rho_xc_external rho structure for GAPW_XC
!> \param local_rho_set ...
!> \param pw_env_external    external plane wave environment
!> \param task_list_external external task list (use for default and GAPW)
!> \param task_list_external_soft external task list (soft density for GAPW_XC)
!> \param para_env_external ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE qs_rho_update_rho_low(rho_struct, qs_env, rho_xc_external, &
                                    local_rho_set, pw_env_external, &
                                    task_list_external, task_list_external_soft, &
                                    para_env_external)
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_struct
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), OPTIONAL, POINTER               :: rho_xc_external
      TYPE(local_rho_type), OPTIONAL, POINTER            :: local_rho_set
      TYPE(pw_env_type), OPTIONAL, POINTER               :: pw_env_external
      TYPE(task_list_type), OPTIONAL, POINTER            :: task_list_external, &
                                                            task_list_external_soft
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_external

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_update_rho_low'

      INTEGER                                            :: handle, img, ispin, nimg, nspins
      LOGICAL                                            :: gapw, gapw_xc
      REAL(KIND=dp)                                      :: dum
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r, tot_rho_r_xc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp, rho_xc_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab
      TYPE(oce_matrix_type), POINTER                     :: oce
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_type), DIMENSION(:), POINTER               :: rho_g, rho_r, rho_xc_g, rho_xc_r, tau_g, &
                                                            tau_r, tau_xc_g, tau_xc_r
      TYPE(pw_type), DIMENSION(:, :), POINTER            :: drho_g, drho_r, drho_xc_g, drho_xc_r
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_xc
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      TYPE(task_list_type), POINTER                      :: task_list

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, rho_xc, ks_env, rho_ao, rho_r, rho_g, drho_r, drho_g, tau_r, tau_g)
      NULLIFY (rho_xc_ao, rho_xc_g, rho_xc_r, drho_xc_g, tau_xc_r, tau_xc_g, tot_rho_r, tot_rho_r_xc)
      NULLIFY (para_env, pw_env, atomic_kind_set)

      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      dft_control=dft_control, &
                      atomic_kind_set=atomic_kind_set)

      CALL qs_rho_get(rho_struct, &
                      rho_r=rho_r, &
                      rho_g=rho_g, &
                      tot_rho_r=tot_rho_r, &
                      drho_r=drho_r, &
                      drho_g=drho_g, &
                      tau_r=tau_r, &
                      tau_g=tau_g)

      CALL get_qs_env(qs_env, task_list=task_list, &
                      para_env=para_env, pw_env=pw_env)
      IF (PRESENT(pw_env_external)) pw_env => pw_env_external
      IF (PRESENT(task_list_external)) task_list => task_list_external
      IF (PRESENT(para_env_external)) para_env => para_env_external

      nspins = dft_control%nspins
      nimg = dft_control%nimages
      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc

      CALL qs_rho_get(rho_struct, rho_ao_kp=rho_ao_kp)
      DO ispin = 1, nspins
         rho_ao => rho_ao_kp(ispin, :)
         CALL calculate_rho_elec(matrix_p_kp=rho_ao, &
                                 rho=rho_r(ispin), &
                                 rho_gspace=rho_g(ispin), &
                                 total_rho=tot_rho_r(ispin), &
                                 ks_env=ks_env, soft_valid=gapw, &
                                 task_list_external=task_list_external, &
                                 pw_env_external=pw_env_external)
      END DO
      CALL qs_rho_set(rho_struct, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)

      IF (gapw_xc) THEN
         IF (PRESENT(rho_xc_external)) THEN
            rho_xc => rho_xc_external
         ELSE
            CALL get_qs_env(qs_env=qs_env, rho_xc=rho_xc)
         END IF
         CALL qs_rho_get(rho_xc, &
                         rho_ao_kp=rho_xc_ao, &
                         rho_r=rho_xc_r, &
                         rho_g=rho_xc_g, &
                         tot_rho_r=tot_rho_r_xc)
         ! copy rho_ao into rho_xc_ao
         DO ispin = 1, nspins
            DO img = 1, nimg
               CALL dbcsr_copy(rho_xc_ao(ispin, img)%matrix, rho_ao_kp(ispin, img)%matrix)
            END DO
         END DO
         DO ispin = 1, nspins
            rho_ao => rho_xc_ao(ispin, :)
            CALL calculate_rho_elec(matrix_p_kp=rho_ao, &
                                    rho=rho_xc_r(ispin), &
                                    rho_gspace=rho_xc_g(ispin), &
                                    total_rho=tot_rho_r_xc(ispin), &
                                    ks_env=ks_env, soft_valid=gapw_xc, &
                                    task_list_external=task_list_external_soft, &
                                    pw_env_external=pw_env_external)
         END DO
         CALL qs_rho_set(rho_xc, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)
      END IF

      ! GAPW o GAPW_XC require the calculation of hard and soft local densities
      IF (gapw .OR. gapw_xc) THEN
         CALL get_qs_env(qs_env=qs_env, &
                         rho_atom_set=rho_atom_set, &
                         qs_kind_set=qs_kind_set, &
                         oce=oce, sab_orb=sab)
         IF (PRESENT(local_rho_set)) rho_atom_set => local_rho_set%rho_atom_set
         CPASSERT(ASSOCIATED(rho_atom_set))
         CALL qs_rho_get(rho_struct, rho_ao_kp=rho_ao_kp)
         CALL calculate_rho_atom_coeff(qs_env, rho_ao_kp, rho_atom_set, qs_kind_set, oce, sab, para_env)
      END IF

      IF (.NOT. gapw_xc) THEN
         ! if needed compute also the gradient of the density
         IF (dft_control%drho_by_collocation) THEN
            CALL qs_rho_get(rho_struct, rho_ao_kp=rho_ao_kp)
            CPASSERT(.NOT. PRESENT(task_list_external))
            DO ispin = 1, nspins
               rho_ao => rho_ao_kp(ispin, :)
               CALL calculate_drho_elec(matrix_p_kp=rho_ao, &
                                        drho=drho_r(:, ispin), &
                                        drho_gspace=drho_g(:, ispin), &
                                        qs_env=qs_env, soft_valid=gapw)
            END DO
            CALL qs_rho_set(rho_struct, drho_r_valid=.TRUE., drho_g_valid=.TRUE.)
         END IF
         ! if needed compute also the kinetic energy density
         IF (dft_control%use_kinetic_energy_density) THEN
            CALL qs_rho_get(rho_struct, rho_ao_kp=rho_ao_kp)
            DO ispin = 1, nspins
               rho_ao => rho_ao_kp(ispin, :)
               CALL calculate_rho_elec(matrix_p_kp=rho_ao, &
                                       rho=tau_r(ispin), &
                                       rho_gspace=tau_g(ispin), &
                                       total_rho=dum, & ! presumably not meaningful
                                       ks_env=ks_env, soft_valid=gapw, &
                                       compute_tau=.TRUE., &
                                       task_list_external=task_list_external, &
                                       pw_env_external=pw_env_external)
            END DO
            CALL qs_rho_set(rho_struct, tau_r_valid=.TRUE., tau_g_valid=.TRUE.)
         END IF
      ELSE
         CALL qs_rho_get(rho_xc, &
                         drho_r=drho_xc_r, &
                         drho_g=drho_xc_g, &
                         tau_r=tau_xc_r, &
                         tau_g=tau_xc_g)
         ! if needed compute also the gradient of the density
         IF (dft_control%drho_by_collocation) THEN
            CPASSERT(.NOT. PRESENT(task_list_external))
            DO ispin = 1, nspins
               rho_ao => rho_xc_ao(ispin, :)
               CALL calculate_drho_elec(matrix_p_kp=rho_ao, &
                                        drho=drho_xc_r(:, ispin), &
                                        drho_gspace=drho_xc_g(:, ispin), &
                                        qs_env=qs_env, soft_valid=gapw_xc)
            END DO
            CALL qs_rho_set(rho_xc, drho_r_valid=.TRUE., drho_g_valid=.TRUE.)
         END IF
         ! if needed compute also the kinetic energy density
         IF (dft_control%use_kinetic_energy_density) THEN
            DO ispin = 1, nspins
               rho_ao => rho_xc_ao(ispin, :)
               CALL calculate_rho_elec(matrix_p_kp=rho_ao, &
                                       rho=tau_xc_r(ispin), &
                                       rho_gspace=tau_xc_g(ispin), &
                                       ks_env=ks_env, soft_valid=gapw_xc, &
                                       compute_tau=.TRUE., &
                                       task_list_external=task_list_external_soft, &
                                       pw_env_external=pw_env_external)
            END DO
            CALL qs_rho_set(rho_xc, tau_r_valid=.TRUE., tau_g_valid=.TRUE.)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_rho_update_rho_low

! **************************************************************************************************
!> \brief updates rho_r and rho_g to the rho%rho_ao.
!>      if use_kinetic_energy_density also computes tau_r and tau_g
!> \param rho_struct the rho structure that should be updated
!> \param qs_env the qs_env rho_struct refers to
!>        the integrated charge in r space
!> \param pw_env_external    external plane wave environment
!> \param task_list_external external task list
!> \param para_env_external ...
!> \param tddfpt_lri_env ...
!> \param tddfpt_lri_density ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE qs_rho_update_tddfpt(rho_struct, qs_env, pw_env_external, task_list_external, &
                                   para_env_external, tddfpt_lri_env, tddfpt_lri_density)
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_struct
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_env_type), OPTIONAL, POINTER               :: pw_env_external
      TYPE(task_list_type), OPTIONAL, POINTER            :: task_list_external
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_external
      TYPE(lri_environment_type), OPTIONAL, POINTER      :: tddfpt_lri_env
      TYPE(lri_density_type), OPTIONAL, POINTER          :: tddfpt_lri_density

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_update_tddfpt'

      INTEGER                                            :: handle, ispin, nspins
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: lri_response
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_type), DIMENSION(:), POINTER               :: rho_g, rho_r
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(task_list_type), POINTER                      :: task_list

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      dft_control=dft_control, &
                      atomic_kind_set=atomic_kind_set, &
                      task_list=task_list, &
                      para_env=para_env, &
                      pw_env=pw_env)
      IF (PRESENT(pw_env_external)) pw_env => pw_env_external
      IF (PRESENT(task_list_external)) task_list => task_list_external
      IF (PRESENT(para_env_external)) para_env => para_env_external

      CALL qs_rho_get(rho_struct, &
                      rho_r=rho_r, &
                      rho_g=rho_g, &
                      tot_rho_r=tot_rho_r)

      nspins = dft_control%nspins

      lri_response = PRESENT(tddfpt_lri_env)
      IF (lri_response) THEN
         CPASSERT(PRESENT(tddfpt_lri_density))
      END IF

      CPASSERT(.NOT. dft_control%drho_by_collocation)
      CPASSERT(.NOT. dft_control%use_kinetic_energy_density)
      CPASSERT(.NOT. dft_control%qs_control%gapw)
      CPASSERT(.NOT. dft_control%qs_control%gapw_xc)

      IF (lri_response) THEN
         CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
         CALL qs_rho_get(rho_struct, rho_ao_kp=rho_ao_kp)
         CALL calculate_lri_densities(tddfpt_lri_env, tddfpt_lri_density, qs_env, rho_ao_kp, cell_to_index, &
                                      lri_rho_struct=rho_struct, &
                                      atomic_kind_set=atomic_kind_set, &
                                      para_env=para_env, &
                                      response_density=lri_response)
         CALL qs_rho_set(rho_struct, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)
      ELSE
         CALL qs_rho_get(rho_struct, rho_ao_kp=rho_ao_kp)
         DO ispin = 1, nspins
            rho_ao => rho_ao_kp(ispin, :)
            CALL calculate_rho_elec(matrix_p_kp=rho_ao, &
                                    rho=rho_r(ispin), &
                                    rho_gspace=rho_g(ispin), &
                                    total_rho=tot_rho_r(ispin), &
                                    ks_env=ks_env, &
                                    task_list_external=task_list_external, &
                                    pw_env_external=pw_env_external)
         END DO
         CALL qs_rho_set(rho_struct, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_rho_update_tddfpt

! **************************************************************************************************
!> \brief Allocate a density structure and fill it with data from an input structure
!>        SIZE(rho_input) == mspin == 1  direct copy
!>        SIZE(rho_input) == mspin == 2  direct copy of alpha and beta spin
!>        SIZE(rho_input) == 1 AND mspin == 2  copy rho/2 into alpha and beta spin
!> \param rho_input ...
!> \param rho_output ...
!> \param auxbas_pw_pool ...
!> \param mspin ...
! **************************************************************************************************
   SUBROUTINE qs_rho_copy(rho_input, rho_output, auxbas_pw_pool, mspin)

      TYPE(qs_rho_type), INTENT(IN)                      :: rho_input
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_output
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      INTEGER, INTENT(IN)                                :: mspin

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_rho_copy'

      INTEGER                                            :: handle, i, j, nspins
      LOGICAL :: complex_rho_ao, drho_g_valid_in, drho_r_valid_in, rho_g_valid_in, rho_r_valid_in, &
         soft_valid_in, tau_g_valid_in, tau_r_valid_in
      REAL(KIND=dp)                                      :: ospin
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_g_in, tot_rho_g_out, &
                                                            tot_rho_r_in, tot_rho_r_out
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao_im_in, rho_ao_im_out, rho_ao_in, &
                                                            rho_ao_out
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp_in
      TYPE(pw_type), DIMENSION(:), POINTER               :: rho_g_in, rho_g_out, rho_r_in, &
                                                            rho_r_out, tau_g_in, tau_g_out, &
                                                            tau_r_in, tau_r_out
      TYPE(pw_type), DIMENSION(:, :), POINTER            :: drho_g_in, drho_g_out, drho_r_in, &
                                                            drho_r_out
      TYPE(pw_type), POINTER                             :: rho_r_sccs_in, rho_r_sccs_out

      CALL timeset(routineN, handle)

      CPASSERT(mspin == 1 .OR. mspin == 2)
      ospin = 1._dp/REAL(mspin, KIND=dp)

      CALL qs_rho_clear(rho_output)

      NULLIFY (rho_ao_in, rho_ao_kp_in, rho_ao_im_in, rho_r_in, rho_g_in, drho_r_in, &
               drho_g_in, tau_r_in, tau_g_in, tot_rho_r_in, tot_rho_g_in, rho_r_sccs_in)

      CALL qs_rho_get(rho_input, &
                      rho_ao=rho_ao_in, &
                      rho_ao_kp=rho_ao_kp_in, &
                      rho_ao_im=rho_ao_im_in, &
                      rho_r=rho_r_in, &
                      rho_g=rho_g_in, &
                      drho_r=drho_r_in, &
                      drho_g=drho_g_in, &
                      tau_r=tau_r_in, &
                      tau_g=tau_g_in, &
                      tot_rho_r=tot_rho_r_in, &
                      tot_rho_g=tot_rho_g_in, &
                      rho_g_valid=rho_g_valid_in, &
                      rho_r_valid=rho_r_valid_in, &
                      drho_g_valid=drho_g_valid_in, &
                      drho_r_valid=drho_r_valid_in, &
                      tau_r_valid=tau_r_valid_in, &
                      tau_g_valid=tau_g_valid_in, &
                      rho_r_sccs=rho_r_sccs_in, &
                      soft_valid=soft_valid_in, &
                      complex_rho_ao=complex_rho_ao)

      NULLIFY (rho_ao_out, rho_ao_im_out, rho_r_out, rho_g_out, drho_r_out, &
               drho_g_out, tau_r_out, tau_g_out, tot_rho_r_out, tot_rho_g_out, rho_r_sccs_out)
      ! rho_ao
      IF (ASSOCIATED(rho_ao_in)) THEN
         nspins = SIZE(rho_ao_in)
         CPASSERT(mspin >= nspins)
         CALL dbcsr_allocate_matrix_set(rho_ao_out, mspin)
         CALL qs_rho_set(rho_output, rho_ao=rho_ao_out)
         IF (mspin > nspins) THEN
            DO i = 1, mspin
               ALLOCATE (rho_ao_out(i)%matrix)
               CALL dbcsr_copy(rho_ao_out(i)%matrix, rho_ao_in(1)%matrix, name="RHO copy")
               CALL dbcsr_scale(rho_ao_out(i)%matrix, ospin)
            END DO
         ELSE
            DO i = 1, nspins
               ALLOCATE (rho_ao_out(i)%matrix)
               CALL dbcsr_copy(rho_ao_out(i)%matrix, rho_ao_in(i)%matrix, name="RHO copy")
            END DO
         END IF
      END IF

      ! rho_ao_kp
      ! only for non-kp, we could probably just copy this pointer, should work also for non-kp?
      !IF (ASSOCIATED(rho_ao_kp_in)) THEN
      !   CPABORT("Copy not available")
      !END IF

      ! rho_ao_im
      IF (ASSOCIATED(rho_ao_im_in)) THEN
         nspins = SIZE(rho_ao_im_in)
         CPASSERT(mspin >= nspins)
         CALL dbcsr_allocate_matrix_set(rho_ao_im_out, mspin)
         CALL qs_rho_set(rho_output, rho_ao_im=rho_ao_im_out)
         IF (mspin > nspins) THEN
            DO i = 1, mspin
               ALLOCATE (rho_ao_im_out(i)%matrix)
               CALL dbcsr_copy(rho_ao_im_out(i)%matrix, rho_ao_im_in(1)%matrix, name="RHO copy")
               CALL dbcsr_scale(rho_ao_im_out(i)%matrix, ospin)
            END DO
         ELSE
            DO i = 1, nspins
               ALLOCATE (rho_ao_im_out(i)%matrix)
               CALL dbcsr_copy(rho_ao_im_out(i)%matrix, rho_ao_im_in(i)%matrix, name="RHO copy")
            END DO
         END IF
      END IF

      ! rho_r
      IF (ASSOCIATED(rho_r_in)) THEN
         nspins = SIZE(rho_r_in)
         CPASSERT(mspin >= nspins)
         ALLOCATE (rho_r_out(mspin))
         CALL qs_rho_set(rho_output, rho_r=rho_r_out)
         IF (mspin > nspins) THEN
            DO i = 1, mspin
               CALL auxbas_pw_pool%create_pw(rho_r_out(i), &
                                             use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_copy(rho_r_in(1), rho_r_out(i))
               CALL pw_scale(rho_r_out(i), ospin)
            END DO
         ELSE
            DO i = 1, nspins
               CALL auxbas_pw_pool%create_pw(rho_r_out(i), &
                                             use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_copy(rho_r_in(i), rho_r_out(i))
            END DO
         END IF
      END IF

      ! rho_g
      IF (ASSOCIATED(rho_g_in)) THEN
         nspins = SIZE(rho_g_in)
         CPASSERT(mspin >= nspins)
         ALLOCATE (rho_g_out(mspin))
         CALL qs_rho_set(rho_output, rho_g=rho_g_out)
         IF (mspin > nspins) THEN
            DO i = 1, mspin
               CALL auxbas_pw_pool%create_pw(rho_g_out(i), &
                                             use_data=COMPLEXDATA1D, &
                                             in_space=RECIPROCALSPACE)
               CALL pw_copy(rho_g_in(1), rho_g_out(i))
               CALL pw_scale(rho_g_out(i), ospin)
            END DO
         ELSE
            DO i = 1, nspins
               CALL auxbas_pw_pool%create_pw(rho_g_out(i), &
                                             use_data=COMPLEXDATA1D, &
                                             in_space=RECIPROCALSPACE)
               CALL pw_copy(rho_g_in(i), rho_g_out(i))
            END DO
         END IF
      END IF

      ! SCCS
      IF (ASSOCIATED(rho_r_sccs_in)) THEN
         CALL qs_rho_set(rho_output, rho_r_sccs=rho_r_sccs_out)
         CALL auxbas_pw_pool%create_pw(rho_r_sccs_out, &
                                       in_space=REALSPACE, &
                                       use_data=REALDATA3D)
         CALL pw_copy(rho_r_sccs_in, rho_r_sccs_out)
      END IF

      ! drho_r
      IF (ASSOCIATED(drho_r_in)) THEN
         nspins = SIZE(drho_r_in)
         CPASSERT(mspin >= nspins)
         ALLOCATE (drho_r_out(3, mspin))
         CALL qs_rho_set(rho_output, drho_r=drho_r_out)
         IF (mspin > nspins) THEN
            DO j = 1, mspin
               DO i = 1, 3
                  CALL auxbas_pw_pool%create_pw(drho_r_out(i, j), &
                                                use_data=REALDATA3D, in_space=REALSPACE)
                  CALL pw_copy(drho_r_in(i, 1), drho_r_out(i, j))
                  CALL pw_scale(drho_r_out(i, j), ospin)
               END DO
            END DO
         ELSE
            DO j = 1, nspins
               DO i = 1, 3
                  CALL auxbas_pw_pool%create_pw(drho_r_out(i, j), &
                                                use_data=REALDATA3D, in_space=REALSPACE)
                  CALL pw_copy(drho_r_in(i, j), drho_r_out(i, j))
               END DO
            END DO
         END IF
      END IF

      ! drho_g
      IF (ASSOCIATED(drho_g_in)) THEN
         nspins = SIZE(drho_g_in)
         CPASSERT(mspin >= nspins)
         ALLOCATE (drho_g_out(3, mspin))
         CALL qs_rho_set(rho_output, drho_g=drho_g_out)
         IF (mspin > nspins) THEN
            DO j = 1, mspin
               DO i = 1, 3
                  CALL auxbas_pw_pool%create_pw(drho_g_out(i, j), &
                                                use_data=COMPLEXDATA1D, &
                                                in_space=RECIPROCALSPACE)
                  CALL pw_copy(drho_g_in(i, 1), drho_g_out(i, j))
                  CALL pw_scale(drho_g_out(i, j), ospin)
               END DO
            END DO
         ELSE
            DO j = 1, nspins
               DO i = 1, 3
                  CALL auxbas_pw_pool%create_pw(drho_g_out(i, j), &
                                                use_data=COMPLEXDATA1D, &
                                                in_space=RECIPROCALSPACE)
                  CALL pw_copy(drho_g_in(i, j), drho_g_out(i, j))
               END DO
            END DO
         END IF
      END IF

      ! tau_r
      IF (ASSOCIATED(tau_r_in)) THEN
         nspins = SIZE(tau_r_in)
         CPASSERT(mspin >= nspins)
         ALLOCATE (tau_r_out(mspin))
         CALL qs_rho_set(rho_output, tau_r=tau_r_out)
         IF (mspin > nspins) THEN
            DO i = 1, mspin
               CALL auxbas_pw_pool%create_pw(tau_r_out(i), &
                                             use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_copy(tau_r_in(1), tau_r_out(i))
               CALL pw_scale(tau_r_out(i), ospin)
            END DO
         ELSE
            DO i = 1, nspins
               CALL auxbas_pw_pool%create_pw(tau_r_out(i), &
                                             use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_copy(tau_r_in(i), tau_r_out(i))
            END DO
         END IF
      END IF

      ! tau_g
      IF (ASSOCIATED(tau_g_in)) THEN
         nspins = SIZE(tau_g_in)
         CPASSERT(mspin >= nspins)
         ALLOCATE (tau_g_out(mspin))
         CALL qs_rho_set(rho_output, tau_g=tau_g_out)
         IF (mspin > nspins) THEN
            DO i = 1, mspin
               CALL auxbas_pw_pool%create_pw(tau_g_out(i), &
                                             use_data=COMPLEXDATA1D, &
                                             in_space=RECIPROCALSPACE)
               CALL pw_copy(tau_g_in(1), tau_g_out(i))
               CALL pw_scale(tau_g_out(i), ospin)
            END DO
         ELSE
            DO i = 1, nspins
               CALL auxbas_pw_pool%create_pw(tau_g_out(i), &
                                             use_data=COMPLEXDATA1D, &
                                             in_space=RECIPROCALSPACE)
               CALL pw_copy(tau_g_in(i), tau_g_out(i))
            END DO
         END IF
      END IF

      ! tot_rho_r
      IF (ASSOCIATED(tot_rho_r_in)) THEN
         nspins = SIZE(tot_rho_r_in)
         CPASSERT(mspin >= nspins)
         ALLOCATE (tot_rho_r_out(mspin))
         CALL qs_rho_set(rho_output, tot_rho_r=tot_rho_r_out)
         IF (mspin > nspins) THEN
            DO i = 1, mspin
               tot_rho_r_out(i) = tot_rho_r_in(1)*ospin
            END DO
         ELSE
            DO i = 1, nspins
               tot_rho_r_out(i) = tot_rho_r_in(i)
            END DO
         END IF
      END IF

      ! tot_rho_g
      IF (ASSOCIATED(tot_rho_g_in)) THEN
         nspins = SIZE(tot_rho_g_in)
         CPASSERT(mspin >= nspins)
         ALLOCATE (tot_rho_g_out(mspin))
         CALL qs_rho_set(rho_output, tot_rho_g=tot_rho_g_out)
         IF (mspin > nspins) THEN
            DO i = 1, mspin
               tot_rho_g_out(i) = tot_rho_g_in(1)*ospin
            END DO
         ELSE
            DO i = 1, nspins
               tot_rho_g_out(i) = tot_rho_g_in(i)
            END DO
         END IF
      END IF

      CALL qs_rho_set(rho_output, &
                      rho_g_valid=rho_g_valid_in, &
                      rho_r_valid=rho_r_valid_in, &
                      drho_g_valid=drho_g_valid_in, &
                      drho_r_valid=drho_r_valid_in, &
                      tau_r_valid=tau_r_valid_in, &
                      tau_g_valid=tau_g_valid_in, &
                      soft_valid=soft_valid_in, &
                      complex_rho_ao=complex_rho_ao)

      CALL timestop(handle)

   END SUBROUTINE qs_rho_copy

! **************************************************************************************************
!> \brief rhoa = alpha*rhoa+beta*rhob
!> \param rhoa ...
!> \param rhob ...
!> \param alpha ...
!> \param beta ...
! **************************************************************************************************
   SUBROUTINE qs_rho_scale_and_add(rhoa, rhob, alpha, beta)

      TYPE(qs_rho_type), INTENT(IN)                      :: rhoa, rhob
      REAL(KIND=dp), INTENT(IN)                          :: alpha, beta

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_rho_scale_and_add'

      INTEGER                                            :: handle, i, j, nspina, nspinb, nspins
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_g_a, tot_rho_g_b, tot_rho_r_a, &
                                                            tot_rho_r_b
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao_a, rho_ao_b, rho_ao_im_a, &
                                                            rho_ao_im_b
      TYPE(pw_type), DIMENSION(:), POINTER               :: rho_g_a, rho_g_b, rho_r_a, rho_r_b, &
                                                            tau_g_a, tau_g_b, tau_r_a, tau_r_b
      TYPE(pw_type), DIMENSION(:, :), POINTER            :: drho_g_a, drho_g_b, drho_r_a, drho_r_b
      TYPE(pw_type), POINTER                             :: rho_r_sccs_a, rho_r_sccs_b

      CALL timeset(routineN, handle)

      NULLIFY (rho_ao_a, rho_ao_im_a, rho_r_a, rho_g_a, drho_r_a, &
               drho_g_a, tau_r_a, tau_g_a, tot_rho_r_a, tot_rho_g_a, rho_r_sccs_a)

      CALL qs_rho_get(rhoa, &
                      rho_ao=rho_ao_a, &
                      rho_ao_im=rho_ao_im_a, &
                      rho_r=rho_r_a, &
                      rho_g=rho_g_a, &
                      drho_r=drho_r_a, &
                      drho_g=drho_g_a, &
                      tau_r=tau_r_a, &
                      tau_g=tau_g_a, &
                      tot_rho_r=tot_rho_r_a, &
                      tot_rho_g=tot_rho_g_a, &
                      rho_r_sccs=rho_r_sccs_a)

      NULLIFY (rho_ao_b, rho_ao_im_b, rho_r_b, rho_g_b, drho_r_b, &
               drho_g_b, tau_r_b, tau_g_b, tot_rho_r_b, tot_rho_g_b, rho_r_sccs_b)

      CALL qs_rho_get(rhob, &
                      rho_ao=rho_ao_b, &
                      rho_ao_im=rho_ao_im_b, &
                      rho_r=rho_r_b, &
                      rho_g=rho_g_b, &
                      drho_r=drho_r_b, &
                      drho_g=drho_g_b, &
                      tau_r=tau_r_b, &
                      tau_g=tau_g_b, &
                      tot_rho_r=tot_rho_r_b, &
                      tot_rho_g=tot_rho_g_b, &
                      rho_r_sccs=rho_r_sccs_b)
      ! rho_ao
      IF (ASSOCIATED(rho_ao_a) .AND. ASSOCIATED(rho_ao_b)) THEN
         nspina = SIZE(rho_ao_a)
         nspinb = SIZE(rho_ao_b)
         nspins = MIN(nspina, nspinb)
         DO i = 1, nspins
            CALL dbcsr_add(rho_ao_a(i)%matrix, rho_ao_b(i)%matrix, alpha, beta)
         END DO
      END IF

      ! rho_ao_im
      IF (ASSOCIATED(rho_ao_im_a) .AND. ASSOCIATED(rho_ao_im_b)) THEN
         nspina = SIZE(rho_ao_im_a)
         nspinb = SIZE(rho_ao_im_b)
         nspins = MIN(nspina, nspinb)
         DO i = 1, nspins
            CALL dbcsr_add(rho_ao_im_a(i)%matrix, rho_ao_im_b(i)%matrix, alpha, beta)
         END DO
      END IF

      ! rho_r
      IF (ASSOCIATED(rho_r_a) .AND. ASSOCIATED(rho_r_b)) THEN
         nspina = SIZE(rho_ao_a)
         nspinb = SIZE(rho_ao_b)
         nspins = MIN(nspina, nspinb)
         DO i = 1, nspins
            CALL pw_axpy(rho_r_b(i), rho_r_a(i), beta, alpha)
         END DO
      END IF

      ! rho_g
      IF (ASSOCIATED(rho_g_a) .AND. ASSOCIATED(rho_g_b)) THEN
         nspina = SIZE(rho_ao_a)
         nspinb = SIZE(rho_ao_b)
         nspins = MIN(nspina, nspinb)
         DO i = 1, nspins
            CALL pw_axpy(rho_g_b(i), rho_g_a(i), beta, alpha)
         END DO
      END IF

      ! SCCS
      IF (ASSOCIATED(rho_r_sccs_a) .AND. ASSOCIATED(rho_r_sccs_b)) THEN
         CALL pw_axpy(rho_r_sccs_b, rho_r_sccs_a, beta, alpha)
      END IF

      ! drho_r
      IF (ASSOCIATED(drho_r_a) .AND. ASSOCIATED(drho_r_b)) THEN
         CPASSERT(ALL(SHAPE(drho_r_a) == SHAPE(drho_r_b))) ! not implemented
         DO j = 1, SIZE(drho_r_a, 2)
            DO i = 1, SIZE(drho_r_a, 1)
               CALL pw_axpy(drho_r_b(i, j), drho_r_a(i, j), beta, alpha)
            END DO
         END DO
      END IF

      ! drho_g
      IF (ASSOCIATED(drho_g_a) .AND. ASSOCIATED(drho_g_b)) THEN
         CPASSERT(ALL(SHAPE(drho_g_a) == SHAPE(drho_g_b))) ! not implemented
         DO j = 1, SIZE(drho_g_a, 2)
            DO i = 1, SIZE(drho_g_a, 1)
               CALL pw_axpy(drho_g_b(i, j), drho_g_a(i, j), beta, alpha)
            END DO
         END DO
      END IF

      ! tau_r
      IF (ASSOCIATED(tau_r_a) .AND. ASSOCIATED(tau_r_b)) THEN
         nspina = SIZE(rho_ao_a)
         nspinb = SIZE(rho_ao_b)
         nspins = MIN(nspina, nspinb)
         DO i = 1, nspins
            CALL pw_axpy(tau_r_b(i), tau_r_a(i), beta, alpha)
         END DO
      END IF

      ! tau_g
      IF (ASSOCIATED(tau_g_a) .AND. ASSOCIATED(tau_g_b)) THEN
         nspina = SIZE(rho_ao_a)
         nspinb = SIZE(rho_ao_b)
         nspins = MIN(nspina, nspinb)
         DO i = 1, nspins
            CALL pw_axpy(tau_g_b(i), tau_g_a(i), beta, alpha)
         END DO
      END IF

      ! tot_rho_r
      IF (ASSOCIATED(tot_rho_r_a) .AND. ASSOCIATED(tot_rho_r_b)) THEN
         nspina = SIZE(rho_ao_a)
         nspinb = SIZE(rho_ao_b)
         nspins = MIN(nspina, nspinb)
         DO i = 1, nspins
            tot_rho_r_a(i) = alpha*tot_rho_r_a(i) + beta*tot_rho_r_b(i)
         END DO
      END IF

      ! tot_rho_g
      IF (ASSOCIATED(tot_rho_g_a) .AND. ASSOCIATED(tot_rho_g_b)) THEN
         nspina = SIZE(rho_ao_a)
         nspinb = SIZE(rho_ao_b)
         nspins = MIN(nspina, nspinb)
         DO i = 1, nspins
            tot_rho_g_a(i) = alpha*tot_rho_g_a(i) + beta*tot_rho_g_b(i)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_rho_scale_and_add

! **************************************************************************************************
!> \brief Duplicates a pointer physically
!> \param rho_input The rho structure to be duplicated
!> \param rho_output The duplicate rho structure
!> \param qs_env The QS environment from which the auxiliary PW basis-set
!>                pool is taken
!> \par History
!>      07.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!> \note
!>      Associated pointers are deallocated, nullified pointers are NOT accepted!
! **************************************************************************************************
   SUBROUTINE duplicate_rho_type(rho_input, rho_output, qs_env)

      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_input, rho_output
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'duplicate_rho_type'

      INTEGER                                            :: handle, i, j, nspins
      LOGICAL :: complex_rho_ao_in, drho_g_valid_in, drho_r_valid_in, rho_g_valid_in, &
         rho_r_valid_in, soft_valid_in, tau_g_valid_in, tau_r_valid_in
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_g_in, tot_rho_g_out, &
                                                            tot_rho_r_in, tot_rho_r_out
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao_im_in, rho_ao_im_out, rho_ao_in, &
                                                            rho_ao_out
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type), DIMENSION(:), POINTER               :: rho_g_in, rho_g_out, rho_r_in, &
                                                            rho_r_out, tau_g_in, tau_g_out, &
                                                            tau_r_in, tau_r_out
      TYPE(pw_type), DIMENSION(:, :), POINTER            :: drho_g_in, drho_g_out, drho_r_in, &
                                                            drho_r_out
      TYPE(pw_type), POINTER                             :: rho_r_sccs_in, rho_r_sccs_out

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, pw_env, auxbas_pw_pool)
      NULLIFY (rho_ao_in, rho_ao_out, rho_ao_im_in, rho_ao_im_out)
      NULLIFY (rho_r_in, rho_r_out, rho_g_in, rho_g_out, drho_r_in, drho_r_out)
      NULLIFY (drho_g_in, drho_g_out, tau_r_in, tau_r_out, tau_g_in, tau_g_out)
      NULLIFY (tot_rho_r_in, tot_rho_r_out, tot_rho_g_in, tot_rho_g_out)
      NULLIFY (rho_r_sccs_in, rho_r_sccs_out)

      CPASSERT(ASSOCIATED(qs_env))

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, dft_control=dft_control)
      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool)
      nspins = dft_control%nspins

      CALL qs_rho_clear(rho_output)

      CALL qs_rho_get(rho_input, &
                      rho_ao=rho_ao_in, &
                      rho_ao_im=rho_ao_im_in, &
                      rho_r=rho_r_in, &
                      rho_g=rho_g_in, &
                      drho_r=drho_r_in, &
                      drho_g=drho_g_in, &
                      tau_r=tau_r_in, &
                      tau_g=tau_g_in, &
                      tot_rho_r=tot_rho_r_in, &
                      tot_rho_g=tot_rho_g_in, &
                      rho_g_valid=rho_g_valid_in, &
                      rho_r_valid=rho_r_valid_in, &
                      drho_g_valid=drho_g_valid_in, &
                      drho_r_valid=drho_r_valid_in, &
                      tau_r_valid=tau_r_valid_in, &
                      tau_g_valid=tau_g_valid_in, &
                      rho_r_sccs=rho_r_sccs_in, &
                      soft_valid=soft_valid_in, &
                      complex_rho_ao=complex_rho_ao_in)

      ! rho_ao
      IF (ASSOCIATED(rho_ao_in)) THEN
         CALL dbcsr_allocate_matrix_set(rho_ao_out, nspins)
         CALL qs_rho_set(rho_output, rho_ao=rho_ao_out)
         DO i = 1, nspins
            ALLOCATE (rho_ao_out(i)%matrix)
            CALL dbcsr_copy(rho_ao_out(i)%matrix, rho_ao_in(i)%matrix, &
                            name="myDensityMatrix_for_Spin_"//TRIM(ADJUSTL(cp_to_string(i))))
            CALL dbcsr_set(rho_ao_out(i)%matrix, 0.0_dp)
         END DO
      END IF

      ! rho_ao_im
      IF (ASSOCIATED(rho_ao_im_in)) THEN
         CALL dbcsr_allocate_matrix_set(rho_ao_im_out, nspins)
         CALL qs_rho_set(rho_output, rho_ao=rho_ao_im_out)
         DO i = 1, nspins
            ALLOCATE (rho_ao_im_out(i)%matrix)
            CALL dbcsr_copy(rho_ao_im_out(i)%matrix, rho_ao_im_in(i)%matrix, &
                            name="myImagDensityMatrix_for_Spin_"//TRIM(ADJUSTL(cp_to_string(i))))
            CALL dbcsr_set(rho_ao_im_out(i)%matrix, 0.0_dp)
         END DO
      END IF

      ! rho_r
      IF (ASSOCIATED(rho_r_in)) THEN
         ALLOCATE (rho_r_out(nspins))
         CALL qs_rho_set(rho_output, rho_r=rho_r_out)
         DO i = 1, nspins
            CALL auxbas_pw_pool%create_pw(rho_r_out(i), &
                                          use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_copy(rho_r_in(i), rho_r_out(i))
         END DO
      END IF

      ! rho_g
      IF (ASSOCIATED(rho_g_in)) THEN
         ALLOCATE (rho_g_out(nspins))
         CALL qs_rho_set(rho_output, rho_g=rho_g_out)
         DO i = 1, nspins
            CALL auxbas_pw_pool%create_pw(rho_g_out(i), &
                                          use_data=COMPLEXDATA1D, &
                                          in_space=RECIPROCALSPACE)
            CALL pw_copy(rho_g_in(i), rho_g_out(i))
         END DO
      END IF

      ! SCCS
      IF (ASSOCIATED(rho_r_sccs_in)) THEN
         CALL qs_rho_set(rho_output, rho_r_sccs=rho_r_sccs_out)
         CALL auxbas_pw_pool%create_pw(rho_r_sccs_out, &
                                       in_space=REALSPACE, &
                                       use_data=REALDATA3D)
         CALL pw_copy(rho_r_sccs_in, rho_r_sccs_out)
      END IF

      ! drho_r and drho_g are only needed if calculated by collocation
      IF (dft_control%drho_by_collocation) THEN
         ! drho_r
         IF (ASSOCIATED(drho_r_in)) THEN
            ALLOCATE (drho_r_out(3, nspins))
            CALL qs_rho_set(rho_output, drho_r=drho_r_out)
            DO j = 1, nspins
               DO i = 1, 3
                  CALL auxbas_pw_pool%create_pw(drho_r_out(i, j), &
                                                use_data=REALDATA3D, in_space=REALSPACE)
                  CALL pw_copy(drho_r_in(i, j), drho_r_out(i, j))
               END DO
            END DO
         END IF

         ! drho_g
         IF (ASSOCIATED(drho_g_in)) THEN
            ALLOCATE (drho_g_out(3, nspins))
            CALL qs_rho_set(rho_output, drho_g=drho_g_out)
            DO j = 1, nspins
               DO i = 1, 3
                  CALL auxbas_pw_pool%create_pw(drho_g_out(i, j), &
                                                use_data=COMPLEXDATA1D, &
                                                in_space=RECIPROCALSPACE)
                  CALL pw_copy(drho_g_in(i, j), drho_g_out(i, j))
               END DO
            END DO
         END IF
      END IF

      ! tau_r and tau_g are only needed in the case of Meta-GGA XC-functionals
      ! are used. Therefore they are only allocated if
      ! dft_control%use_kinetic_energy_density is true
      IF (dft_control%use_kinetic_energy_density) THEN
         ! tau_r
         IF (ASSOCIATED(tau_r_in)) THEN
            ALLOCATE (tau_r_out(nspins))
            CALL qs_rho_set(rho_output, tau_r=tau_r_out)
            DO i = 1, nspins
               CALL auxbas_pw_pool%create_pw(tau_r_out(i), &
                                             use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_copy(tau_r_in(i), tau_r_out(i))
            END DO
         END IF

         ! tau_g
         IF (ASSOCIATED(tau_g_in)) THEN
            ALLOCATE (tau_g_out(nspins))
            CALL qs_rho_set(rho_output, tau_g=tau_g_out)
            DO i = 1, nspins
               CALL auxbas_pw_pool%create_pw(tau_g_out(i), &
                                             use_data=COMPLEXDATA1D, &
                                             in_space=RECIPROCALSPACE)
               CALL pw_copy(tau_g_in(i), tau_g_out(i))
            END DO
         END IF
      END IF

      CALL qs_rho_set(rho_output, &
                      rho_g_valid=rho_g_valid_in, &
                      rho_r_valid=rho_r_valid_in, &
                      drho_g_valid=drho_g_valid_in, &
                      drho_r_valid=drho_r_valid_in, &
                      tau_r_valid=tau_r_valid_in, &
                      tau_g_valid=tau_g_valid_in, &
                      soft_valid=soft_valid_in, &
                      complex_rho_ao=complex_rho_ao_in)

      ! tot_rho_r
      IF (ASSOCIATED(tot_rho_r_in)) THEN
         ALLOCATE (tot_rho_r_out(nspins))
         CALL qs_rho_set(rho_output, tot_rho_r=tot_rho_r_out)
         DO i = 1, nspins
            tot_rho_r_out(i) = tot_rho_r_in(i)
         END DO
      END IF

      ! tot_rho_g
      IF (ASSOCIATED(tot_rho_g_in)) THEN
         ALLOCATE (tot_rho_g_out(nspins))
         CALL qs_rho_set(rho_output, tot_rho_g=tot_rho_g_out)
         DO i = 1, nspins
            tot_rho_g_out(i) = tot_rho_g_in(i)
         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE duplicate_rho_type

! **************************************************************************************************
!> \brief (Re-)allocates rho_ao_im from real part rho_ao
!> \param rho ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE allocate_rho_ao_imag_from_real(rho, qs_env)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=default_string_length)               :: headline
      INTEGER                                            :: i, ic, nimages, nspins
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_im_kp, rho_ao_kp
      TYPE(dbcsr_type), POINTER                          :: template
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      NULLIFY (rho_ao_im_kp, rho_ao_kp, dft_control, template, sab_orb)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      sab_orb=sab_orb)

      CALL qs_rho_get(rho, rho_ao_im_kp=rho_ao_im_kp, rho_ao_kp=rho_ao_kp)

      nspins = dft_control%nspins
      nimages = dft_control%nimages

      CPASSERT(nspins .EQ. SIZE(rho_ao_kp, 1))
      CPASSERT(nimages .EQ. SIZE(rho_ao_kp, 2))

      CALL dbcsr_allocate_matrix_set(rho_ao_im_kp, nspins, nimages)
      CALL qs_rho_set(rho, rho_ao_im_kp=rho_ao_im_kp)
      DO i = 1, nspins
         DO ic = 1, nimages
            IF (nspins > 1) THEN
               IF (i == 1) THEN
                  headline = "IMAGINARY PART OF DENSITY MATRIX FOR ALPHA SPIN"
               ELSE
                  headline = "IMAGINARY PART OF DENSITY MATRIX FOR BETA SPIN"
               END IF
            ELSE
               headline = "IMAGINARY PART OF DENSITY MATRIX"
            END IF
            ALLOCATE (rho_ao_im_kp(i, ic)%matrix)
            template => rho_ao_kp(i, ic)%matrix ! base on real part, but anti-symmetric
            CALL dbcsr_create(matrix=rho_ao_im_kp(i, ic)%matrix, template=template, &
                              name=TRIM(headline), matrix_type=dbcsr_type_antisymmetric, nze=0)
            CALL cp_dbcsr_alloc_block_from_nbl(rho_ao_im_kp(i, ic)%matrix, sab_orb)
            CALL dbcsr_set(rho_ao_im_kp(i, ic)%matrix, 0.0_dp)
         END DO
      END DO

   END SUBROUTINE allocate_rho_ao_imag_from_real

END MODULE qs_rho_methods
