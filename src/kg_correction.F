!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for a Kim-Gordon-like partitioning into molecular subunits
!> \par History
!>       2012.06 created [Martin Haeufel]
!> \author Martin Haeufel and Florian Schiffmann
! **************************************************************************************************
MODULE kg_correction
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_dot,&
                                              dbcsr_p_type
   USE input_constants,                 ONLY: kg_tnadd_atomic,&
                                              kg_tnadd_embed,&
                                              kg_tnadd_embed_ri,&
                                              kg_tnadd_none
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kg_environment_types,            ONLY: kg_environment_type
   USE kinds,                           ONLY: dp
   USE lri_environment_methods,         ONLY: calculate_lri_densities,&
                                              lri_kg_rho_update
   USE lri_environment_types,           ONLY: lri_density_type,&
                                              lri_environment_type,&
                                              lri_kind_type
   USE lri_forces,                      ONLY: calculate_lri_forces
   USE lri_ks_methods,                  ONLY: calculate_lri_ks_matrix
   USE message_passing,                 ONLY: mp_sum
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace,&
                                              integrate_v_rspace_one_center
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild,&
                                              qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE virial_types,                    ONLY: virial_type
   USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                              xc_prep_2nd_deriv
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_release
   USE xc_rho_set_types,                ONLY: xc_rho_set_release,&
                                              xc_rho_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_correction'

   PUBLIC :: kg_ekin_subset

CONTAINS

! **************************************************************************************************
!> \brief Calculates the subsystem Hohenberg-Kohn kinetic energy and the forces
!> \param qs_env ...
!> \param ks_matrix ...
!> \param ekin_mol ...
!> \param calc_force ...
!> \param do_kernel Contribution of kinetic energy functional to kernel in response calculation
!> \param pmat_ext Response density used to fold 2nd deriv or to integrate kinetic energy functional
!> \param alpha Multiplicative factor required in linear response equations
!> \par History
!>       2012.06 created [Martin Haeufel]
!>       2014.01 added atomic potential option [JGH]
!>       2020.01 Added KG contribution to linear response [fbelle]
!> \author Martin Haeufel and Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE kg_ekin_subset(qs_env, ks_matrix, ekin_mol, calc_force, do_kernel, pmat_ext, alpha)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_matrix
      REAL(KIND=dp), INTENT(out)                         :: ekin_mol
      LOGICAL, INTENT(IN)                                :: calc_force, do_kernel
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: pmat_ext
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: alpha

      LOGICAL                                            :: lrigpw
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kg_environment_type), POINTER                 :: kg_env

      CALL get_qs_env(qs_env, kg_env=kg_env, dft_control=dft_control)
      lrigpw = dft_control%qs_control%lrigpw

      IF (kg_env%tnadd_method == kg_tnadd_embed) THEN
         IF (lrigpw) THEN
            CALL kg_ekin_embed_lri(qs_env, kg_env, ks_matrix, ekin_mol, calc_force)
         ELSE
            CALL kg_ekin_embed(qs_env, kg_env, ks_matrix, ekin_mol, calc_force, &
                               do_kernel, pmat_ext, alpha)
         END IF
      ELSE IF (kg_env%tnadd_method == kg_tnadd_embed_ri) THEN
         CALL kg_ekin_ri_embed(qs_env, kg_env, ks_matrix, ekin_mol, calc_force)
      ELSE IF (kg_env%tnadd_method == kg_tnadd_atomic) THEN
         CALL kg_ekin_atomic(qs_env, ks_matrix, ekin_mol)
      ELSE IF (kg_env%tnadd_method == kg_tnadd_none) THEN
         ekin_mol = 0.0_dp
      ELSE
         CPABORT("Unknown KG embedding method")
      END IF

   END SUBROUTINE kg_ekin_subset

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kg_env ...
!> \param ks_matrix ...
!> \param ekin_mol ...
!> \param calc_force ...
!> \param do_kernel Contribution of kinetic energy functional to kernel in response calculation
!> \param pmat_ext Response density used to fold 2nd deriv or to integrate kinetic energy functional
!> \param alpha Multiplicative factor required in linear response equations
! **************************************************************************************************
   SUBROUTINE kg_ekin_embed(qs_env, kg_env, ks_matrix, ekin_mol, calc_force, do_kernel, pmat_ext, alpha)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_matrix
      REAL(KIND=dp), INTENT(out)                         :: ekin_mol
      LOGICAL, INTENT(IN)                                :: calc_force, do_kernel
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: pmat_ext
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: alpha

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'kg_ekin_embed'

      INTEGER                                            :: handle, ispin, isub, natom, nspins
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: ekin_imol, my_alpha
      REAL(KIND=dp), DIMENSION(3, 3)                     :: xcvirial
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: density_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho1_g, rho1_r, rho_r, vxc_rho, vxc_tau
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: old_rho, rho1, rho_struct
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (vxc_rho, vxc_tau, old_rho, rho_struct, ks_env)

      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      rho=old_rho, &
                      natom=natom, &
                      dft_control=dft_control, &
                      virial=virial, &
                      pw_env=pw_env)
      nspins = dft_control%nspins
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      use_virial = use_virial .AND. calc_force

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      ! get the density matrix
      CALL qs_rho_get(old_rho, rho_ao=density_matrix)
      ! allocate and initialize the density
      CALL qs_rho_create(rho_struct)
      ! set the density matrix to the blocked matrix
      CALL qs_rho_set(rho_struct, rho_ao=density_matrix) ! blocked_matrix
      CALL qs_rho_rebuild(rho_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)
      ! full density kinetic energy term
      CALL qs_rho_update_rho(rho_struct, qs_env)
      ! get blocked density that has been put on grid
      CALL qs_rho_get(rho_struct, rho_r=rho_r)

      ! If external density associated then it is needed either for
      ! 1) folding of second derivative while partially integrating, or
      ! 2) integration of response forces
      NULLIFY (rho1, rho1_r, rho1_g)
      IF (PRESENT(pmat_ext)) THEN
         CALL qs_rho_create(rho1)
         CALL qs_rho_set(rho1, rho_ao=pmat_ext)
         CALL qs_rho_rebuild(rho1, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)
         CALL qs_rho_update_rho(rho1, qs_env)
         CALL qs_rho_get(rho1, rho_r=rho1_r, rho_g=rho1_g)
      END IF

      ! Response solver requires factor (spin) for kernel matrix
      ! alpha = 2 closed-shell
      ! alpha = 1 open-shell
      my_alpha = 1.0_dp
      IF (PRESENT(alpha)) my_alpha = alpha

      ekin_imol = 0.0_dp

      ! calculate xc potential or kernel
      IF (do_kernel) THEN
         ! derivation wrt to rho_struct and evaluation at rho_struct
         CALL create_kernel(qs_env, &
                            vxc=vxc_rho, &
                            rho=rho_struct, &
                            rho1=rho1, &
                            xc_section=kg_env%xc_section_kg)
      ELSE
         CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_struct, xc_section=kg_env%xc_section_kg, &
                            vxc_rho=vxc_rho, vxc_tau=vxc_tau, exc=ekin_imol)
         IF (ASSOCIATED(vxc_tau)) THEN
            CPABORT(" KG with meta-kinetic energy functionals not implemented")
         END IF
      END IF

      DO ispin = 1, nspins
         vxc_rho(ispin)%pw%cr3d = my_alpha*vxc_rho(ispin)%pw%cr3d*vxc_rho(ispin)%pw%pw_grid%dvol
         ! Integrate xc-potential with external density for outer response forces
         IF (PRESENT(pmat_ext) .AND. .NOT. do_kernel) THEN
            CALL qs_rho_get(rho1, rho_ao=density_matrix)
         END IF
         CALL integrate_v_rspace(v_rspace=vxc_rho(ispin), &
                                 pmat=density_matrix(ispin), hmat=ks_matrix(ispin), &
                                 qs_env=qs_env, calculate_forces=calc_force)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_rho(ispin)%pw)
      END DO
      DEALLOCATE (vxc_rho)
      ekin_mol = -ekin_imol
      xcvirial(1:3, 1:3) = 0.0_dp
      IF (use_virial) xcvirial(1:3, 1:3) = xcvirial(1:3, 1:3) - virial%pv_xc(1:3, 1:3)

      ! loop over all subsets
      DO isub = 1, kg_env%nsubsets
         ! calculate the densities for the given blocked density matrix - pass the subset task_list
         CALL qs_rho_update_rho(rho_struct, qs_env, &
                                task_list_external=kg_env%subset(isub)%task_list)
         ! Same for external (response) density if present
         IF (PRESENT(pmat_ext)) THEN
            CALL qs_rho_update_rho(rho1, qs_env, &
                                   task_list_external=kg_env%subset(isub)%task_list)
         END IF

         ekin_imol = 0.0_dp
         ! calc Hohenberg-Kohn kin. energy of the density corresp. to the remaining molecular block(s)
         ! info per block in rho_struct now

         ! calculate kernel
         IF (do_kernel) THEN
            CALL create_kernel(qs_env, &
                               vxc=vxc_rho, &
                               rho=rho_struct, &
                               rho1=rho1, &
                               xc_section=kg_env%xc_section_kg)
         ELSE
            CALL qs_vxc_create(ks_env=ks_env, &
                               rho_struct=rho_struct, &
                               xc_section=kg_env%xc_section_kg, &
                               vxc_rho=vxc_rho, &
                               vxc_tau=vxc_tau, &
                               exc=ekin_imol)
         END IF

         ekin_mol = ekin_mol + ekin_imol

         DO ispin = 1, nspins
            vxc_rho(ispin)%pw%cr3d = -my_alpha*vxc_rho(ispin)%pw%cr3d*vxc_rho(ispin)%pw%pw_grid%dvol
            ! Integrate with response density for outer response forces
            IF (PRESENT(pmat_ext) .AND. .NOT. do_kernel) THEN
               CALL qs_rho_get(rho1, rho_ao=density_matrix)
            END IF
            CALL integrate_v_rspace(v_rspace=vxc_rho(ispin), &
                                    pmat=density_matrix(ispin), &
                                    hmat=ks_matrix(ispin), &
                                    qs_env=qs_env, &
                                    calculate_forces=calc_force, &
                                    task_list_external=kg_env%subset(isub)%task_list)
            ! clean up vxc_rho
            CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_rho(ispin)%pw)
         END DO
         DEALLOCATE (vxc_rho)

         IF (use_virial) THEN
            xcvirial(1:3, 1:3) = xcvirial(1:3, 1:3) + virial%pv_xc(1:3, 1:3)
         END IF

      END DO

      IF (use_virial) THEN
         virial%pv_xc(1:3, 1:3) = xcvirial(1:3, 1:3)
      END IF

      ! clean up rho_struct
      CALL qs_rho_set(rho_struct, rho_ao=Null())
      CALL qs_rho_release(rho_struct)
      IF (PRESENT(pmat_ext)) THEN
         CALL qs_rho_set(rho1, rho_ao_kp=Null())
         CALL qs_rho_release(rho1)
      END IF

      CALL timestop(handle)

   END SUBROUTINE kg_ekin_embed

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kg_env ...
!> \param ks_matrix ...
!> \param ekin_mol ...
!> \param calc_force ...
! **************************************************************************************************
   SUBROUTINE kg_ekin_embed_lri(qs_env, kg_env, ks_matrix, ekin_mol, calc_force)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_matrix
      REAL(KIND=dp), INTENT(out)                         :: ekin_mol
      LOGICAL                                            :: calc_force

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'kg_ekin_embed_lri'

      INTEGER                                            :: color, handle, iatom, ikind, imol, &
                                                            ispin, isub, natom, nkind, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atomlist
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: ekin_imol
      REAL(KIND=dp), DIMENSION(3, 3)                     :: xcvirial
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: density_matrix, ksmat
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: pmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_v_int
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rho, vxc_tau
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: old_rho, rho_struct
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (vxc_rho, vxc_tau, old_rho, rho_struct, ks_env)

      CALL get_qs_env(qs_env, dft_control=dft_control)

      ! get set of molecules, natom, dft_control, pw_env
      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      rho=old_rho, &
                      natom=natom, &
                      dft_control=dft_control, &
                      virial=virial, &
                      para_env=para_env, &
                      pw_env=pw_env)

      nspins = dft_control%nspins
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      use_virial = use_virial .AND. calc_force

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      ! get the density matrix
      CALL qs_rho_get(old_rho, rho_ao=density_matrix)
      ! allocate and initialize the density
      CALL qs_rho_create(rho_struct)
      ! set the density matrix to the blocked matrix
      CALL qs_rho_set(rho_struct, rho_ao=density_matrix) ! blocked_matrix
      CALL qs_rho_rebuild(rho_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)

      CALL get_qs_env(qs_env, lri_env=lri_env, lri_density=lri_density, nkind=nkind)
      IF (lri_env%exact_1c_terms) THEN
         CPABORT(" KG with LRI and exact one-center terms not implemented")
      END IF
      ALLOCATE (atomlist(natom))
      DO ispin = 1, nspins
         lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
         DO ikind = 1, nkind
            lri_v_int(ikind)%v_int = 0.0_dp
            IF (calc_force) THEN
               lri_v_int(ikind)%v_dadr = 0.0_dp
               lri_v_int(ikind)%v_dfdr = 0.0_dp
            END IF
         END DO
      END DO

      ! full density kinetic energy term
      atomlist = 1
      CALL lri_kg_rho_update(rho_struct, qs_env, lri_env, lri_density, atomlist)
      ekin_imol = 0.0_dp
      CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_struct, xc_section=kg_env%xc_section_kg, &
                         vxc_rho=vxc_rho, vxc_tau=vxc_tau, exc=ekin_imol)
      IF (ASSOCIATED(vxc_tau)) THEN
         CPABORT(" KG with meta-kinetic energy functionals not implemented")
      END IF
      DO ispin = 1, nspins
         vxc_rho(ispin)%pw%cr3d = vxc_rho(ispin)%pw%cr3d*vxc_rho(ispin)%pw%pw_grid%dvol
         lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
         CALL integrate_v_rspace_one_center(vxc_rho(ispin), qs_env, lri_v_int, calc_force, "LRI_AUX")
         CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_rho(ispin)%pw)
      END DO
      DEALLOCATE (vxc_rho)
      ekin_mol = -ekin_imol
      xcvirial(1:3, 1:3) = 0.0_dp
      IF (use_virial) xcvirial(1:3, 1:3) = xcvirial(1:3, 1:3) - virial%pv_xc(1:3, 1:3)

      ! loop over all subsets
      DO isub = 1, kg_env%nsubsets
         atomlist = 0
         DO iatom = 1, natom
            imol = kg_env%atom_to_molecule(iatom)
            color = kg_env%subset_of_mol(imol)
            IF (color == isub) atomlist(iatom) = 1
         END DO
         CALL lri_kg_rho_update(rho_struct, qs_env, lri_env, lri_density, atomlist)

         ekin_imol = 0.0_dp
         ! calc Hohenberg-Kohn kin. energy of the density corresp. to the remaining molecular block(s)
         CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_struct, xc_section=kg_env%xc_section_kg, &
                            vxc_rho=vxc_rho, vxc_tau=vxc_tau, exc=ekin_imol)
         ekin_mol = ekin_mol + ekin_imol

         DO ispin = 1, nspins
            vxc_rho(ispin)%pw%cr3d = -vxc_rho(ispin)%pw%cr3d*vxc_rho(ispin)%pw%pw_grid%dvol
            lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
            CALL integrate_v_rspace_one_center(vxc_rho(ispin), qs_env, &
                                               lri_v_int, calc_force, &
                                               "LRI_AUX", atomlist=atomlist)
            ! clean up vxc_rho
            CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_rho(ispin)%pw)
         END DO
         DEALLOCATE (vxc_rho)

         IF (use_virial) THEN
            xcvirial(1:3, 1:3) = xcvirial(1:3, 1:3) + virial%pv_xc(1:3, 1:3)
         END IF

      END DO

      IF (use_virial) THEN
         virial%pv_xc(1:3, 1:3) = xcvirial(1:3, 1:3)
      END IF

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
      ALLOCATE (ksmat(1))
      DO ispin = 1, nspins
         lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
         DO ikind = 1, nkind
            CALL mp_sum(lri_v_int(ikind)%v_int, para_env%group)
         END DO
         ksmat(1)%matrix => ks_matrix(ispin)%matrix
         CALL calculate_lri_ks_matrix(lri_env, lri_v_int, ksmat, atomic_kind_set)
      END DO
      IF (calc_force) THEN
         pmat(1:nspins, 1:1) => density_matrix(1:nspins)
         CALL calculate_lri_forces(lri_env, lri_density, qs_env, pmat, atomic_kind_set)
      ENDIF
      DEALLOCATE (atomlist, ksmat)

      ! clean up rho_struct
      CALL qs_rho_set(rho_struct, rho_ao=Null())
      CALL qs_rho_release(rho_struct)

      CALL timestop(handle)

   END SUBROUTINE kg_ekin_embed_lri

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kg_env ...
!> \param ks_matrix ...
!> \param ekin_mol ...
!> \param calc_force ...
! **************************************************************************************************
   SUBROUTINE kg_ekin_ri_embed(qs_env, kg_env, ks_matrix, ekin_mol, calc_force)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_matrix
      REAL(KIND=dp), INTENT(out)                         :: ekin_mol
      LOGICAL                                            :: calc_force

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'kg_ekin_ri_embed'

      INTEGER                                            :: color, handle, iatom, ikind, imol, &
                                                            ispin, isub, natom, nkind, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atomlist
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: ekin_imol
      REAL(KIND=dp), DIMENSION(3, 3)                     :: xcvirial
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ksmat
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: density_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(lri_density_type), POINTER                    :: lri_density
      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_v_int
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rho, vxc_tau
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_struct
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      rho=rho, &
                      natom=natom, &
                      nkind=nkind, &
                      dft_control=dft_control, &
                      virial=virial, &
                      para_env=para_env, &
                      pw_env=pw_env)

      nspins = dft_control%nspins
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      use_virial = use_virial .AND. calc_force

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      ! get the density matrix
      CALL qs_rho_get(rho, rho_ao_kp=density_matrix)
      ! allocate and initialize the density
      NULLIFY (rho_struct)
      CALL qs_rho_create(rho_struct)
      ! set the density matrix to the blocked matrix
      CALL qs_rho_set(rho_struct, rho_ao_kp=density_matrix)
      CALL qs_rho_rebuild(rho_struct, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
      ALLOCATE (cell_to_index(1, 1, 1))
      cell_to_index(1, 1, 1) = 1
      lri_env => kg_env%lri_env
      lri_density => kg_env%lri_density
      CALL calculate_lri_densities(lri_env, lri_density, qs_env, density_matrix, cell_to_index, &
                                   rho_struct, atomic_kind_set, para_env)
      kg_env%lri_density => lri_density
      ! full density kinetic energy term
      ekin_imol = 0.0_dp
      NULLIFY (vxc_rho, vxc_tau)
      CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_struct, xc_section=kg_env%xc_section_kg, &
                         vxc_rho=vxc_rho, vxc_tau=vxc_tau, exc=ekin_imol)
      IF (ASSOCIATED(vxc_tau)) THEN
         CPABORT(" KG with meta-kinetic energy functionals not implemented")
      END IF
      DO ispin = 1, nspins
         vxc_rho(ispin)%pw%cr3d = vxc_rho(ispin)%pw%cr3d*vxc_rho(ispin)%pw%pw_grid%dvol
         lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
         CALL integrate_v_rspace_one_center(vxc_rho(ispin), qs_env, lri_v_int, calc_force, "LRI_AUX")
         CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_rho(ispin)%pw)
      END DO
      DEALLOCATE (vxc_rho)
      ekin_mol = -ekin_imol
      xcvirial(1:3, 1:3) = 0.0_dp
      IF (use_virial) xcvirial(1:3, 1:3) = xcvirial(1:3, 1:3) - virial%pv_xc(1:3, 1:3)
!deb
!     WRITE(6,*) " E KIN (full)       ",-ekin_mol
!deb

      ! loop over all subsets
      ALLOCATE (atomlist(natom))
      DO isub = 1, kg_env%nsubsets
         atomlist = 0
         DO iatom = 1, natom
            imol = kg_env%atom_to_molecule(iatom)
            color = kg_env%subset_of_mol(imol)
            IF (color == isub) atomlist(iatom) = 1
         END DO
         CALL lri_kg_rho_update(rho_struct, qs_env, lri_env, lri_density, atomlist)

         ekin_imol = 0.0_dp
         ! calc Hohenberg-Kohn kin. energy of the density corresp. to the remaining molecular block(s)
         NULLIFY (vxc_rho, vxc_tau)
         CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_struct, xc_section=kg_env%xc_section_kg, &
                            vxc_rho=vxc_rho, vxc_tau=vxc_tau, exc=ekin_imol)
         ekin_mol = ekin_mol + ekin_imol
!deb
!     WRITE(6,*) " E KIN (molecule)    ",isub,ekin_imol
!deb

         DO ispin = 1, nspins
            vxc_rho(ispin)%pw%cr3d = -vxc_rho(ispin)%pw%cr3d*vxc_rho(ispin)%pw%pw_grid%dvol
            lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
            CALL integrate_v_rspace_one_center(vxc_rho(ispin), qs_env, &
                                               lri_v_int, calc_force, &
                                               "LRI_AUX", atomlist=atomlist)
            ! clean up vxc_rho
            CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_rho(ispin)%pw)
         END DO
         DEALLOCATE (vxc_rho)

         IF (use_virial) THEN
            xcvirial(1:3, 1:3) = xcvirial(1:3, 1:3) + virial%pv_xc(1:3, 1:3)
         END IF

      END DO

      IF (use_virial) THEN
         virial%pv_xc(1:3, 1:3) = xcvirial(1:3, 1:3)
      END IF

      ALLOCATE (ksmat(1))
      DO ispin = 1, nspins
         lri_v_int => lri_density%lri_coefs(ispin)%lri_kinds
         DO ikind = 1, nkind
            CALL mp_sum(lri_v_int(ikind)%v_int, para_env%group)
         END DO
         ksmat(1)%matrix => ks_matrix(ispin)%matrix
         CALL calculate_lri_ks_matrix(lri_env, lri_v_int, ksmat, atomic_kind_set)
      END DO
      IF (calc_force) THEN
         CALL calculate_lri_forces(lri_env, lri_density, qs_env, density_matrix, atomic_kind_set)
      ENDIF
      DEALLOCATE (atomlist, ksmat)

      ! clean up rho_struct
      CALL qs_rho_set(rho_struct, rho_ao=Null())
      CALL qs_rho_release(rho_struct)
      DEALLOCATE (cell_to_index)

      CALL timestop(handle)

   END SUBROUTINE kg_ekin_ri_embed

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ks_matrix ...
!> \param ekin_mol ...
! **************************************************************************************************
   SUBROUTINE kg_ekin_atomic(qs_env, ks_matrix, ekin_mol)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_matrix
      REAL(KIND=dp), INTENT(out)                         :: ekin_mol

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'kg_ekin_atomic'

      INTEGER                                            :: handle, ispin, nspins
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: density_matrix, tnadd_matrix
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(qs_rho_type), POINTER                         :: rho

      NULLIFY (rho, kg_env, density_matrix, tnadd_matrix)

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, kg_env=kg_env, rho=rho)

      nspins = SIZE(ks_matrix)
      ! get the density matrix
      CALL qs_rho_get(rho, rho_ao=density_matrix)
      ! get the tnadd matrix
      tnadd_matrix => kg_env%tnadd_mat

      ekin_mol = 0.0_dp
      DO ispin = 1, nspins
         CALL dbcsr_dot(tnadd_matrix(1)%matrix, density_matrix(ispin)%matrix, ekin_mol)
         CALL dbcsr_add(ks_matrix(ispin)%matrix, tnadd_matrix(1)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END DO
      ! definition is inverted (see qs_ks_methods)
      ekin_mol = -ekin_mol

      CALL timestop(handle)

   END SUBROUTINE kg_ekin_atomic

! **************************************************************************************************
!> \brief   Creation of second derivative xc-potential
!> \param   qs_env ...
!> \param   vxc will contain the partially integrated second derivative
!>          taken with respect to rho, evaluated in rho and folded with rho1
!>          vxc is allocated here and needs to be deallocated by the caller.
!> \param rho density at which derivatives were calculated
!> \param rho1 density with which to fold
!> \param xc_section XC parameters
!> \date    11.2019
!> \author  fbelle
! **************************************************************************************************
   SUBROUTINE create_kernel(qs_env, vxc, rho, rho1, xc_section)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), DIMENSION(:), INTENT(OUT), &
         POINTER                                         :: vxc
      TYPE(qs_rho_type), POINTER                         :: rho, rho1
      TYPE(section_vals_type), POINTER                   :: xc_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_kernel'

      INTEGER                                            :: handle, ispin, nspins
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho1_g, rho1_r, rho_r, vxc_rho
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(xc_rho_set_type), POINTER                     :: rho_set

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(vxc))
      NULLIFY (auxbas_pw_pool, pw_env, vxc_rho)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      pw_env=pw_env)
      nspins = dft_control%nspins

      CALL pw_env_get(pw_env=pw_env, &
                      auxbas_pw_pool=auxbas_pw_pool)

      ! Get grid based density
      NULLIFY (rho_r, rho1_r, rho1_g)
      CALL qs_rho_get(rho, rho_r=rho_r)
      CALL qs_rho_get(rho1, rho_r=rho1_r, rho_g=rho1_g)

      NULLIFY (deriv_set, rho_set)

      ! main ingredient is xc_rho_set_and_dset_create
      CALL xc_prep_2nd_deriv(deriv_set=deriv_set, &    ! containing potentials
                             rho_set=rho_set, &        ! density at which derivs are calculated
                             rho_r=rho_r, &            ! place where derivative is evaluated
                             pw_pool=auxbas_pw_pool, & ! pool for grids
                             xc_section=xc_section)

      ! evaluation of 2nd deriv in rho_set density
      ! folding of second deriv with density in rho1_set
      CALL xc_calc_2nd_deriv(v_xc=vxc_rho, &           ! XC-potential
                             deriv_set=deriv_set, &    ! deriv of xc-potential
                             rho_set=rho_set, &        ! density at which deriv are calculated
                             rho1_r=rho1_r, &      ! density with which to fold
                             rho1_g=rho1_g, &
                             pw_pool=auxbas_pw_pool, & ! pool for grids
                             xc_section=xc_section, &
                             gapw=.FALSE.)

      ! Release second deriv stuff
      CALL xc_dset_release(deriv_set)
      CALL xc_rho_set_release(rho_set=rho_set, pw_pool=auxbas_pw_pool)

      ! export vxc
      ALLOCATE (vxc(nspins))
      DO ispin = 1, nspins
         vxc(ispin)%pw => vxc_rho(ispin)%pw
      END DO
      DEALLOCATE (vxc_rho)

      CALL timestop(handle)

   END SUBROUTINE

END MODULE kg_correction
