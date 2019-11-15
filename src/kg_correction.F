!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
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
   USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: pw_p_type
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
!> \par History
!>       2012.06 created [Martin Haeufel]
!>       2014.01 added atomic potential option [JGH]
!> \author Martin Haeufel and Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE kg_ekin_subset(qs_env, ks_matrix, ekin_mol, calc_force)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_matrix
      REAL(KIND=dp), INTENT(out)                         :: ekin_mol
      LOGICAL                                            :: calc_force

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_ekin_subset', routineP = moduleN//':'//routineN

      LOGICAL                                            :: lrigpw
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kg_environment_type), POINTER                 :: kg_env

      CALL get_qs_env(qs_env, kg_env=kg_env, dft_control=dft_control)
      lrigpw = dft_control%qs_control%lrigpw

      IF (kg_env%tnadd_method == kg_tnadd_embed) THEN
         IF (lrigpw) THEN
            CALL kg_ekin_embed_lri(qs_env, kg_env, ks_matrix, ekin_mol, calc_force)
         ELSE
            CALL kg_ekin_embed(qs_env, kg_env, ks_matrix, ekin_mol, calc_force)
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
! **************************************************************************************************
   SUBROUTINE kg_ekin_embed(qs_env, kg_env, ks_matrix, ekin_mol, calc_force)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_matrix
      REAL(KIND=dp), INTENT(out)                         :: ekin_mol
      LOGICAL                                            :: calc_force

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_ekin_embed', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, isub, natom, nspins
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: ekin_imol
      REAL(KIND=dp), DIMENSION(3, 3)                     :: xcvirial
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: density_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rho, vxc_tau
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: old_rho, rho_struct
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (vxc_rho, vxc_tau, old_rho, rho_struct, ks_env)

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

      ! full density kinetic energy term
      CALL qs_rho_update_rho(rho_struct, qs_env)
      ekin_imol = 0.0_dp
      CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_struct, xc_section=kg_env%xc_section_kg, &
                         vxc_rho=vxc_rho, vxc_tau=vxc_tau, exc=ekin_imol)
      IF (ASSOCIATED(vxc_tau)) THEN
         CPABORT(" KG with meta-kinetic energy functionals not implemented")
      END IF
      DO ispin = 1, nspins
         vxc_rho(ispin)%pw%cr3d = vxc_rho(ispin)%pw%cr3d*vxc_rho(ispin)%pw%pw_grid%dvol
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

         ekin_imol = 0.0_dp
         ! calc Hohenberg-Kohn kin. energy of the density corresp. to the remaining molecular block(s)
         CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_struct, xc_section=kg_env%xc_section_kg, &
                            vxc_rho=vxc_rho, vxc_tau=vxc_tau, exc=ekin_imol)
         ekin_mol = ekin_mol + ekin_imol

         DO ispin = 1, nspins
            vxc_rho(ispin)%pw%cr3d = -vxc_rho(ispin)%pw%cr3d*vxc_rho(ispin)%pw%pw_grid%dvol
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_ekin_embed_lri', &
         routineP = moduleN//':'//routineN

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_ekin_ri_embed', &
         routineP = moduleN//':'//routineN

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_ekin_atomic', routineP = moduleN//':'//routineN

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

END MODULE kg_correction
