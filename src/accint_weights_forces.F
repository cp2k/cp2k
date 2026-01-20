!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!> \author JGH (01.2026)
! **************************************************************************************************
MODULE accint_weights_forces
   USE ao_util,                         ONLY: exp_radius_very_extended
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE grid_api,                        ONLY: integrate_pgf_product
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE memory_utilities,                ONLY: reallocate
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_scale
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_type
   USE qs_vxc,                          ONLY: qs_xc_energy
   USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                              transfer_pw2rs
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'accint_weights_forces'

   PUBLIC :: accint_weight_force

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho ...
!> \param xc_section ...
! **************************************************************************************************
   SUBROUTINE accint_weight_force(qs_env, rho, xc_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: xc_section

      CHARACTER(len=*), PARAMETER :: routineN = 'accint_weight_force'

      INTEGER                                            :: atom_a, handle, iatom, ikind, natom, &
                                                            natom_of_kind, nkind
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: use_virial
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: calpha, cvalue
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aforce
      REAL(KIND=dp), DIMENSION(3, 3)                     :: avirial
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: e_rspace
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)

      IF (dft_control%qs_control%gapw_control%accurate_xcint) THEN

         CALL get_qs_env(qs_env=qs_env, force=force, virial=virial)
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

         CALL get_qs_env(qs_env, natom=natom, nkind=nkind)
         ALLOCATE (aforce(3, natom))
         ALLOCATE (calpha(nkind), cvalue(nkind))
         cvalue = 1.0_dp
         calpha(1:nkind) = dft_control%qs_control%gapw_control%aw(1:nkind)

         CALL get_qs_env(qs_env, ks_env=ks_env, pw_env=pw_env)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL auxbas_pw_pool%create_pw(e_rspace)

         CALL qs_xc_energy(ks_env, rho, xc_section, e_rspace)
         CALL pw_scale(e_rspace, e_rspace%pw_grid%dvol)

         CALL gauss_grid_force(e_rspace, qs_env, calpha, cvalue, aforce, avirial)

         CALL auxbas_pw_pool%give_back_pw(e_rspace)

         CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom_of_kind, atom_list=atom_list)
            DO iatom = 1, natom_of_kind
               atom_a = atom_list(iatom)
               force(ikind)%rho_elec(1:3, iatom) = &
                  force(ikind)%rho_elec(1:3, iatom) + aforce(1:3, atom_a)
            END DO
         END DO
         IF (use_virial) THEN
            virial%pv_exc = virial%pv_exc + avirial
            virial%pv_virial = virial%pv_virial + avirial
         END IF

         DEALLOCATE (aforce)

      END IF

      CALL timestop(handle)

   END SUBROUTINE accint_weight_force

! **************************************************************************************************
!> \brief computes the forces/virial due to atomic centered Gaussian functions
!> \param e_rspace Energy density
!> \param qs_env ...
!> \param calpha ...
!> \param cvalue ...
!> \param aforce ...
!> \param avirial ...
! **************************************************************************************************
   SUBROUTINE gauss_grid_force(e_rspace, qs_env, calpha, cvalue, aforce, avirial)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: e_rspace
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: calpha, cvalue
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: aforce
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: avirial

      CHARACTER(len=*), PARAMETER                        :: routineN = 'gauss_grid_force'

      INTEGER                                            :: atom_a, handle, iatom, ikind, j, &
                                                            natom_of_kind, npme
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: alpha, eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, ra
      REAL(KIND=dp), DIMENSION(3, 3)                     :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hab, pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(realspace_grid_type), POINTER                 :: rs_v

      CALL timeset(routineN, handle)

      ALLOCATE (cores(1))
      ALLOCATE (hab(1, 1))
      ALLOCATE (pab(1, 1))

      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_v)

      CALL transfer_pw2rs(rs_v, e_rspace)

      CALL get_qs_env(qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set)

      use_virial = .TRUE.
      avirial = 0.0_dp
      aforce = 0.0_dp

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      DO ikind = 1, SIZE(atomic_kind_set)

         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom_of_kind, atom_list=atom_list)

         alpha = calpha(ikind)
         pab(1, 1) = -cvalue(ikind)
         IF (alpha == 0.0_dp .OR. pab(1, 1) == 0.0_dp) CYCLE

         CALL reallocate(cores, 1, natom_of_kind)
         npme = 0
         cores = 0

         DO iatom = 1, natom_of_kind
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r, cell)
            IF (rs_v%desc%parallel .AND. .NOT. rs_v%desc%distributed) THEN
               ! replicated realspace grid, split the atoms up between procs
               IF (MODULO(iatom, rs_v%desc%group_size) == rs_v%desc%my_pos) THEN
                  npme = npme + 1
                  cores(npme) = iatom
               END IF
            ELSE
               npme = npme + 1
               cores(npme) = iatom
            END IF
         END DO

         DO j = 1, npme

            iatom = cores(j)
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r, cell)
            hab(1, 1) = 0.0_dp
            force_a(:) = 0.0_dp
            force_b(:) = 0.0_dp
            my_virial_a = 0.0_dp
            my_virial_b = 0.0_dp

            radius = exp_radius_very_extended(la_min=0, la_max=0, lb_min=0, lb_max=0, &
                                              ra=ra, rb=ra, rp=ra, &
                                              zetp=alpha, eps=eps_rho_rspace, &
                                              pab=pab, o1=0, o2=0, &
                                              prefactor=1.0_dp, cutoff=1.0_dp)

            CALL integrate_pgf_product(0, alpha, 0, &
                                       0, 0.0_dp, 0, ra, [0.0_dp, 0.0_dp, 0.0_dp], &
                                       rs_v, hab, pab=pab, o1=0, o2=0, &
                                       radius=radius, &
                                       calculate_forces=.TRUE., force_a=force_a, &
                                       force_b=force_b, use_virial=use_virial, my_virial_a=my_virial_a, &
                                       my_virial_b=my_virial_b, use_subpatch=.TRUE., subpatch_pattern=0)

            aforce(1:3, atom_a) = aforce(1:3, atom_a) + force_a(1:3)
            avirial = avirial + my_virial_a

         END DO

      END DO

      DEALLOCATE (hab, pab, cores)

      CALL timestop(handle)

   END SUBROUTINE gauss_grid_force

END MODULE accint_weights_forces
