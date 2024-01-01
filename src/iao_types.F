!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate ntrinsic atomic orbitals and analyze wavefunctions
!> \par History
!>      03.2023 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE iao_types
   USE cell_types,                      ONLY: cell_type
   USE input_constants,                 ONLY: do_iaoloc_enone,&
                                              do_iaoloc_pm2
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC ::  iao_env_type, iao_read_input, iao_set_default

! **************************************************************************************************
   TYPE iao_env_type
      LOGICAL                                :: do_iao = .FALSE.
      !
      REAL(KIND=dp)                          :: eps_svd = 0.0_dp
      REAL(KIND=dp)                          :: eps_occ = 0.0_dp
      ! chages
      LOGICAL                                :: do_charges = .FALSE.
      ! one-center expansion
      LOGICAL                                :: do_oce = .FALSE.
      INTEGER                                :: lmax_oce = 0
      INTEGER                                :: nbas_oce = 0
      LOGICAL                                :: append_oce = .FALSE.
      ! Bond orbitals
      LOGICAL                                :: do_bondorbitals = .FALSE.
      ! Wannier centers
      LOGICAL                                :: do_center = .FALSE.
      LOGICAL                                :: pos_periodic = .FALSE.
      INTEGER                                :: loc_operator = 0
      INTEGER                                :: eloc_function = 0
      REAL(KIND=dp)                          :: eloc_weight = 0.0_dp
      ! Molden
      LOGICAL                                :: molden_iao = .FALSE.
      LOGICAL                                :: molden_ibo = .FALSE.
      ! CUBE files
      LOGICAL                                :: cubes_iao = .FALSE.
      LOGICAL                                :: cubes_ibo = .FALSE.
      ! Input sections
      TYPE(section_vals_type), POINTER       :: iao_cubes_section => NULL(), &
                                                iao_molden_section => NULL(), &
                                                ibo_cubes_section => NULL(), &
                                                ibo_molden_section => NULL(), &
                                                ibo_cc_section => NULL()
   END TYPE iao_env_type

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param iao_env ...
! **************************************************************************************************
   SUBROUTINE iao_set_default(iao_env)
      TYPE(iao_env_type), INTENT(INOUT)                  :: iao_env

      !iao
      iao_env%do_iao = .FALSE.
      iao_env%eps_svd = 0.0_dp
      iao_env%eps_occ = 0.0_dp
      ! charges
      iao_env%do_charges = .FALSE.
      ! one-center expansion
      iao_env%do_oce = .FALSE.
      iao_env%lmax_oce = 3
      iao_env%nbas_oce = 10
      iao_env%append_oce = .FALSE.
      ! Bond orbitals
      iao_env%do_bondorbitals = .FALSE.
      ! Wannier centers
      iao_env%do_center = .FALSE.
      iao_env%pos_periodic = .FALSE.
      iao_env%loc_operator = do_iaoloc_pm2
      iao_env%eloc_function = do_iaoloc_enone
      iao_env%eloc_weight = 0.0_dp
      ! i/o
      iao_env%molden_iao = .FALSE.
      iao_env%molden_ibo = .FALSE.
      iao_env%cubes_iao = .FALSE.
      iao_env%cubes_ibo = .FALSE.
      ! Input sections
      NULLIFY (iao_env%iao_cubes_section, iao_env%iao_molden_section)
      NULLIFY (iao_env%ibo_cubes_section, iao_env%ibo_molden_section)
      NULLIFY (iao_env%ibo_cc_section)

   END SUBROUTINE iao_set_default

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param iao_env ...
!> \param iao_section ...
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE iao_read_input(iao_env, iao_section, cell)
      TYPE(iao_env_type), INTENT(INOUT)                  :: iao_env
      TYPE(section_vals_type), POINTER                   :: iao_section
      TYPE(cell_type), OPTIONAL                          :: cell

      LOGICAL                                            :: explicit, iao_explicit
      TYPE(section_vals_type), POINTER                   :: subsection

      CALL iao_set_default(iao_env)

      CALL section_vals_get(iao_section, explicit=iao_explicit)
      IF (iao_explicit) THEN
         iao_env%do_iao = .TRUE.
         ! input options
         CALL section_vals_val_get(iao_section, "EPS_SVD", r_val=iao_env%eps_svd)
         CALL section_vals_val_get(iao_section, "EPS_OCC", r_val=iao_env%eps_occ)
         CALL section_vals_val_get(iao_section, "ATOMIC_CHARGES", l_val=iao_env%do_charges)
         iao_env%iao_molden_section => section_vals_get_subs_vals(iao_section, "IAO_MOLDEN")
         CALL section_vals_get(iao_env%iao_molden_section, explicit=iao_env%molden_iao)
         iao_env%iao_cubes_section => section_vals_get_subs_vals(iao_section, "IAO_CUBES")
         CALL section_vals_get(iao_env%iao_cubes_section, explicit=iao_env%cubes_iao)
         subsection => section_vals_get_subs_vals(iao_section, "ONE_CENTER_EXPANSION")
         CALL section_vals_get(subsection, explicit=iao_env%do_oce)
         IF (iao_env%do_oce) THEN
            subsection => section_vals_get_subs_vals(iao_section, "ONE_CENTER_EXPANSION")
            CALL section_vals_val_get(subsection, "LMAX", i_val=iao_env%lmax_oce)
            CALL section_vals_val_get(subsection, "NBAS", i_val=iao_env%nbas_oce)
            CALL section_vals_val_get(subsection, "APPEND", l_val=iao_env%append_oce)
         END IF
         subsection => section_vals_get_subs_vals(iao_section, "BOND_ORBITALS")
         CALL section_vals_get(subsection, explicit=iao_env%do_bondorbitals)
         IF (iao_env%do_bondorbitals) THEN
            subsection => section_vals_get_subs_vals(iao_section, "BOND_ORBITALS")
            CALL section_vals_val_get(subsection, "LOCALIZATION_OPERATOR", i_val=iao_env%loc_operator)
            CALL section_vals_val_get(subsection, "ENERGY_LOCALIZATION_FUNCTION", i_val=iao_env%eloc_function)
            CALL section_vals_val_get(subsection, "ENERGY_LOCALIZATION_WEIGHT", r_val=iao_env%eloc_weight)
            iao_env%ibo_molden_section => section_vals_get_subs_vals(subsection, "IBO_MOLDEN")
            CALL section_vals_get(iao_env%ibo_molden_section, explicit=iao_env%molden_ibo)
            iao_env%ibo_cubes_section => section_vals_get_subs_vals(subsection, "IBO_CUBES")
            CALL section_vals_get(iao_env%ibo_cubes_section, explicit=iao_env%cubes_ibo)
            iao_env%ibo_cc_section => section_vals_get_subs_vals(subsection, "CHARGE_CENTER")
            CALL section_vals_get(iao_env%ibo_cc_section, explicit=iao_env%do_center)
            IF (iao_env%do_center) THEN
               CALL section_vals_val_get(iao_env%ibo_cc_section, "POSITION_OPERATOR_BERRY", &
                                         l_val=iao_env%pos_periodic, explicit=explicit)
               IF (.NOT. explicit) THEN
                  ! set default according to cell periodicity
                  iao_env%pos_periodic = .TRUE.
                  IF (PRESENT(cell)) THEN
                     IF (ALL(cell%perd == 0)) iao_env%pos_periodic = .FALSE.
                  END IF
               END IF
            END IF
         END IF
      END IF

   END SUBROUTINE iao_read_input

! **************************************************************************************************

END MODULE iao_types
