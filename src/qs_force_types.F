!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      Add CP2K error reporting, new add_force routine [07.2014,JGH]
!> \author MK (03.06.2002)
! **************************************************************************************************
MODULE qs_force_types

   !USE cp_control_types,                ONLY: qs_control_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
#include "./base/base_uses.f90"

   IMPLICIT NONE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_force_types'
   PRIVATE

   TYPE qs_force_type
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: all_potential, &
                                                 core_overlap, &
                                                 gth_ppl, &
                                                 gth_nlcc, &
                                                 gth_ppnl, &
                                                 kinetic, &
                                                 overlap, &
                                                 overlap_admm, &
                                                 rho_core, &
                                                 rho_elec, &
                                                 rho_lri_elec, &
                                                 vhxc_atom, &
                                                 g0s_Vh_elec, &
                                                 repulsive, &
                                                 dispersion, &
                                                 gcp, &
                                                 other, &
                                                 ch_pulay, &
                                                 fock_4c, &
                                                 ehrenfest, &
                                                 efield, &
                                                 eev, &
                                                 mp2_non_sep, &
                                                 total
   END TYPE qs_force_type

   PUBLIC :: qs_force_type

   PUBLIC :: allocate_qs_force, &
             add_qs_force, &
             deallocate_qs_force, &
             replicate_qs_force, &
             sum_qs_force, &
             get_qs_force, &
             put_qs_force, &
             total_qs_force, &
             zero_qs_force

CONTAINS

! **************************************************************************************************
!> \brief   Allocate a Quickstep force data structure.
!> \param qs_force ...
!> \param natom_of_kind ...
!> \date    05.06.2002
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE allocate_qs_force(qs_force, natom_of_kind)

      TYPE(qs_force_type), DIMENSION(:), POINTER         :: qs_force
      INTEGER, DIMENSION(:), INTENT(IN)                  :: natom_of_kind

      INTEGER                                            :: ikind, n, nkind

      IF (ASSOCIATED(qs_force)) CALL deallocate_qs_force(qs_force)

      nkind = SIZE(natom_of_kind)

      ALLOCATE (qs_force(nkind))

      DO ikind = 1, nkind
         n = natom_of_kind(ikind)
         ALLOCATE (qs_force(ikind)%all_potential(3, n))
         ALLOCATE (qs_force(ikind)%core_overlap(3, n))
         ALLOCATE (qs_force(ikind)%gth_ppl(3, n))
         ALLOCATE (qs_force(ikind)%gth_nlcc(3, n))
         ALLOCATE (qs_force(ikind)%gth_ppnl(3, n))
         ALLOCATE (qs_force(ikind)%kinetic(3, n))
         ALLOCATE (qs_force(ikind)%overlap(3, n))
         ALLOCATE (qs_force(ikind)%overlap_admm(3, n))
         ALLOCATE (qs_force(ikind)%rho_core(3, n))
         ALLOCATE (qs_force(ikind)%rho_elec(3, n))
         ALLOCATE (qs_force(ikind)%rho_lri_elec(3, n))
         ALLOCATE (qs_force(ikind)%vhxc_atom(3, n))
         ALLOCATE (qs_force(ikind)%g0s_Vh_elec(3, n))
         ALLOCATE (qs_force(ikind)%repulsive(3, n))
         ALLOCATE (qs_force(ikind)%dispersion(3, n))
         ALLOCATE (qs_force(ikind)%gcp(3, n))
         ALLOCATE (qs_force(ikind)%other(3, n))
         ALLOCATE (qs_force(ikind)%ch_pulay(3, n))
         ALLOCATE (qs_force(ikind)%ehrenfest(3, n))
         ALLOCATE (qs_force(ikind)%efield(3, n))
         ALLOCATE (qs_force(ikind)%eev(3, n))
         ! Always initialize ch_pulay to zero..
         qs_force(ikind)%ch_pulay = 0.0_dp
         ALLOCATE (qs_force(ikind)%fock_4c(3, n))
         ALLOCATE (qs_force(ikind)%mp2_non_sep(3, n))
         ALLOCATE (qs_force(ikind)%total(3, n))
      END DO

   END SUBROUTINE allocate_qs_force

! **************************************************************************************************
!> \brief   Deallocate a Quickstep force data structure.
!> \param qs_force ...
!> \date    05.06.2002
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE deallocate_qs_force(qs_force)

      TYPE(qs_force_type), DIMENSION(:), POINTER         :: qs_force

      INTEGER                                            :: ikind, nkind

      CPASSERT(ASSOCIATED(qs_force))

      nkind = SIZE(qs_force)

      DO ikind = 1, nkind

         IF (ASSOCIATED(qs_force(ikind)%all_potential)) THEN
            DEALLOCATE (qs_force(ikind)%all_potential)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%core_overlap)) THEN
            DEALLOCATE (qs_force(ikind)%core_overlap)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%gth_ppl)) THEN
            DEALLOCATE (qs_force(ikind)%gth_ppl)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%gth_nlcc)) THEN
            DEALLOCATE (qs_force(ikind)%gth_nlcc)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%gth_ppnl)) THEN
            DEALLOCATE (qs_force(ikind)%gth_ppnl)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%kinetic)) THEN
            DEALLOCATE (qs_force(ikind)%kinetic)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%overlap)) THEN
            DEALLOCATE (qs_force(ikind)%overlap)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%overlap_admm)) THEN
            DEALLOCATE (qs_force(ikind)%overlap_admm)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%rho_core)) THEN
            DEALLOCATE (qs_force(ikind)%rho_core)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%rho_elec)) THEN
            DEALLOCATE (qs_force(ikind)%rho_elec)
         END IF
         IF (ASSOCIATED(qs_force(ikind)%rho_lri_elec)) THEN
            DEALLOCATE (qs_force(ikind)%rho_lri_elec)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%vhxc_atom)) THEN
            DEALLOCATE (qs_force(ikind)%vhxc_atom)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%g0s_Vh_elec)) THEN
            DEALLOCATE (qs_force(ikind)%g0s_Vh_elec)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%repulsive)) THEN
            DEALLOCATE (qs_force(ikind)%repulsive)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%dispersion)) THEN
            DEALLOCATE (qs_force(ikind)%dispersion)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%gcp)) THEN
            DEALLOCATE (qs_force(ikind)%gcp)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%other)) THEN
            DEALLOCATE (qs_force(ikind)%other)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%total)) THEN
            DEALLOCATE (qs_force(ikind)%total)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%ch_pulay)) THEN
            DEALLOCATE (qs_force(ikind)%ch_pulay)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%fock_4c)) THEN
            DEALLOCATE (qs_force(ikind)%fock_4c)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%mp2_non_sep)) THEN
            DEALLOCATE (qs_force(ikind)%mp2_non_sep)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%ehrenfest)) THEN
            DEALLOCATE (qs_force(ikind)%ehrenfest)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%efield)) THEN
            DEALLOCATE (qs_force(ikind)%efield)
         END IF

         IF (ASSOCIATED(qs_force(ikind)%eev)) THEN
            DEALLOCATE (qs_force(ikind)%eev)
         END IF
      END DO

      DEALLOCATE (qs_force)

   END SUBROUTINE deallocate_qs_force

! **************************************************************************************************
!> \brief    Initialize a Quickstep force data structure.
!> \param qs_force ...
!> \date    15.07.2002
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE zero_qs_force(qs_force)

      TYPE(qs_force_type), DIMENSION(:), POINTER         :: qs_force

      INTEGER                                            :: ikind

      CPASSERT(ASSOCIATED(qs_force))

      DO ikind = 1, SIZE(qs_force)
         qs_force(ikind)%all_potential(:, :) = 0.0_dp
         qs_force(ikind)%core_overlap(:, :) = 0.0_dp
         qs_force(ikind)%gth_ppl(:, :) = 0.0_dp
         qs_force(ikind)%gth_nlcc(:, :) = 0.0_dp
         qs_force(ikind)%gth_ppnl(:, :) = 0.0_dp
         qs_force(ikind)%kinetic(:, :) = 0.0_dp
         qs_force(ikind)%overlap(:, :) = 0.0_dp
         qs_force(ikind)%overlap_admm(:, :) = 0.0_dp
         qs_force(ikind)%rho_core(:, :) = 0.0_dp
         qs_force(ikind)%rho_elec(:, :) = 0.0_dp
         qs_force(ikind)%rho_lri_elec(:, :) = 0.0_dp
         qs_force(ikind)%vhxc_atom(:, :) = 0.0_dp
         qs_force(ikind)%g0s_Vh_elec(:, :) = 0.0_dp
         qs_force(ikind)%repulsive(:, :) = 0.0_dp
         qs_force(ikind)%dispersion(:, :) = 0.0_dp
         qs_force(ikind)%gcp(:, :) = 0.0_dp
         qs_force(ikind)%other(:, :) = 0.0_dp
         qs_force(ikind)%fock_4c(:, :) = 0.0_dp
         qs_force(ikind)%ehrenfest(:, :) = 0.0_dp
         qs_force(ikind)%efield(:, :) = 0.0_dp
         qs_force(ikind)%eev(:, :) = 0.0_dp
         qs_force(ikind)%mp2_non_sep(:, :) = 0.0_dp
         qs_force(ikind)%total(:, :) = 0.0_dp
      END DO

   END SUBROUTINE zero_qs_force

! **************************************************************************************************
!> \brief    Sum up two qs_force entities qs_force_out = qs_force_out + qs_force_in
!> \param qs_force_out ...
!> \param qs_force_in ...
!> \author  JGH
! **************************************************************************************************
   SUBROUTINE sum_qs_force(qs_force_out, qs_force_in)

      TYPE(qs_force_type), DIMENSION(:), POINTER         :: qs_force_out, qs_force_in

      INTEGER                                            :: ikind

      CPASSERT(ASSOCIATED(qs_force_out))
      CPASSERT(ASSOCIATED(qs_force_in))

      DO ikind = 1, SIZE(qs_force_out)
         qs_force_out(ikind)%all_potential(:, :) = qs_force_out(ikind)%all_potential(:, :) + &
                                                   qs_force_in(ikind)%all_potential(:, :)
         qs_force_out(ikind)%core_overlap(:, :) = qs_force_out(ikind)%core_overlap(:, :) + &
                                                  qs_force_in(ikind)%core_overlap(:, :)
         qs_force_out(ikind)%gth_ppl(:, :) = qs_force_out(ikind)%gth_ppl(:, :) + &
                                             qs_force_in(ikind)%gth_ppl(:, :)
         qs_force_out(ikind)%gth_nlcc(:, :) = qs_force_out(ikind)%gth_nlcc(:, :) + &
                                              qs_force_in(ikind)%gth_nlcc(:, :)
         qs_force_out(ikind)%gth_ppnl(:, :) = qs_force_out(ikind)%gth_ppnl(:, :) + &
                                              qs_force_in(ikind)%gth_ppnl(:, :)
         qs_force_out(ikind)%kinetic(:, :) = qs_force_out(ikind)%kinetic(:, :) + &
                                             qs_force_in(ikind)%kinetic(:, :)
         qs_force_out(ikind)%overlap(:, :) = qs_force_out(ikind)%overlap(:, :) + &
                                             qs_force_in(ikind)%overlap(:, :)
         qs_force_out(ikind)%overlap_admm(:, :) = qs_force_out(ikind)%overlap_admm(:, :) + &
                                                  qs_force_in(ikind)%overlap_admm(:, :)
         qs_force_out(ikind)%rho_core(:, :) = qs_force_out(ikind)%rho_core(:, :) + &
                                              qs_force_in(ikind)%rho_core(:, :)
         qs_force_out(ikind)%rho_elec(:, :) = qs_force_out(ikind)%rho_elec(:, :) + &
                                              qs_force_in(ikind)%rho_elec(:, :)
         qs_force_out(ikind)%rho_lri_elec(:, :) = qs_force_out(ikind)%rho_lri_elec(:, :) + &
                                                  qs_force_in(ikind)%rho_lri_elec(:, :)
         qs_force_out(ikind)%vhxc_atom(:, :) = qs_force_out(ikind)%vhxc_atom(:, :) + &
                                               qs_force_in(ikind)%vhxc_atom(:, :)
         qs_force_out(ikind)%g0s_Vh_elec(:, :) = qs_force_out(ikind)%g0s_Vh_elec(:, :) + &
                                                 qs_force_in(ikind)%g0s_Vh_elec(:, :)
         qs_force_out(ikind)%repulsive(:, :) = qs_force_out(ikind)%repulsive(:, :) + &
                                               qs_force_in(ikind)%repulsive(:, :)
         qs_force_out(ikind)%dispersion(:, :) = qs_force_out(ikind)%dispersion(:, :) + &
                                                qs_force_in(ikind)%dispersion(:, :)
         qs_force_out(ikind)%gcp(:, :) = qs_force_out(ikind)%gcp(:, :) + &
                                         qs_force_in(ikind)%gcp(:, :)
         qs_force_out(ikind)%other(:, :) = qs_force_out(ikind)%other(:, :) + &
                                           qs_force_in(ikind)%other(:, :)
         qs_force_out(ikind)%fock_4c(:, :) = qs_force_out(ikind)%fock_4c(:, :) + &
                                             qs_force_in(ikind)%fock_4c(:, :)
         qs_force_out(ikind)%ehrenfest(:, :) = qs_force_out(ikind)%ehrenfest(:, :) + &
                                               qs_force_in(ikind)%ehrenfest(:, :)
         qs_force_out(ikind)%efield(:, :) = qs_force_out(ikind)%efield(:, :) + &
                                            qs_force_in(ikind)%efield(:, :)
         qs_force_out(ikind)%eev(:, :) = qs_force_out(ikind)%eev(:, :) + &
                                         qs_force_in(ikind)%eev(:, :)
         qs_force_out(ikind)%mp2_non_sep(:, :) = qs_force_out(ikind)%mp2_non_sep(:, :) + &
                                                 qs_force_in(ikind)%mp2_non_sep(:, :)
         qs_force_out(ikind)%total(:, :) = qs_force_out(ikind)%total(:, :) + &
                                           qs_force_in(ikind)%total(:, :)
      END DO

   END SUBROUTINE sum_qs_force

! **************************************************************************************************
!> \brief    Replicate and sum up the force
!> \param qs_force ...
!> \param para_env ...
!> \date    25.05.2016
!> \author  JHU
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE replicate_qs_force(qs_force, para_env)

      TYPE(qs_force_type), DIMENSION(:), POINTER         :: qs_force
      TYPE(cp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: ikind

      !  *** replicate forces ***
      DO ikind = 1, SIZE(qs_force)
         CALL mp_sum(qs_force(ikind)%overlap, para_env%group)
         CALL mp_sum(qs_force(ikind)%overlap_admm, para_env%group)
         CALL mp_sum(qs_force(ikind)%kinetic, para_env%group)
         CALL mp_sum(qs_force(ikind)%gth_ppl, para_env%group)
         CALL mp_sum(qs_force(ikind)%gth_nlcc, para_env%group)
         CALL mp_sum(qs_force(ikind)%gth_ppnl, para_env%group)
         CALL mp_sum(qs_force(ikind)%all_potential, para_env%group)
         CALL mp_sum(qs_force(ikind)%core_overlap, para_env%group)
         CALL mp_sum(qs_force(ikind)%rho_core, para_env%group)
         CALL mp_sum(qs_force(ikind)%rho_elec, para_env%group)
         CALL mp_sum(qs_force(ikind)%rho_lri_elec, para_env%group)
         CALL mp_sum(qs_force(ikind)%vhxc_atom, para_env%group)
         CALL mp_sum(qs_force(ikind)%g0s_Vh_elec, para_env%group)
         CALL mp_sum(qs_force(ikind)%fock_4c, para_env%group)
         CALL mp_sum(qs_force(ikind)%mp2_non_sep, para_env%group)
         CALL mp_sum(qs_force(ikind)%repulsive, para_env%group)
         CALL mp_sum(qs_force(ikind)%dispersion, para_env%group)
         CALL mp_sum(qs_force(ikind)%gcp, para_env%group)
         CALL mp_sum(qs_force(ikind)%ehrenfest, para_env%group)

         qs_force(ikind)%total(:, :) = qs_force(ikind)%total(:, :) + &
                                       qs_force(ikind)%core_overlap(:, :) + &
                                       qs_force(ikind)%gth_ppl(:, :) + &
                                       qs_force(ikind)%gth_nlcc(:, :) + &
                                       qs_force(ikind)%gth_ppnl(:, :) + &
                                       qs_force(ikind)%all_potential(:, :) + &
                                       qs_force(ikind)%kinetic(:, :) + &
                                       qs_force(ikind)%overlap(:, :) + &
                                       qs_force(ikind)%overlap_admm(:, :) + &
                                       qs_force(ikind)%rho_core(:, :) + &
                                       qs_force(ikind)%rho_elec(:, :) + &
                                       qs_force(ikind)%rho_lri_elec(:, :) + &
                                       qs_force(ikind)%vhxc_atom(:, :) + &
                                       qs_force(ikind)%g0s_Vh_elec(:, :) + &
                                       qs_force(ikind)%fock_4c(:, :) + &
                                       qs_force(ikind)%mp2_non_sep(:, :) + &
                                       qs_force(ikind)%repulsive(:, :) + &
                                       qs_force(ikind)%dispersion(:, :) + &
                                       qs_force(ikind)%gcp(:, :) + &
                                       qs_force(ikind)%ehrenfest(:, :) + &
                                       qs_force(ikind)%efield(:, :) + &
                                       qs_force(ikind)%eev(:, :)
      END DO

   END SUBROUTINE replicate_qs_force

! **************************************************************************************************
!> \brief Add force to a force_type  variable.
!> \param force Input force, dimension (3,natom)
!> \param qs_force The force type variable to be used
!> \param forcetype ...
!> \param atomic_kind_set ...
!> \par History
!>      07.2014 JGH
!> \author JGH
! **************************************************************************************************
   SUBROUTINE add_qs_force(force, qs_force, forcetype, atomic_kind_set)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: force
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: qs_force
      CHARACTER(LEN=*), INTENT(IN)                       :: forcetype
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      INTEGER                                            :: ia, iatom, ikind, natom_kind
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind

!   ------------------------------------------------------------------------

      CPASSERT(ASSOCIATED(qs_force))

      SELECT CASE (forcetype)
      CASE ("overlap_admm")
         DO ikind = 1, SIZE(atomic_kind_set, 1)
            atomic_kind => atomic_kind_set(ikind)
            CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom_kind)
            DO ia = 1, natom_kind
               iatom = atomic_kind%atom_list(ia)
               qs_force(ikind)%overlap_admm(:, ia) = qs_force(ikind)%overlap_admm(:, ia) + force(:, iatom)
            END DO
         END DO
      CASE DEFAULT
         CPABORT("")
      END SELECT

   END SUBROUTINE add_qs_force

! **************************************************************************************************
!> \brief Put force to a force_type  variable.
!> \param force Input force, dimension (3,natom)
!> \param qs_force The force type variable to be used
!> \param forcetype ...
!> \param atomic_kind_set ...
!> \par History
!>      09.2019 JGH
!> \author JGH
! **************************************************************************************************
   SUBROUTINE put_qs_force(force, qs_force, forcetype, atomic_kind_set)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: force
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: qs_force
      CHARACTER(LEN=*), INTENT(IN)                       :: forcetype
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      INTEGER                                            :: ia, iatom, ikind, natom_kind
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind

!   ------------------------------------------------------------------------

      SELECT CASE (forcetype)
      CASE ("dispersion")
         DO ikind = 1, SIZE(atomic_kind_set, 1)
            atomic_kind => atomic_kind_set(ikind)
            CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom_kind)
            DO ia = 1, natom_kind
               iatom = atomic_kind%atom_list(ia)
               qs_force(ikind)%dispersion(:, ia) = force(:, iatom)
            END DO
         END DO
      CASE DEFAULT
         CPABORT("")
      END SELECT

   END SUBROUTINE put_qs_force

! **************************************************************************************************
!> \brief Get force from a force_type  variable.
!> \param force Input force, dimension (3,natom)
!> \param qs_force The force type variable to be used
!> \param forcetype ...
!> \param atomic_kind_set ...
!> \par History
!>      09.2019 JGH
!> \author JGH
! **************************************************************************************************
   SUBROUTINE get_qs_force(force, qs_force, forcetype, atomic_kind_set)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: force
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: qs_force
      CHARACTER(LEN=*), INTENT(IN)                       :: forcetype
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      INTEGER                                            :: ia, iatom, ikind, natom_kind
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind

!   ------------------------------------------------------------------------

      SELECT CASE (forcetype)
      CASE ("dispersion")
         DO ikind = 1, SIZE(atomic_kind_set, 1)
            atomic_kind => atomic_kind_set(ikind)
            CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom_kind)
            DO ia = 1, natom_kind
               iatom = atomic_kind%atom_list(ia)
               force(:, iatom) = qs_force(ikind)%dispersion(:, ia)
            END DO
         END DO
      CASE DEFAULT
         CPABORT("")
      END SELECT

   END SUBROUTINE get_qs_force

! **************************************************************************************************
!> \brief Get current total force
!> \param force Input force, dimension (3,natom)
!> \param qs_force The force type variable to be used
!> \param atomic_kind_set ...
!> \par History
!>      09.2019 JGH
!> \author JGH
! **************************************************************************************************
   SUBROUTINE total_qs_force(force, qs_force, atomic_kind_set)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: force
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: qs_force
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

      INTEGER                                            :: ia, iatom, ikind, natom_kind
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind

!   ------------------------------------------------------------------------

      force(:, :) = 0.0_dp
      DO ikind = 1, SIZE(atomic_kind_set, 1)
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom_kind)
         DO ia = 1, natom_kind
            iatom = atomic_kind%atom_list(ia)
            force(:, iatom) = qs_force(ikind)%core_overlap(:, ia) + &
                              qs_force(ikind)%gth_ppl(:, ia) + &
                              qs_force(ikind)%gth_nlcc(:, ia) + &
                              qs_force(ikind)%gth_ppnl(:, ia) + &
                              qs_force(ikind)%all_potential(:, ia) + &
                              qs_force(ikind)%kinetic(:, ia) + &
                              qs_force(ikind)%overlap(:, ia) + &
                              qs_force(ikind)%overlap_admm(:, ia) + &
                              qs_force(ikind)%rho_core(:, ia) + &
                              qs_force(ikind)%rho_elec(:, ia) + &
                              qs_force(ikind)%rho_lri_elec(:, ia) + &
                              qs_force(ikind)%vhxc_atom(:, ia) + &
                              qs_force(ikind)%g0s_Vh_elec(:, ia) + &
                              qs_force(ikind)%fock_4c(:, ia) + &
                              qs_force(ikind)%mp2_non_sep(:, ia) + &
                              qs_force(ikind)%repulsive(:, ia) + &
                              qs_force(ikind)%dispersion(:, ia) + &
                              qs_force(ikind)%gcp(:, ia) + &
                              qs_force(ikind)%ehrenfest(:, ia) + &
                              qs_force(ikind)%efield(:, ia) + &
                              qs_force(ikind)%eev(:, ia)
         END DO
      END DO

   END SUBROUTINE total_qs_force

END MODULE qs_force_types
