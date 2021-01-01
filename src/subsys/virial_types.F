!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      - Refactoring (4.4.2007, JGH)
!>      - Revise virial components (16.10.2020, MK)
! **************************************************************************************************
MODULE virial_types

   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'virial_types'

   INTEGER, PRIVATE, SAVE :: last_virial_id_nr = 0

   PUBLIC :: virial_type, virial_p_type

   TYPE virial_type
      INTEGER                        :: ref_count, id_nr
      REAL(KIND=dp), DIMENSION(3, 3) :: pv_total, &
                                        pv_kinetic, &
                                        pv_virial, &
                                        pv_xc, &
                                        pv_fock_4c, &
                                        pv_constraint
      REAL(KIND=dp), DIMENSION(3, 3) :: pv_overlap, &
                                        pv_ekinetic, &
                                        pv_ppl, &
                                        pv_ppnl, &
                                        pv_ecore_overlap, &
                                        pv_ehartree, &
                                        pv_exc, &
                                        pv_exx, &
                                        pv_vdw, &
                                        pv_mp2, &
                                        pv_nlcc, &
                                        pv_gapw, &
                                        pv_lrigpw
      LOGICAL                        :: pv_availability, &
                                        pv_calculate, &
                                        pv_numer, &
                                        pv_diagonal
   END TYPE virial_type

   TYPE virial_p_type
      TYPE(virial_type), POINTER     :: virial
   END TYPE virial_p_type

   PUBLIC :: cp_virial, virial_create, virial_release, virial_set, &
             symmetrize_virial, zero_virial

CONTAINS

! **************************************************************************************************
!> \brief   copy virial_in into virial_out
!> \param virial_in ...
!> \param virial_out ...
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE cp_virial(virial_in, virial_out)
      TYPE(virial_type), POINTER                         :: virial_in, virial_out

      virial_out%pv_total = virial_in%pv_total
      virial_out%pv_kinetic = virial_in%pv_kinetic
      virial_out%pv_virial = virial_in%pv_virial
      virial_out%pv_xc = virial_in%pv_xc
      virial_out%pv_fock_4c = virial_in%pv_fock_4c
      virial_out%pv_constraint = virial_in%pv_constraint

      virial_out%pv_overlap = virial_in%pv_overlap
      virial_out%pv_ekinetic = virial_in%pv_ekinetic
      virial_out%pv_ppl = virial_in%pv_ppl
      virial_out%pv_ppnl = virial_in%pv_ppnl
      virial_out%pv_ecore_overlap = virial_in%pv_ecore_overlap
      virial_out%pv_ehartree = virial_in%pv_ehartree
      virial_out%pv_exc = virial_in%pv_exc
      virial_out%pv_exx = virial_in%pv_exx
      virial_out%pv_vdw = virial_in%pv_vdw
      virial_out%pv_mp2 = virial_in%pv_mp2
      virial_out%pv_nlcc = virial_in%pv_nlcc
      virial_out%pv_gapw = virial_in%pv_gapw
      virial_out%pv_lrigpw = virial_in%pv_lrigpw

      virial_out%pv_availability = virial_in%pv_availability
      virial_out%pv_calculate = virial_in%pv_calculate
      virial_out%pv_numer = virial_in%pv_numer
      virial_out%pv_diagonal = virial_in%pv_diagonal

   END SUBROUTINE cp_virial

! **************************************************************************************************
!> \brief   Symmetrize the virial components
!> \param virial ...
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE symmetrize_virial(virial)
      TYPE(virial_type), POINTER                         :: virial

      INTEGER                                            :: i, j

      DO i = 1, 3
         DO j = 1, i - 1
            virial%pv_total(j, i) = 0.5_dp*(virial%pv_total(i, j) + virial%pv_total(j, i))
            virial%pv_total(i, j) = virial%pv_total(j, i)
            virial%pv_kinetic(j, i) = 0.5_dp*(virial%pv_kinetic(i, j) + virial%pv_kinetic(j, i))
            virial%pv_kinetic(i, j) = virial%pv_kinetic(j, i)
            virial%pv_virial(j, i) = 0.5_dp*(virial%pv_virial(i, j) + virial%pv_virial(j, i))
            virial%pv_virial(i, j) = virial%pv_virial(j, i)
            virial%pv_xc(j, i) = 0.5_dp*(virial%pv_xc(i, j) + virial%pv_xc(j, i))
            virial%pv_xc(i, j) = virial%pv_xc(j, i)
            virial%pv_fock_4c(j, i) = 0.5_dp*(virial%pv_fock_4c(i, j) + virial%pv_fock_4c(j, i))
            virial%pv_fock_4c(i, j) = virial%pv_fock_4c(j, i)
            virial%pv_constraint(j, i) = 0.5_dp*(virial%pv_constraint(i, j) + virial%pv_constraint(j, i))
            virial%pv_constraint(i, j) = virial%pv_constraint(j, i)
            ! Virial components
            virial%pv_overlap(j, i) = 0.5_dp*(virial%pv_overlap(i, j) + virial%pv_overlap(j, i))
            virial%pv_overlap(i, j) = virial%pv_overlap(j, i)
            virial%pv_ekinetic(j, i) = 0.5_dp*(virial%pv_ekinetic(i, j) + virial%pv_ekinetic(j, i))
            virial%pv_ekinetic(i, j) = virial%pv_ekinetic(j, i)
            virial%pv_ppl(j, i) = 0.5_dp*(virial%pv_ppl(i, j) + virial%pv_ppl(j, i))
            virial%pv_ppl(i, j) = virial%pv_ppl(j, i)
            virial%pv_ppnl(j, i) = 0.5_dp*(virial%pv_ppnl(i, j) + virial%pv_ppnl(j, i))
            virial%pv_ppnl(i, j) = virial%pv_ppnl(j, i)
            virial%pv_ecore_overlap(j, i) = 0.5_dp*(virial%pv_ecore_overlap(i, j) + virial%pv_ecore_overlap(j, i))
            virial%pv_ecore_overlap(i, j) = virial%pv_ecore_overlap(j, i)
            virial%pv_ehartree(j, i) = 0.5_dp*(virial%pv_ehartree(i, j) + virial%pv_ehartree(j, i))
            virial%pv_ehartree(i, j) = virial%pv_ehartree(j, i)
            virial%pv_exc(j, i) = 0.5_dp*(virial%pv_exc(i, j) + virial%pv_exc(j, i))
            virial%pv_exc(i, j) = virial%pv_exc(j, i)
            virial%pv_exx(j, i) = 0.5_dp*(virial%pv_exx(i, j) + virial%pv_exx(j, i))
            virial%pv_exx(i, j) = virial%pv_exx(j, i)
            virial%pv_vdw(j, i) = 0.5_dp*(virial%pv_vdw(i, j) + virial%pv_vdw(j, i))
            virial%pv_vdw(i, j) = virial%pv_vdw(j, i)
            virial%pv_mp2(j, i) = 0.5_dp*(virial%pv_mp2(i, j) + virial%pv_mp2(j, i))
            virial%pv_mp2(i, j) = virial%pv_mp2(j, i)
            virial%pv_nlcc(j, i) = 0.5_dp*(virial%pv_nlcc(i, j) + virial%pv_nlcc(j, i))
            virial%pv_nlcc(i, j) = virial%pv_nlcc(j, i)
            virial%pv_gapw(j, i) = 0.5_dp*(virial%pv_gapw(i, j) + virial%pv_gapw(j, i))
            virial%pv_gapw(i, j) = virial%pv_gapw(j, i)
            virial%pv_lrigpw(j, i) = 0.5_dp*(virial%pv_lrigpw(i, j) + virial%pv_lrigpw(j, i))
            virial%pv_lrigpw(i, j) = virial%pv_lrigpw(j, i)
         END DO
      END DO

   END SUBROUTINE symmetrize_virial

! **************************************************************************************************
!> \brief ...
!> \param virial ...
!> \param reset ...
! **************************************************************************************************
   SUBROUTINE zero_virial(virial, reset)
      TYPE(virial_type), POINTER                         :: virial
      LOGICAL, INTENT(IN), OPTIONAL                      :: reset

      LOGICAL                                            :: my_reset

      my_reset = .TRUE.
      IF (PRESENT(reset)) my_reset = reset

      virial%pv_total = 0.0_dp
      virial%pv_kinetic = 0.0_dp
      virial%pv_virial = 0.0_dp
      virial%pv_xc = 0.0_dp
      virial%pv_fock_4c = 0.0_dp
      virial%pv_constraint = 0.0_dp

      virial%pv_overlap = 0.0_dp
      virial%pv_ekinetic = 0.0_dp
      virial%pv_ppl = 0.0_dp
      virial%pv_ppnl = 0.0_dp
      virial%pv_ecore_overlap = 0.0_dp
      virial%pv_ehartree = 0.0_dp
      virial%pv_exc = 0.0_dp
      virial%pv_exx = 0.0_dp
      virial%pv_vdw = 0.0_dp
      virial%pv_mp2 = 0.0_dp
      virial%pv_nlcc = 0.0_dp
      virial%pv_gapw = 0.0_dp
      virial%pv_lrigpw = 0.0_dp

      IF (my_reset) THEN
         virial%pv_availability = .FALSE.
         virial%pv_calculate = .FALSE.
         virial%pv_numer = .FALSE.
         virial%pv_diagonal = .FALSE.
      END IF

   END SUBROUTINE zero_virial

! **************************************************************************************************
!> \brief ...
!> \param virial ...
!> \param pv_total ...
!> \param pv_kinetic ...
!> \param pv_virial ...
!> \param pv_xc ...
!> \param pv_fock_4c ...
!> \param pv_constraint ...
!> \param pv_overlap ...
!> \param pv_ekinetic ...
!> \param pv_ppl ...
!> \param pv_ppnl ...
!> \param pv_ecore_overlap ...
!> \param pv_ehartree ...
!> \param pv_exc ...
!> \param pv_exx ...
!> \param pv_vdw ...
!> \param pv_mp2 ...
!> \param pv_nlcc ...
!> \param pv_gapw ...
!> \param pv_lrigpw ...
!> \param pv_availability ...
!> \param pv_calculate ...
!> \param pv_numer ...
!> \param pv_diagonal ...
! **************************************************************************************************
   SUBROUTINE virial_set(virial, pv_total, pv_kinetic, pv_virial, pv_xc, pv_fock_4c, pv_constraint, &
                         pv_overlap, pv_ekinetic, pv_ppl, pv_ppnl, pv_ecore_overlap, pv_ehartree, &
                         pv_exc, pv_exx, pv_vdw, pv_mp2, pv_nlcc, pv_gapw, pv_lrigpw, &
                         pv_availability, pv_calculate, pv_numer, pv_diagonal)

      TYPE(virial_type), POINTER                         :: virial
      REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL :: pv_total, pv_kinetic, pv_virial, pv_xc, &
         pv_fock_4c, pv_constraint, pv_overlap, pv_ekinetic, pv_ppl, pv_ppnl, pv_ecore_overlap, &
         pv_ehartree, pv_exc, pv_exx, pv_vdw, pv_mp2, pv_nlcc, pv_gapw, pv_lrigpw
      LOGICAL, OPTIONAL                                  :: pv_availability, pv_calculate, pv_numer, &
                                                            pv_diagonal

      IF (PRESENT(pv_total)) virial%pv_total = pv_total
      IF (PRESENT(pv_kinetic)) virial%pv_kinetic = pv_kinetic
      IF (PRESENT(pv_virial)) virial%pv_virial = pv_virial
      IF (PRESENT(pv_xc)) virial%pv_xc = pv_xc
      IF (PRESENT(pv_fock_4c)) virial%pv_fock_4c = pv_fock_4c
      IF (PRESENT(pv_constraint)) virial%pv_constraint = pv_constraint

      IF (PRESENT(pv_overlap)) virial%pv_overlap = pv_overlap
      IF (PRESENT(pv_ekinetic)) virial%pv_ekinetic = pv_ekinetic
      IF (PRESENT(pv_ppl)) virial%pv_ppl = pv_ppl
      IF (PRESENT(pv_ppnl)) virial%pv_ppnl = pv_ppnl
      IF (PRESENT(pv_ecore_overlap)) virial%pv_ecore_overlap = pv_ecore_overlap
      IF (PRESENT(pv_ehartree)) virial%pv_ehartree = pv_ehartree
      IF (PRESENT(pv_exc)) virial%pv_exc = pv_exc
      IF (PRESENT(pv_exx)) virial%pv_exx = pv_exx
      IF (PRESENT(pv_vdw)) virial%pv_vdw = pv_vdw
      IF (PRESENT(pv_mp2)) virial%pv_mp2 = pv_mp2
      IF (PRESENT(pv_nlcc)) virial%pv_nlcc = pv_nlcc
      IF (PRESENT(pv_gapw)) virial%pv_gapw = pv_gapw
      IF (PRESENT(pv_lrigpw)) virial%pv_lrigpw = pv_lrigpw

      IF (PRESENT(pv_availability)) virial%pv_availability = pv_availability
      IF (PRESENT(pv_calculate)) virial%pv_calculate = pv_calculate
      IF (PRESENT(pv_numer)) virial%pv_numer = pv_numer
      IF (PRESENT(pv_diagonal)) virial%pv_diagonal = pv_diagonal

   END SUBROUTINE virial_set

! **************************************************************************************************
!> \brief ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE virial_create(virial)
      TYPE(virial_type), POINTER                         :: virial

      ALLOCATE (virial)
      CALL zero_virial(virial)
      last_virial_id_nr = last_virial_id_nr + 1
      virial%id_nr = last_virial_id_nr
      virial%ref_count = 1

   END SUBROUTINE virial_create

! **************************************************************************************************
!> \brief releases the given virial_type
!> \param virial the virial_type to release
!> \par History
!>      Created (04.2003, fawzi)
!> \author fawzi
!> \note
!>      see doc/ReferenceCounting.html
! **************************************************************************************************
   SUBROUTINE virial_release(virial)
      TYPE(virial_type), POINTER                         :: virial

      IF (ASSOCIATED(virial)) THEN
         CPASSERT(virial%ref_count > 0)
         virial%ref_count = virial%ref_count - 1
         IF (virial%ref_count == 0) THEN
            DEALLOCATE (virial)
         ENDIF
         NULLIFY (virial)
      END IF

   END SUBROUTINE virial_release

END MODULE virial_types
