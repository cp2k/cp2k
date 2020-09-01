!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      JGH [04042007] refactoring
! **************************************************************************************************
MODULE virial_types

   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'virial_types'
   INTEGER, PRIVATE, SAVE :: last_virial_id_nr = 0

   PUBLIC :: virial_type, virial_p_type, virial_set, cp_virial, &
             virial_create, virial_release, zero_virial, sym_virial

! **************************************************************************************************
   TYPE virial_type
      INTEGER                            :: ref_count, id_nr
      REAL(KIND=dp), DIMENSION(3, 3) :: pv_total, &
                                        pv_virial, &
                                        pv_kinetic, &
                                        pv_xc, &
                                        pv_constraint
      REAL(KIND=dp), DIMENSION(3, 3) :: pv_ekin, &
                                        pv_overlap, &
                                        pv_hartree, &
                                        pv_ppl, &
                                        pv_ppnl, &
                                        pv_exc, &
                                        pv_vdw, &
                                        pv_fock_4c
      LOGICAL                        :: pv_availability, &
                                        pv_calculate, &
                                        pv_numer, &
                                        pv_diagonal
   END TYPE virial_type

! **************************************************************************************************
   TYPE virial_p_type
      TYPE(virial_type), POINTER        :: virial
   END TYPE virial_p_type

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
      virial_out%pv_constraint = virial_in%pv_constraint
      virial_out%pv_availability = virial_in%pv_availability
      virial_out%pv_calculate = virial_in%pv_calculate
      virial_out%pv_numer = virial_in%pv_numer
      virial_out%pv_diagonal = virial_in%pv_diagonal

      virial_out%pv_ekin = virial_in%pv_ekin
      virial_out%pv_hartree = virial_in%pv_hartree
      virial_out%pv_overlap = virial_in%pv_overlap
      virial_out%pv_ppl = virial_in%pv_ppl
      virial_out%pv_ppnl = virial_in%pv_ppnl
      virial_out%pv_exc = virial_in%pv_exc
      virial_out%pv_vdw = virial_in%pv_vdw
      virial_out%pv_fock_4c = virial_in%pv_fock_4c

   END SUBROUTINE cp_virial

! **************************************************************************************************
!> \brief   symmetrize the virial components
!> \param virial ...
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE sym_virial(virial)
      TYPE(virial_type), POINTER                         :: virial

      INTEGER                                            :: i, j

      DO i = 1, 3
         DO j = 1, i - 1
            ! Symmetrize total
            virial%pv_total(j, i) = (virial%pv_total(i, j) + virial%pv_total(j, i))*0.5_dp
            virial%pv_total(i, j) = virial%pv_total(j, i)
            ! Symmetrize Kinetic
            virial%pv_kinetic(j, i) = (virial%pv_kinetic(i, j) + virial%pv_kinetic(j, i))*0.5_dp
            virial%pv_kinetic(i, j) = virial%pv_kinetic(j, i)
            ! Symmetrize Virial
            virial%pv_virial(j, i) = (virial%pv_virial(i, j) + virial%pv_virial(j, i))*0.5_dp
            virial%pv_virial(i, j) = virial%pv_virial(j, i)
            ! Symmetrize XC
            virial%pv_xc(j, i) = (virial%pv_xc(i, j) + virial%pv_xc(j, i))*0.5_dp
            virial%pv_xc(i, j) = virial%pv_xc(j, i)

            ! Symmetrize tensor parts
            virial%pv_ekin(j, i) = (virial%pv_ekin(i, j) + virial%pv_ekin(j, i))*0.5_dp
            virial%pv_ekin(i, j) = virial%pv_ekin(j, i)
            !
            virial%pv_overlap(j, i) = (virial%pv_overlap(i, j) + virial%pv_overlap(j, i))*0.5_dp
            virial%pv_overlap(i, j) = virial%pv_overlap(j, i)
            !
            virial%pv_hartree(j, i) = (virial%pv_hartree(i, j) + virial%pv_hartree(j, i))*0.5_dp
            virial%pv_hartree(i, j) = virial%pv_hartree(j, i)
            !
            virial%pv_exc(j, i) = (virial%pv_exc(i, j) + virial%pv_exc(j, i))*0.5_dp
            virial%pv_exc(i, j) = virial%pv_exc(j, i)
            !
            virial%pv_vdw(j, i) = (virial%pv_vdw(i, j) + virial%pv_vdw(j, i))*0.5_dp
            virial%pv_vdw(i, j) = virial%pv_vdw(j, i)
            !
            virial%pv_ppl(j, i) = (virial%pv_ppl(i, j) + virial%pv_ppl(j, i))*0.5_dp
            virial%pv_ppl(i, j) = virial%pv_ppl(j, i)
            !
            virial%pv_ppnl(j, i) = (virial%pv_ppnl(i, j) + virial%pv_ppnl(j, i))*0.5_dp
            virial%pv_ppnl(i, j) = virial%pv_ppnl(j, i)
            !
            virial%pv_fock_4c(j, i) = (virial%pv_fock_4c(i, j) + virial%pv_fock_4c(j, i))*0.5_dp
            virial%pv_fock_4c(i, j) = virial%pv_fock_4c(j, i)

            ! Symmetrize constraints
            virial%pv_constraint(j, i) = (virial%pv_constraint(i, j) + virial%pv_constraint(j, i))*0.5_dp
            virial%pv_constraint(i, j) = virial%pv_constraint(j, i)
         END DO
      END DO

   END SUBROUTINE sym_virial

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
      virial%pv_constraint = 0.0_dp
      virial%pv_ekin = 0.0_dp
      virial%pv_overlap = 0.0_dp
      virial%pv_hartree = 0.0_dp
      virial%pv_exc = 0.0_dp
      virial%pv_vdw = 0.0_dp
      virial%pv_ppl = 0.0_dp
      virial%pv_ppnl = 0.0_dp
      virial%pv_fock_4c = 0.0_dp
      IF (my_reset) THEN
         virial%pv_availability = .FALSE.
         virial%pv_numer = .FALSE.
         virial%pv_calculate = .FALSE.
         virial%pv_diagonal = .FALSE.
      END IF
   END SUBROUTINE zero_virial

! **************************************************************************************************
!> \brief ...
!> \param virial ...
!> \param pv_total ...
!> \param pv_virial ...
!> \param pv_xc ...
!> \param pv_constraint ...
!> \param pv_kinetic ...
!> \param pv_ekin ...
!> \param pv_overlap ...
!> \param pv_hartree ...
!> \param pv_exc ...
!> \param pv_vdw ...
!> \param pv_ppl ...
!> \param pv_ppnl ...
!> \param pv_fock_4c ...
!> \param pv_availability ...
!> \param pv_calculate ...
!> \param pv_numer ...
!> \param pv_diagonal ...
! **************************************************************************************************
   SUBROUTINE virial_set(virial, pv_total, pv_virial, pv_xc, pv_constraint, pv_kinetic, &
                         pv_ekin, pv_overlap, pv_hartree, pv_exc, pv_vdw, pv_ppl, pv_ppnl, pv_fock_4c, &
                         pv_availability, pv_calculate, pv_numer, pv_diagonal)

      TYPE(virial_type), POINTER                         :: virial
      REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL :: pv_total, pv_virial, pv_xc, pv_constraint, &
         pv_kinetic, pv_ekin, pv_overlap, pv_hartree, pv_exc, pv_vdw, pv_ppl, pv_ppnl, pv_fock_4c
      LOGICAL, OPTIONAL                                  :: pv_availability, pv_calculate, pv_numer, &
                                                            pv_diagonal

      IF (PRESENT(pv_total)) virial%pv_total = pv_total
      IF (PRESENT(pv_virial)) virial%pv_virial = pv_virial
      IF (PRESENT(pv_xc)) virial%pv_xc = pv_xc
      IF (PRESENT(pv_constraint)) virial%pv_constraint = pv_constraint
      IF (PRESENT(pv_kinetic)) virial%pv_kinetic = pv_kinetic
      IF (PRESENT(pv_ekin)) virial%pv_ekin = pv_ekin
      IF (PRESENT(pv_overlap)) virial%pv_overlap = pv_overlap
      IF (PRESENT(pv_hartree)) virial%pv_hartree = pv_hartree
      IF (PRESENT(pv_exc)) virial%pv_exc = pv_exc
      IF (PRESENT(pv_vdw)) virial%pv_vdw = pv_vdw
      IF (PRESENT(pv_ppl)) virial%pv_ppl = pv_ppl
      IF (PRESENT(pv_ppnl)) virial%pv_ppnl = pv_ppnl
      IF (PRESENT(pv_fock_4c)) virial%pv_fock_4c = pv_fock_4c
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
!>      04.2003 created [fawzi]
!> \author fawzi
!> \note
!>      see doc/ReferenceCounting.html
! **************************************************************************************************
   SUBROUTINE virial_release(virial)
      TYPE(virial_type), POINTER                         :: virial

      IF (ASSOCIATED(virial)) THEN
         CPASSERT(virial%ref_count > 0)
         virial%ref_count = virial%ref_count - 1
         IF (virial%ref_count .EQ. 0) THEN
            DEALLOCATE (virial)
         ENDIF
         NULLIFY (virial)
      END IF
   END SUBROUTINE virial_release

END MODULE virial_types

