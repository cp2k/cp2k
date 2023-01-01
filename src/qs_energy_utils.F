!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility subroutine for qs energy calculation
!> \par History
!>      none
!> \author MK (29.10.2002)
! **************************************************************************************************
MODULE qs_energy_utils
   USE atprop_types,                    ONLY: atprop_array_add,&
                                              atprop_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_control_utils,                ONLY: read_ddapc_section
   USE dbcsr_api,                       ONLY: &
        dbcsr_get_block_p, dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, dbcsr_type
   USE et_coupling,                     ONLY: calc_et_coupling
   USE et_coupling_proj,                ONLY: calc_et_coupling_proj
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_type
   USE pw_types,                        ONLY: pw_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_core_rspace
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_linres_module,                ONLY: linres_calculation_low
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf,                          ONLY: scf
   USE qs_tddfpt2_methods,              ONLY: tddfpt
   USE tip_scan_methods,                ONLY: tip_scanning
   USE xas_methods,                     ONLY: xas
   USE xas_tdp_methods,                 ONLY: xas_tdp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy_utils'

   PUBLIC :: qs_energies_properties

CONTAINS

! **************************************************************************************************
!> \brief Refactoring of qs_energies_scf. Moves computation of properties
!>        into separate subroutine
!> \param qs_env ...
!> \param calc_forces ...
!> \par History
!>      05.2013 created [Florian Schiffmann]
! **************************************************************************************************

   SUBROUTINE qs_energies_properties(qs_env, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies_properties'

      INTEGER                                            :: handle
      LOGICAL                                            :: do_et, do_et_proj, do_tip_scan
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_type), POINTER                             :: v_hartree_rspace
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: input, proj_section, rest_b_section, &
                                                            tip_scan_section

      NULLIFY (atprop, energy, pw_env)
      CALL timeset(routineN, handle)

      ! atomic energies using Mulliken partition
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      input=input, &
                      atprop=atprop, &
                      energy=energy, &
                      v_hartree_rspace=v_hartree_rspace, &
                      pw_env=pw_env)
      IF (atprop%energy) THEN
         CALL qs_energies_mulliken(qs_env)
         IF (.NOT. dft_control%qs_control%semi_empirical .AND. &
             .NOT. dft_control%qs_control%xtb .AND. &
             .NOT. dft_control%qs_control%dftb) THEN
            ! Nuclear charge correction
            CALL integrate_v_core_rspace(v_hartree_rspace, qs_env)
            ! Kohn-Sham Functional corrections
         END IF
         CALL atprop_array_add(atprop%atener, atprop%ateb)
         CALL atprop_array_add(atprop%atener, atprop%ateself)
         CALL atprop_array_add(atprop%atener, atprop%atexc)
         CALL atprop_array_add(atprop%atener, atprop%atecoul)
         CALL atprop_array_add(atprop%atener, atprop%atevdw)
         CALL atprop_array_add(atprop%atener, atprop%ategcp)
         CALL atprop_array_add(atprop%atener, atprop%atecc)
         CALL atprop_array_add(atprop%atener, atprop%ate1c)
      END IF

      ! ET coupling - projection-operator approach
      NULLIFY (proj_section)
      proj_section => &
         section_vals_get_subs_vals(qs_env%input, "PROPERTIES%ET_COUPLING%PROJECTION")
      CALL section_vals_get(proj_section, explicit=do_et_proj)
      IF (do_et_proj) THEN
         CALL calc_et_coupling_proj(qs_env)
      END IF

      ! **********  Calculate the electron transfer coupling elements********
      do_et = .FALSE.
      do_et = dft_control%qs_control%et_coupling_calc
      IF (do_et) THEN
         qs_env%et_coupling%energy = energy%total
         qs_env%et_coupling%keep_matrix = .TRUE.
         qs_env%et_coupling%first_run = .TRUE.
         CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., just_energy=.TRUE.)
         qs_env%et_coupling%first_run = .FALSE.
         IF (dft_control%qs_control%ddapc_restraint) THEN
            rest_b_section => section_vals_get_subs_vals(input, "PROPERTIES%ET_COUPLING%DDAPC_RESTRAINT_B")
            CALL read_ddapc_section(qs_control=dft_control%qs_control, &
                                    ddapc_restraint_section=rest_b_section)
         END IF
         CALL scf(qs_env=qs_env)
         qs_env%et_coupling%keep_matrix = .TRUE.

         CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., just_energy=.TRUE.)
         CALL calc_et_coupling(qs_env)
      END IF

      !Properties
      IF (dft_control%do_xas_calculation) THEN
         CALL xas(qs_env, dft_control)
      END IF

      IF (dft_control%do_xas_tdp_calculation) THEN
         CALL xas_tdp(qs_env)
      END IF

      ! Compute Linear Response properties as post-scf
      IF (.NOT. qs_env%linres_run) THEN
         CALL linres_calculation_low(qs_env)
      END IF

      IF (dft_control%tddfpt2_control%enabled) THEN
         CALL tddfpt(qs_env, calc_forces)
      END IF

      ! tip scan
      NULLIFY (tip_scan_section)
      tip_scan_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%TIP_SCAN")
      CALL section_vals_get(tip_scan_section, explicit=do_tip_scan)
      IF (do_tip_scan) THEN
         CALL tip_scanning(qs_env, tip_scan_section)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_energies_properties

! **************************************************************************************************
!> \brief   Use a simple Mulliken-like energy decomposition
!> \param qs_env ...
!> \date    07.2011
!> \author  JHU
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE qs_energies_mulliken(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: ispin
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_h, matrix_ks, rho_ao
      TYPE(qs_rho_type), POINTER                         :: rho

      NULLIFY (atprop, matrix_h, matrix_ks, rho, rho_ao)
      CALL get_qs_env(qs_env=qs_env, matrix_ks=matrix_ks, matrix_h=matrix_h, &
                      rho=rho, atprop=atprop)
      CALL qs_rho_get(rho, rho_ao=rho_ao)

      IF (atprop%energy) THEN
         ! E = 0.5*Tr(H*P+F*P)
         atprop%atener = 0._dp
         DO ispin = 1, SIZE(rho_ao)
            CALL atom_trace(matrix_h(1)%matrix, rho_ao(ispin)%matrix, &
                            0.5_dp, atprop%atener)
            CALL atom_trace(matrix_ks(ispin)%matrix, rho_ao(ispin)%matrix, &
                            0.5_dp, atprop%atener)
         END DO
      END IF

   END SUBROUTINE qs_energies_mulliken

! **************************************************************************************************
!> \brief Compute partial trace of product of two matrices
!> \param amat ...
!> \param bmat ...
!> \param factor ...
!> \param atrace ...
!> \par History
!>      06.2004 created [Joost VandeVondele]
!> \note
!>      charges are computed per spin in the LSD case
! **************************************************************************************************
   SUBROUTINE atom_trace(amat, bmat, factor, atrace)
      TYPE(dbcsr_type), POINTER                          :: amat, bmat
      REAL(kind=dp), INTENT(IN)                          :: factor
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atrace

      INTEGER                                            :: blk, iblock_col, iblock_row, nblock
      LOGICAL                                            :: found
      REAL(kind=dp)                                      :: btr, mult
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: a_block, b_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_get_info(bmat, nblkrows_total=nblock)
      CPASSERT(nblock == SIZE(atrace))

      CALL dbcsr_iterator_start(iter, bmat)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, b_block, blk)
         CALL dbcsr_get_block_p(matrix=amat, &
                                row=iblock_row, col=iblock_col, BLOCK=a_block, found=found)

         ! we can cycle if a block is not present
         IF (.NOT. (ASSOCIATED(b_block) .AND. ASSOCIATED(a_block))) CYCLE

         IF (iblock_row .EQ. iblock_col) THEN
            mult = 0.5_dp ! avoid double counting of diagonal blocks
         ELSE
            mult = 1.0_dp
         END IF
         btr = factor*mult*SUM(a_block*b_block)
         atrace(iblock_row) = atrace(iblock_row) + btr
         atrace(iblock_col) = atrace(iblock_col) + btr

      END DO
      CALL dbcsr_iterator_stop(iter)

   END SUBROUTINE atom_trace

END MODULE qs_energy_utils
