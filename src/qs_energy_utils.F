!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
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
   USE cp_control_utils,                ONLY: read_becke_section,&
                                              read_ddapc_section
   USE cp_external_control,             ONLY: external_control
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_get_block_p, dbcsr_get_info, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_set, dbcsr_type
   USE et_coupling,                     ONLY: calc_et_coupling
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_density_matrices,&
                                              kpoint_density_transform
   USE kpoint_types,                    ONLY: kpoint_type
   USE mp2,                             ONLY: mp2_main
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: pw_p_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integrate_potential,          ONLY: integrate_v_core_rspace
   USE qs_ks_methods,                   ONLY: calculate_w_matrix,&
                                              calculate_w_matrix_ot,&
                                              qs_ks_update_qs_env
   USE qs_linres_module,                ONLY: linres_calculation_low
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf,                          ONLY: scf
   USE qs_tddfpt2_methods,              ONLY: tddfpt
   USE scf_control_types,               ONLY: scf_control_type
   USE xas_methods,                     ONLY: xas
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy_utils'

   PUBLIC :: qs_energies_compute_matrix_w, &
             qs_energies_properties, &
             qs_energies_mp2

CONTAINS

! **************************************************************************************************
!> \brief Refactoring of qs_energies_scf. Moves computation of matrix_w
!>        into separate subroutine
!> \param qs_env ...
!> \param calc_forces ...
!> \par History
!>      05.2013 created [Florian Schiffmann]
! **************************************************************************************************

   SUBROUTINE qs_energies_compute_matrix_w(qs_env, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies_compute_matrix_w', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, is, ispin, nao, nspin
      LOGICAL                                            :: do_kpoints, has_unit_metric
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: fmwork
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_fmstruct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s, matrix_w, &
                                                            matrix_w_mp2, mo_derivs, rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, matrix_w_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      ! if calculate forces, time to compute the w matrix
      CALL get_qs_env(qs_env, has_unit_metric=has_unit_metric)

      IF (calc_forces .AND. .NOT. has_unit_metric) THEN
         CALL get_qs_env(qs_env, do_kpoints=do_kpoints)

         IF (do_kpoints) THEN

            CALL get_qs_env(qs_env, &
                            matrix_w_kp=matrix_w_kp, &
                            matrix_s_kp=matrix_s_kp, &
                            sab_orb=sab_nl, &
                            mos=mos, &
                            kpoints=kpoints)

            CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao)
            CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct, nrow_global=nao, ncol_global=nao, &
                                     template_fmstruct=mo_coeff%matrix_struct)

            ALLOCATE (fmwork(2))
            DO is = 1, SIZE(fmwork)
               NULLIFY (fmwork(is)%matrix)
               CALL cp_fm_create(fmwork(is)%matrix, matrix_struct=ao_ao_fmstruct)
            END DO
            CALL cp_fm_struct_release(ao_ao_fmstruct)

            ! energy weighted density matrices in k-space
            CALL kpoint_density_matrices(kpoints, energy_weighted=.TRUE.)
            ! energy weighted density matrices in real space
            CALL kpoint_density_transform(kpoints, matrix_w_kp, .TRUE., &
                                          matrix_s_kp(1, 1)%matrix, sab_nl, fmwork)

            DO is = 1, SIZE(fmwork)
               CALL cp_fm_release(fmwork(is)%matrix)
            END DO
            DEALLOCATE (fmwork)

         ELSE

            NULLIFY (dft_control, rho_ao)
            CALL get_qs_env(qs_env, &
                            matrix_w=matrix_w, &
                            matrix_ks=matrix_ks, &
                            matrix_s=matrix_s, &
                            matrix_w_mp2=matrix_w_mp2, &
                            mo_derivs=mo_derivs, &
                            scf_control=scf_control, &
                            mos=mos, &
                            rho=rho, &
                            dft_control=dft_control)

            CALL qs_rho_get(rho, rho_ao=rho_ao)

            nspin = SIZE(mos)
            DO ispin = 1, nspin
               mo_set => mos(ispin)%mo_set
               IF (dft_control%roks) THEN
                  IF (scf_control%use_ot) THEN
                     IF (ispin > 1) THEN
                        ! not very elegant, indeed ...
                        CALL dbcsr_set(matrix_w(ispin)%matrix, 0.0_dp)
                     ELSE
                        CALL calculate_w_matrix_ot(mo_set, mo_derivs(ispin)%matrix, &
                                                   matrix_w(ispin)%matrix, matrix_s(1)%matrix)
                     END IF
                  ELSE
                     CALL calculate_w_matrix(mo_set=mo_set, &
                                             matrix_ks=matrix_ks(ispin)%matrix, &
                                             matrix_p=rho_ao(ispin)%matrix, &
                                             matrix_w=matrix_w(ispin)%matrix)
                  END IF
               ELSE
                  IF (scf_control%use_ot) THEN
                     CALL calculate_w_matrix_ot(mo_set, mo_derivs(ispin)%matrix, &
                                                matrix_w(ispin)%matrix, matrix_s(1)%matrix)
                  ELSE
                     CALL calculate_w_matrix(mo_set, matrix_w(ispin)%matrix)
                  END IF
               END IF
               ! if MP2 time to update the W matrix with the MP2 contribution
               IF (ASSOCIATED(qs_env%mp2_env)) THEN
                  CALL dbcsr_add(matrix_w(ispin)%matrix, matrix_w_mp2(ispin)%matrix, 1.0_dp, -1.0_dp)
               END IF
            END DO

         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_energies_compute_matrix_w

! **************************************************************************************************
!> \brief Refactoring of qs_energies_scf. Moves computation of properties
!>        into separate subroutine
!> \param qs_env ...
!> \par History
!>      05.2013 created [Florian Schiffmann]
! **************************************************************************************************

   SUBROUTINE qs_energies_properties(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies_properties', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: do_et
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: v_hartree_rspace
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: input, rest_b_section

      NULLIFY (atprop, energy, pw_env)
      CALL timeset(routineN, handle)

      NULLIFY (v_hartree_rspace%pw)

      ! atomic energies using Mulliken partition
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      input=input, &
                      atprop=atprop, &
                      energy=energy, &
                      v_hartree_rspace=v_hartree_rspace%pw, &
                      pw_env=pw_env)
      IF (atprop%energy) THEN
         CALL qs_energies_mulliken(qs_env)
         IF (.NOT. dft_control%qs_control%semi_empirical .AND. &
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
         CALL atprop_array_add(atprop%atener, atprop%atecc)
         CALL atprop_array_add(atprop%atener, atprop%ate1c)
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
         IF (dft_control%qs_control%becke_restraint) THEN
            rest_b_section => section_vals_get_subs_vals(input, "PROPERTIES%ET_COUPLING%BECKE_RESTRAINT_B")
            CALL read_becke_section(qs_control=dft_control%qs_control, &
                                    becke_section=rest_b_section)
         END IF
         CALL scf(qs_env=qs_env)
         qs_env%et_coupling%keep_matrix = .TRUE.

         CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., just_energy=.TRUE.)
         CALL calc_et_coupling(qs_env)
         IF (dft_control%qs_control%becke_restraint) THEN
            CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, &
                                      dft_control%qs_control%becke_control%becke_pot%pw)
            dft_control%qs_control%becke_control%need_pot = .TRUE.
         END IF
      END IF

      !Properties
      IF (dft_control%do_xas_calculation) THEN
         CALL xas(qs_env, dft_control)
      END IF

      ! Compute Linear Response properties as post-scf
      IF (.NOT. qs_env%linres_run) THEN
         CALL linres_calculation_low(qs_env)
      END IF

      IF (dft_control%tddfpt2_control%enabled) THEN
         CALL tddfpt(qs_env)
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

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies_mulliken', &
         routineP = moduleN//':'//routineN

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

      CHARACTER(len=*), PARAMETER :: routineN = 'atom_trace', routineP = moduleN//':'//routineN

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
         ENDIF
         btr = factor*mult*SUM(a_block*b_block)
         atrace(iblock_row) = atrace(iblock_row)+btr
         atrace(iblock_col) = atrace(iblock_col)+btr

      ENDDO
      CALL dbcsr_iterator_stop(iter)

   END SUBROUTINE atom_trace

! **************************************************************************************************
!> \brief Enters the mp2 part of cp2k
!> \param qs_env ...
!> \param calc_forces ...
! **************************************************************************************************

   SUBROUTINE qs_energies_mp2(qs_env, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies_mp2', &
         routineP = moduleN//':'//routineN

      LOGICAL                                            :: should_stop

! Compute MP2 energy

      IF (ASSOCIATED(qs_env%mp2_env)) THEN

         CALL external_control(should_stop, "MP2", target_time=qs_env%target_time, &
                               start_time=qs_env%start_time)

         CALL mp2_main(qs_env=qs_env, calc_forces=calc_forces)
      ENDIF

   END SUBROUTINE qs_energies_mp2

END MODULE qs_energy_utils
