!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Driver for the localization that should be general
!>      for all the methods available and all the definition of the
!>      spread functional
!>      Write centers, spread and cubes only if required and for the
!>      selected states
!>      The localized functions are copied in the standard mos array
!>      for the next use
!> \par History
!>      01.2008 Teodoro Laino [tlaino] - University of Zurich
!>        - Merging the two localization codes and updating to new structures
!>      04.2023 JGH Code isolation and refactoring
!> \author MI (04.2005)
! **************************************************************************************************
MODULE qs_loc_main
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_init_random, cp_fm_release, &
        cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_symmetric
   USE input_constants,                 ONLY: &
        do_loc_cpo_atomic, do_loc_cpo_random, do_loc_cpo_restart, do_loc_cpo_space_nmo, &
        do_loc_cpo_space_wan, op_loc_berry, op_loc_boys, op_loc_pipek, state_loc_list
   USE input_section_types,             ONLY: section_get_lval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_atomic_block,                 ONLY: calculate_atomic_block_dm
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_loc_methods,                  ONLY: optimize_loc_berry,&
                                              optimize_loc_pipek,&
                                              qs_print_cubes
   USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                              localized_wfn_control_type,&
                                              qs_loc_env_type
   USE qs_mo_methods,                   ONLY: make_basis_simple,&
                                              make_basis_sm
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_loc_main'

! *** Public ***
   PUBLIC :: qs_loc_driver

CONTAINS

! **************************************************************************************************
!> \brief set up the calculation of localized orbitals
!> \param qs_env ...
!> \param qs_loc_env ...
!> \param print_loc_section ...
!> \param myspin ...
!> \param ext_mo_coeff ...
!> \par History
!>      04.2005 created [MI]
!>      04.2023 refactored [JGH]
!> \author MI
! **************************************************************************************************
   SUBROUTINE qs_loc_driver(qs_env, qs_loc_env, print_loc_section, myspin, ext_mo_coeff)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_loc_env_type), POINTER                     :: qs_loc_env
      TYPE(section_vals_type), POINTER                   :: print_loc_section
      INTEGER, INTENT(IN)                                :: myspin
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL, TARGET     :: ext_mo_coeff

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_loc_driver'

      INTEGER                                            :: dim_op, handle, i, imo, imoloc, j, lb, &
                                                            loc_method, nao, nmosub, restricted, ub
      INTEGER, DIMENSION(:), POINTER                     :: ivec
      LOGICAL, SAVE                                      :: first_time = .TRUE.
      REAL(dp), DIMENSION(6)                             :: weights
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vecbuffer
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: moloc_coeff
      TYPE(cp_fm_type), DIMENSION(:, :), POINTER         :: op_fm_set
      TYPE(cp_fm_type), POINTER                          :: locorb
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: op_sm_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(localized_wfn_control_type), POINTER          :: localized_wfn_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: input, low_spin_roks_section

      CALL timeset(routineN, handle)
      NULLIFY (para_env, mos, dft_control)
      NULLIFY (cell, localized_wfn_control, moloc_coeff, op_sm_set, op_fm_set)
      qs_loc_env%first_time = first_time
      qs_loc_env%target_time = qs_env%target_time
      qs_loc_env%start_time = qs_env%start_time

      CALL get_qs_loc_env(qs_loc_env=qs_loc_env, &
                          localized_wfn_control=localized_wfn_control, &
                          moloc_coeff=moloc_coeff, op_sm_set=op_sm_set, op_fm_set=op_fm_set, cell=cell, &
                          weights=weights, dim_op=dim_op)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, &
                      para_env=para_env, mos=mos, input=input)

      !calculation of single occupied states to which unitary transformations should not be applied in LOW SPIN ROKS
      IF (dft_control%restricted) THEN
         low_spin_roks_section => section_vals_get_subs_vals(input, "DFT%LOW_SPIN_ROKS")
         CALL section_vals_val_get(low_spin_roks_section, "SPIN_CONFIGURATION", i_rep_val=1, i_vals=ivec)
         restricted = SIZE(ivec)
      ELSE
         restricted = 0
      END IF

      NULLIFY (locorb)
      IF (PRESENT(ext_mo_coeff)) THEN
         locorb => ext_mo_coeff
      ELSE
         CALL get_mo_set(mo_set=mos(myspin), mo_coeff=locorb)
      END IF

      loc_method = localized_wfn_control%localization_method

      nmosub = localized_wfn_control%nloc_states(myspin)
      IF (localized_wfn_control%operator_type == op_loc_berry) THEN
         ! Here we allocate op_fm_set with the RIGHT size for uks
         NULLIFY (tmp_fm_struct)
         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nmosub, &
                                  ncol_global=nmosub, para_env=para_env, &
                                  context=locorb%matrix_struct%context)
         !
         ALLOCATE (op_fm_set(2, dim_op))
         DO i = 1, dim_op
            DO j = 1, SIZE(op_fm_set, 1)
               CALL cp_fm_create(op_fm_set(j, i), tmp_fm_struct)
               CALL cp_fm_get_info(op_fm_set(j, i), nrow_global=nmosub)
               CALL cp_fm_set_all(op_fm_set(j, i), 0.0_dp)
            END DO
         END DO
         CALL cp_fm_struct_release(tmp_fm_struct)
      END IF

      IF (localized_wfn_control%do_mixed) THEN
         CALL loc_mixed_method(qs_env, qs_loc_env, print_loc_section, myspin, op_fm_set)
      ELSE
         SELECT CASE (localized_wfn_control%operator_type)
         CASE (op_loc_berry)
            CALL optimize_loc_berry(loc_method, qs_loc_env, moloc_coeff(myspin), op_sm_set, &
                                    op_fm_set, para_env, cell, weights, myspin, print_loc_section, &
                                    restricted=restricted)
         CASE (op_loc_boys)
            CPABORT("Boys localization not implemented")
         CASE (op_loc_pipek)
            CALL optimize_loc_pipek(qs_env, loc_method, qs_loc_env, moloc_coeff(myspin), &
                                    op_fm_set, myspin, print_loc_section)
         END SELECT
      END IF

      ! Here we dealloctate op_fm_set
      IF (localized_wfn_control%operator_type == op_loc_berry) THEN
         IF (ASSOCIATED(op_fm_set)) THEN
            DO i = 1, dim_op
               DO j = 1, SIZE(op_fm_set, 1)
                  CALL cp_fm_release(op_fm_set(j, i))
               END DO
            END DO
            DEALLOCATE (op_fm_set)
         END IF
      END IF

      ! give back the localized orbitals
      CALL get_mo_set(mo_set=mos(myspin), nao=nao)
      lb = localized_wfn_control%lu_bound_states(1, myspin)
      ub = localized_wfn_control%lu_bound_states(2, myspin)

      IF (localized_wfn_control%set_of_states == state_loc_list) THEN
         ALLOCATE (vecbuffer(1, nao))
         nmosub = SIZE(localized_wfn_control%loc_states, 1)
         imoloc = 0
         DO i = lb, ub
            ! Get the index in the subset
            imoloc = imoloc + 1
            ! Get the index in the full set
            imo = localized_wfn_control%loc_states(i, myspin)

            CALL cp_fm_get_submatrix(moloc_coeff(myspin), vecbuffer, 1, imoloc, &
                                     nao, 1, transpose=.TRUE.)
            CALL cp_fm_set_submatrix(locorb, vecbuffer, 1, imo, nao, 1, transpose=.TRUE.)
         END DO
         DEALLOCATE (vecbuffer)
      ELSE
         nmosub = localized_wfn_control%nloc_states(myspin)
         CALL cp_fm_to_fm(moloc_coeff(myspin), locorb, nmosub, 1, lb)
      END IF

      ! Write cube files if required
      IF (localized_wfn_control%print_cubes) THEN
         CALL loc_print(qs_env, qs_loc_env, moloc_coeff, myspin, print_loc_section)
      END IF
      first_time = .FALSE.

      CALL timestop(handle)

   END SUBROUTINE qs_loc_driver

! **************************************************************************************************
!> \brief set up the calculation of localized orbitals
!> \param qs_env ...
!> \param qs_loc_env ...
!> \param print_loc_section ...
!> \param myspin ...
!> \param op_fm_set ...
!> \par History
!>      04.2023 refactored [JGH]
!> \author MI
! **************************************************************************************************
   SUBROUTINE loc_mixed_method(qs_env, qs_loc_env, print_loc_section, myspin, op_fm_set)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_loc_env_type), POINTER                     :: qs_loc_env
      TYPE(section_vals_type), POINTER                   :: print_loc_section
      INTEGER, INTENT(IN)                                :: myspin
      TYPE(cp_fm_type), DIMENSION(:, :), POINTER         :: op_fm_set

      CHARACTER(len=*), PARAMETER                        :: routineN = 'loc_mixed_method'

      INTEGER                                            :: dim_op, handle, jspin, loc_method, nao, &
                                                            ndummy, nextra, ngextra, nguess, nmo, &
                                                            nmosub, norextra, restricted
      INTEGER, DIMENSION(2)                              :: nelectron_spin
      INTEGER, DIMENSION(:), POINTER                     :: ivec
      LOGICAL                                            :: do_ortho, has_unit_metric, &
                                                            my_guess_atomic, my_guess_wan
      REAL(dp), DIMENSION(6)                             :: weights
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: tmp_mat
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type)                                   :: mos_guess, tmp_fm, tmp_fm_1, vectors_2
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: moloc_coeff
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_rmpv
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, op_sm_set
      TYPE(dbcsr_type), POINTER                          :: refmatrix, tmatrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(localized_wfn_control_type), POINTER          :: localized_wfn_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: input, low_spin_roks_section

      CALL timeset(routineN, handle)

      NULLIFY (moloc_coeff, op_sm_set)
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, para_env=para_env, mos=mos, input=input)

      !calculation of single occupied states to which unitary transformations should not be applied in LOW SPIN ROKS
      IF (dft_control%restricted) THEN
         low_spin_roks_section => section_vals_get_subs_vals(input, "DFT%LOW_SPIN_ROKS")
         CALL section_vals_val_get(low_spin_roks_section, "SPIN_CONFIGURATION", i_rep_val=1, i_vals=ivec)
         restricted = SIZE(ivec)
      ELSE
         restricted = 0
      END IF

      CALL get_qs_loc_env(qs_loc_env=qs_loc_env, &
                          localized_wfn_control=localized_wfn_control, &
                          moloc_coeff=moloc_coeff, op_sm_set=op_sm_set, cell=cell, &
                          weights=weights, dim_op=dim_op)

      CALL get_mo_set(mo_set=mos(myspin), nao=nao, nmo=nmo)
      loc_method = localized_wfn_control%localization_method
      nmosub = localized_wfn_control%nloc_states(myspin)

      CPASSERT(localized_wfn_control%operator_type == op_loc_berry)
      CPASSERT(localized_wfn_control%do_mixed)

      my_guess_atomic = .FALSE.
      ! SGh-wan: if atomic guess and do_mixed and nextra > 0
      ! read CPO_GUESS; CASE ATOMIC / RESTART / RANDOM (0/1/2)
      ! read CPO_GUESS_SPACE if CASE ATOMIC; CASE ALL / WAN
      nextra = localized_wfn_control%nextra
      IF (nextra > 0) THEN
         my_guess_atomic = .TRUE.
         my_guess_wan = .FALSE.
         do_ortho = .TRUE.
         SELECT CASE (localized_wfn_control%coeff_po_guess)

         CASE (do_loc_cpo_atomic)
            my_guess_atomic = .TRUE.
            NULLIFY (atomic_kind_set, qs_kind_set, particle_set, matrix_s_kp, sab_orb, p_rmpv, &
                     refmatrix, tmatrix)
            CALL get_qs_env(qs_env=qs_env, &
                            atomic_kind_set=atomic_kind_set, &
                            qs_kind_set=qs_kind_set, &
                            particle_set=particle_set, &
                            matrix_s_kp=matrix_s_kp, &
                            has_unit_metric=has_unit_metric, &
                            nelectron_spin=nelectron_spin, &
                            sab_orb=sab_orb)

            refmatrix => matrix_s_kp(1, 1)%matrix
            ! create p_rmpv
            CALL dbcsr_allocate_matrix_set(p_rmpv, dft_control%nspins)
            DO jspin = 1, dft_control%nspins
               ALLOCATE (p_rmpv(jspin)%matrix)
               tmatrix => p_rmpv(jspin)%matrix
               CALL dbcsr_create(matrix=tmatrix, template=refmatrix, &
                                 matrix_type=dbcsr_type_symmetric, nze=0)
               CALL cp_dbcsr_alloc_block_from_nbl(tmatrix, sab_orb)
               CALL dbcsr_set(tmatrix, 0.0_dp)
            END DO
            CALL calculate_atomic_block_dm(p_rmpv, refmatrix, atomic_kind_set, qs_kind_set, &
                                           dft_control%nspins, nelectron_spin, 0, para_env)
         CASE (do_loc_cpo_restart)
            my_guess_atomic = .FALSE.
            my_guess_wan = .TRUE.
         CASE (do_loc_cpo_random)
            my_guess_atomic = .FALSE.
         END SELECT

         norextra = nmo - nmosub
         CALL get_mo_set(mo_set=mos(myspin), mo_coeff=mo_coeff)
         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nao, &
                                  ncol_global=norextra, para_env=para_env, context=mo_coeff%matrix_struct%context)
         CALL cp_fm_create(vectors_2, tmp_fm_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
         ALLOCATE (tmp_mat(nao, norextra))
         CALL cp_fm_get_submatrix(mo_coeff, tmp_mat, 1, nmosub + 1)
         CALL cp_fm_set_submatrix(vectors_2, tmp_mat)
         DEALLOCATE (tmp_mat)

         ! if guess "atomic" generate MOs based on atomic densities and
         ! pass on to optimize_loc_berry
         IF (my_guess_atomic .OR. my_guess_wan) THEN

            SELECT CASE (localized_wfn_control%coeff_po_guess_mo_space)

            CASE (do_loc_cpo_space_wan)
               ndummy = nmosub
            CASE (do_loc_cpo_space_nmo)
               ndummy = nmo
               do_ortho = .FALSE.

            END SELECT

            CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nao, &
                                     ncol_global=ndummy, para_env=para_env, &
                                     context=mo_coeff%matrix_struct%context)
            CALL cp_fm_create(mos_guess, tmp_fm_struct)
            CALL cp_fm_set_all(mos_guess, 0.0_dp)

            IF (my_guess_atomic) THEN
               CALL cp_fm_create(tmp_fm, tmp_fm_struct)
               CALL cp_fm_create(tmp_fm_1, tmp_fm_struct)
               CALL cp_fm_set_all(tmp_fm, 0.0_dp)
               CALL cp_fm_set_all(tmp_fm_1, 0.0_dp)
               CALL cp_fm_init_random(tmp_fm, ndummy)
               IF (has_unit_metric) THEN
                  CALL cp_fm_to_fm(tmp_fm, tmp_fm_1)
               ELSE
                  ! PS*C(:,1:nomo)+C(:,nomo+1:nmo) (nomo=NINT(nelectron/maxocc))
                  CALL cp_dbcsr_sm_fm_multiply(refmatrix, tmp_fm, tmp_fm_1, ndummy)
               END IF
               CALL cp_dbcsr_sm_fm_multiply(p_rmpv(myspin)%matrix, tmp_fm_1, mos_guess, ndummy)
               CALL cp_fm_release(tmp_fm)
               CALL cp_fm_release(tmp_fm_1)
               CALL cp_fm_struct_release(tmp_fm_struct)
            ELSEIF (my_guess_wan) THEN
               nguess = localized_wfn_control%nguess(myspin)
               ALLOCATE (tmp_mat(nao, nguess))
               CALL cp_fm_get_submatrix(moloc_coeff(myspin), tmp_mat, 1, 1, nao, nguess)
               CALL cp_fm_set_submatrix(mos_guess, tmp_mat, 1, 1, nao, nguess)
               DEALLOCATE (tmp_mat)
               ngextra = nmosub - nguess
               !WRITE(*,*) 'nguess, ngextra = ', nguess, ngextra
               CALL cp_fm_struct_release(tmp_fm_struct)
               IF (ngextra > 0) THEN
                  CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nao, &
                                           ncol_global=ngextra, para_env=para_env, &
                                           context=mo_coeff%matrix_struct%context)
                  CALL cp_fm_create(tmp_fm, tmp_fm_struct)
                  CALL cp_fm_init_random(tmp_fm, ngextra)
                  ALLOCATE (tmp_mat(nao, ngextra))
                  CALL cp_fm_get_submatrix(tmp_fm, tmp_mat, 1, 1, nao, ngextra)
                  CALL cp_fm_set_submatrix(mos_guess, tmp_mat, 1, nguess + 1, nao, ngextra)
                  DEALLOCATE (tmp_mat)
                  CALL cp_fm_release(tmp_fm)
                  CALL cp_fm_struct_release(tmp_fm_struct)
               ELSE
                  do_ortho = .FALSE.
               END IF
               ALLOCATE (tmp_mat(nao, nmosub))
               CALL cp_fm_get_submatrix(mo_coeff, tmp_mat, 1, 1, nao, nmosub)
               CALL cp_fm_set_submatrix(moloc_coeff(myspin), tmp_mat)
               DEALLOCATE (tmp_mat)
            END IF

            IF (do_ortho) THEN
               IF ((my_guess_atomic) .OR. (my_guess_wan)) THEN
                        !! and ortho the result
                  IF (has_unit_metric) THEN
                     CALL make_basis_simple(mos_guess, ndummy)
                  ELSE
                     CALL make_basis_sm(mos_guess, ndummy, refmatrix)
                  END IF
               END IF
            END IF

            CALL optimize_loc_berry(loc_method, qs_loc_env, moloc_coeff(myspin), op_sm_set, &
                                    op_fm_set, para_env, cell, weights, myspin, print_loc_section, &
                                    restricted=restricted, &
                                    nextra=nextra, nmo=nmo, vectors_2=vectors_2, guess_mos=mos_guess)
            CALL cp_fm_release(mos_guess)
         ELSE
            CALL optimize_loc_berry(loc_method, qs_loc_env, moloc_coeff(myspin), op_sm_set, &
                                    op_fm_set, para_env, cell, weights, myspin, print_loc_section, &
                                    restricted=restricted, &
                                    nextra=nextra, nmo=nmo, vectors_2=vectors_2)
         END IF
         CALL cp_fm_release(vectors_2)
         IF (my_guess_atomic) CALL dbcsr_deallocate_matrix_set(p_rmpv)
      ELSE
         CALL optimize_loc_berry(loc_method, qs_loc_env, moloc_coeff(myspin), op_sm_set, &
                                 op_fm_set, para_env, cell, weights, myspin, print_loc_section, &
                                 restricted=restricted, nextra=0)
      END IF

      CALL timestop(handle)

   END SUBROUTINE loc_mixed_method

! **************************************************************************************************
!> \brief printing of Cube files of localized orbitals
!> \param qs_env ...
!> \param qs_loc_env ...
!> \param moloc_coeff ...
!> \param ispin ...
!> \param print_loc_section ...
! **************************************************************************************************
   SUBROUTINE loc_print(qs_env, qs_loc_env, moloc_coeff, ispin, print_loc_section)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_loc_env_type), POINTER                     :: qs_loc_env
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: moloc_coeff
      INTEGER, INTENT(IN), OPTIONAL                      :: ispin
      TYPE(section_vals_type), POINTER                   :: print_loc_section

      CHARACTER(LEN=default_string_length)               :: my_pos
      INTEGER                                            :: i, ir, istate, j, jstate, n_rep, ncubes, &
                                                            nmo
      INTEGER, DIMENSION(:), POINTER                     :: bounds, list, list_cubes
      LOGICAL                                            :: append_cube, list_cubes_setup
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: centers
      TYPE(localized_wfn_control_type), POINTER          :: localized_wfn_control
      TYPE(section_vals_type), POINTER                   :: print_key

      list_cubes_setup = .FALSE.
      NULLIFY (bounds, list, list_cubes)

      CALL get_qs_loc_env(qs_loc_env=qs_loc_env, &
                          localized_wfn_control=localized_wfn_control)

      ! Provides boundaries of MOs
      CALL section_vals_val_get(print_loc_section, "WANNIER_CUBES%CUBES_LU_BOUNDS", &
                                i_vals=bounds)
      ncubes = bounds(2) - bounds(1) + 1
      IF (ncubes > 0) THEN
         list_cubes_setup = .TRUE.
         ALLOCATE (list_cubes(ncubes))
         DO ir = 1, ncubes
            list_cubes(ir) = bounds(1) + (ir - 1)
         END DO
      END IF

      ! Provides the list of MOs
      CALL section_vals_val_get(print_loc_section, "WANNIER_CUBES%CUBES_LIST", &
                                n_rep_val=n_rep)
      IF (.NOT. list_cubes_setup) THEN
         ncubes = 0
         DO ir = 1, n_rep
            CALL section_vals_val_get(print_loc_section, "WANNIER_CUBES%CUBES_LIST", &
                                      i_rep_val=ir, i_vals=list)
            IF (ASSOCIATED(list)) THEN
               CALL reallocate(list_cubes, 1, ncubes + SIZE(list))
               DO i = 1, SIZE(list)
                  list_cubes(i + ncubes) = list(i)
               END DO
               ncubes = ncubes + SIZE(list)
            END IF
         END DO
         IF (ncubes > 0) list_cubes_setup = .TRUE.
      END IF

      ! Full list of Mos
      IF (.NOT. list_cubes_setup) THEN
         list_cubes_setup = .TRUE.
         ncubes = localized_wfn_control%nloc_states(1)
         IF (ncubes > 0) THEN
            ALLOCATE (list_cubes(ncubes))
         END IF
         DO i = 1, ncubes
            list_cubes(i) = i
         END DO
      END IF

      ncubes = SIZE(list_cubes)
      ncubes = MIN(ncubes, nmo)
      ALLOCATE (centers(6, ncubes))
      DO i = 1, ncubes
         istate = list_cubes(i)
         DO j = 1, localized_wfn_control%nloc_states(ispin)
            jstate = localized_wfn_control%loc_states(j, ispin)
            IF (istate == jstate) THEN
               centers(1:6, i) = localized_wfn_control%centers_set(ispin)%array(1:6, j)
               EXIT
            END IF
         END DO
      END DO ! ncubes

      ! Real call for dumping the cube files
      print_key => section_vals_get_subs_vals(print_loc_section, "WANNIER_CUBES")
      append_cube = section_get_lval(print_loc_section, "WANNIER_CUBES%APPEND")
      my_pos = "REWIND"
      IF (append_cube) THEN
         my_pos = "APPEND"
      END IF

      CALL qs_print_cubes(qs_env, moloc_coeff(ispin), ncubes, list_cubes, centers, &
                          print_key, "loc"//TRIM(ADJUSTL(qs_loc_env%tag_mo)), &
                          ispin=ispin, file_position=my_pos)

      DEALLOCATE (centers)
      DEALLOCATE (list_cubes)

   END SUBROUTINE loc_print

END MODULE qs_loc_main
