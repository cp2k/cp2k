!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Second order perturbation correction to XAS_TDP spectra (i.e. shift)
!> \author A. Bussy (01.2020)
! **************************************************************************************************

MODULE xas_tdp_correction
   USE admm_types,                      ONLY: admm_type
   USE admm_utils,                      ONLY: admm_correct_for_eigenvalues
   USE bibliography,                    ONLY: Bussy2021b,&
                                              Shigeta2001,&
                                              cite_reference
   USE cp_array_utils,                  ONLY: cp_1d_i_p_type,&
                                              cp_1d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_p_type,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_get_info, dbcsr_p_type, &
        dbcsr_release, dbcsr_scalar, dbcsr_type
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_contract, dbcsr_t_copy, dbcsr_t_copy_matrix_to_tensor, dbcsr_t_create, &
        dbcsr_t_default_distvec, dbcsr_t_destroy, dbcsr_t_distribution_destroy, &
        dbcsr_t_distribution_new, dbcsr_t_distribution_type, dbcsr_t_finalize, dbcsr_t_get_block, &
        dbcsr_t_get_info, dbcsr_t_iterator_blocks_left, dbcsr_t_iterator_next_block, &
        dbcsr_t_iterator_start, dbcsr_t_iterator_stop, dbcsr_t_iterator_type, &
        dbcsr_t_pgrid_create, dbcsr_t_pgrid_destroy, dbcsr_t_pgrid_type, dbcsr_t_put_block, &
        dbcsr_t_type
   USE hfx_admm_utils,                  ONLY: create_admm_xc_section
   USE input_section_types,             ONLY: section_vals_create,&
                                              section_vals_get_subs_vals,&
                                              section_vals_release,&
                                              section_vals_retain,&
                                              section_vals_set_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush
   USE mathlib,                         ONLY: complex_diag
   USE message_passing,                 ONLY: mp_sum
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              duplicate_mo_set,&
                                              get_mo_set,&
                                              mo_set_p_type,&
                                              reassign_allocated_mos
   USE util,                            ONLY: get_limit
   USE xas_tdp_kernel,                  ONLY: contract2_AO_to_doMO,&
                                              ri_all_blocks_mm
   USE xas_tdp_types,                   ONLY: donor_state_type,&
                                              xas_tdp_control_type,&
                                              xas_tdp_env_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_correction'

   PUBLIC :: gw2x_shift, get_soc_splitting

CONTAINS

! **************************************************************************************************
!> \brief Computes the ionization potential using the GW2X method of Shigeta et. al. The result cam
!>        be used for XAS correction (shift) or XPS directly.
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE GW2X_shift(donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'GW2X_shift', routineP = moduleN//':'//routineN

      INTEGER :: ex_idx, exat, first_domo(2), handle, i, ido_mo, iloc, ilocat, ispin, jspin, &
         locat, nao, natom, ndo_mo, nhomo(2), nlumo(2), nonloc, nspins, start_sgf
      INTEGER, DIMENSION(:), POINTER                     :: nsgf_blk
      LOGICAL                                            :: pseudo_canonical
      REAL(dp)                                           :: og_hfx_frac
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: contract_coeffs_backup
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_1d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: homo_evals, lumo_evals
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: all_coeffs, homo_coeffs, lumo_coeffs
      TYPE(cp_fm_struct_p_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: all_struct, homo_struct, lumo_struct
      TYPE(cp_fm_struct_type), POINTER                   :: hoho_struct, lulu_struct
      TYPE(cp_fm_type), POINTER                          :: hoho_fock, hoho_work, homo_work, &
                                                            lulu_fock, lulu_work, lumo_work, &
                                                            mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dbcsr_work, fock_matrix, matrix_ks
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:)      :: ja_X, oI_Y
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:, :)   :: mo_template
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: backup_mos, mos
      TYPE(section_vals_type), POINTER                   :: xc_fun_empty, xc_fun_original, xc_section

      NULLIFY (xc_fun_empty, xc_fun_original, xc_section, backup_mos, mos, dft_control, dbcsr_work, &
               fock_matrix, matrix_ks, para_env, hoho_fock, hoho_work, homo_work, lulu_fock, lulu_work, &
               lumo_work, mo_coeff, blacs_env, nsgf_blk)

      CALL cite_reference(Shigeta2001)
      CALL cite_reference(Bussy2021b)

      CALL timeset(routineN, handle)

      !The GW2X correction we want to compute goes like this, where omega is the corrected epsilon_I:
      !omega = eps_I + 0.5 * sum_ajk |<Ia||jk>|^2/(omega + eps_a - eps_j - eps_k)
      !              + 0.5 * sum_jab |<Ij||ab>|^2/(omega + eps_j - eps_a - eps_b)
      ! j,k denote occupied spin-orbitals and a,b denote virtual spin orbitals

      !The strategy is the following (we assume restricted closed-shell):
      !1) Get the LUMOs from xas_tdp_env
      !2) Get the HOMOs from qs_env
      !3) Compute or fetch the generalize Fock matric
      !4) Diagonalize it in the subspace of HOMOs and LUMOs (or just take diagonal matrix elements)
      !5) Build the full HOMO-LUMO basis that we will use and compute eigenvalues
      !6) Iterate over GW2X steps to compute the self energy

      !We implement 2 approaches => diagonal elements of Fock matrix with original MOs and
      !pseudo-canonical MOs
      pseudo_canonical = xas_tdp_control%pseudo_canonical

      !Get donor state info
      ndo_mo = donor_state%ndo_mo
      nspins = 1; IF (xas_tdp_control%do_uks .OR. xas_tdp_control%do_roks) nspins = 2

      !1) Get the LUMO coefficients from the xas_tdp_env, that have been precomputed

      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, mos=mos, para_env=para_env, &
                      blacs_env=blacs_env, natom=natom)

      ALLOCATE (lumo_struct(nspins), lumo_coeffs(nspins))

      DO ispin = 1, nspins
         CALL get_mo_set(mos(ispin)%mo_set, homo=nhomo(ispin), nao=nao)
         nlumo(ispin) = nao - nhomo(ispin)

         CALL cp_fm_struct_create(lumo_struct(ispin)%struct, para_env=para_env, context=blacs_env, &
                                  ncol_global=nlumo(ispin), nrow_global=nao)

         CALL cp_fm_create(lumo_coeffs(ispin)%matrix, lumo_struct(ispin)%struct)
         CALL cp_fm_to_fm(xas_tdp_env%lumo_evecs(ispin)%matrix, lumo_coeffs(ispin)%matrix)
      END DO

      !2) get the HOMO coeffs. Reminder: keep all non-localized MOs + those localized on core atom
      !   For this to work, it is assumed that the LOCALIZE keyword is used
      ALLOCATE (homo_struct(nspins), homo_coeffs(nspins))

      DO ispin = 1, nspins
         nonloc = nhomo(ispin) - xas_tdp_control%n_search
         exat = donor_state%at_index
         ex_idx = MINLOC(ABS(xas_tdp_env%ex_atom_indices - exat), 1)
         locat = COUNT(xas_tdp_env%mos_of_ex_atoms(:, ex_idx, ispin) == 1)

         CALL cp_fm_struct_create(homo_struct(ispin)%struct, para_env=para_env, context=blacs_env, &
                                  ncol_global=locat + nonloc, nrow_global=nao)
         CALL cp_fm_create(homo_coeffs(ispin)%matrix, homo_struct(ispin)%struct)

         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_to_fm_submat(mo_coeff, homo_coeffs(ispin)%matrix, nrow=nao, ncol=nonloc, s_firstrow=1, &
                                 s_firstcol=xas_tdp_control%n_search + 1, t_firstrow=1, t_firstcol=locat + 1)

         !this bit is taken from xas_tdp_methods
         ilocat = 1
         DO iloc = 1, xas_tdp_control%n_search
            IF (xas_tdp_env%mos_of_ex_atoms(iloc, ex_idx, ispin) == -1) CYCLE
            CALL cp_fm_to_fm_submat(mo_coeff, homo_coeffs(ispin)%matrix, nrow=nao, ncol=1, s_firstrow=1, &
                                    s_firstcol=iloc, t_firstrow=1, t_firstcol=ilocat)
            !keep track of donor MO index
            IF (iloc == donor_state%mo_indices(1, ispin)) first_domo(ispin) = ilocat !first donor MO

            ilocat = ilocat + 1
         END DO
         nhomo(ispin) = locat + nonloc
      END DO

      !3) Computing the generalized Fock Matrix, if not there already
      IF (ASSOCIATED(xas_tdp_env%fock_matrix)) THEN
         fock_matrix => xas_tdp_env%fock_matrix
      ELSE

         ALLOCATE (xas_tdp_env%fock_matrix(nspins))
         fock_matrix => xas_tdp_env%fock_matrix

         ! remove the xc_functionals and set HF fraction to 1
         xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
         xc_fun_original => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
         CALL section_vals_retain(xc_fun_original)
         CALL section_vals_create(xc_fun_empty, xc_fun_original%section)
         CALL section_vals_set_subs_vals(xc_section, "XC_FUNCTIONAL", xc_fun_empty)
         CALL section_vals_release(xc_fun_empty)
         og_hfx_frac = qs_env%x_data(1, 1)%general_parameter%fraction
         qs_env%x_data(:, :)%general_parameter%fraction = 1.0_dp

         !In case of ADMM, we need to re-create the admm XC section for the new hfx_fraction
         !We also need to make a backup of the MOs as theiy are modified
         CALL get_qs_env(qs_env, dft_control=dft_control, admm_env=admm_env)
         IF (dft_control%do_admm) THEN
            IF (ASSOCIATED(admm_env%xc_section_primary)) CALL section_vals_release(admm_env%xc_section_primary)
            IF (ASSOCIATED(admm_env%xc_section_aux)) CALL section_vals_release(admm_env%xc_section_aux)
            CALL create_admm_xc_section(qs_env, xc_section, admm_env)

            ALLOCATE (backup_mos(SIZE(mos)))
            DO i = 1, SIZE(mos)
               CALL duplicate_mo_set(backup_mos(i)%mo_set, mos(i)%mo_set)
            END DO
         END IF

         ALLOCATE (dbcsr_work(nspins))
         DO ispin = 1, nspins
            ALLOCATE (dbcsr_work(ispin)%matrix)
            CALL dbcsr_copy(dbcsr_work(ispin)%matrix, matrix_ks(ispin)%matrix)
         END DO

         !both spins treated internally
         CALL qs_ks_build_kohn_sham_matrix(qs_env, calculate_forces=.FALSE., just_energy=.FALSE.)

         DO ispin = 1, nspins
            ALLOCATE (fock_matrix(ispin)%matrix)
            CALL dbcsr_copy(fock_matrix(ispin)%matrix, matrix_ks(ispin)%matrix, name="FOCK MATRIX")
            CALL dbcsr_release(matrix_ks(ispin)%matrix)
            CALL dbcsr_copy(matrix_ks(ispin)%matrix, dbcsr_work(ispin)%matrix)
         END DO
         CALL dbcsr_deallocate_matrix_set(dbcsr_work)

         !In case of ADMM, we want to correct for eigenvalues
         IF (dft_control%do_admm) THEN
            DO ispin = 1, nspins
               CALL admm_correct_for_eigenvalues(ispin, admm_env, fock_matrix(ispin)%matrix)
            END DO
         END IF

         !restore xc and HF fraction
         CALL section_vals_set_subs_vals(xc_section, "XC_FUNCTIONAL", xc_fun_original)
         CALL section_vals_release(xc_fun_original)
         qs_env%x_data(:, :)%general_parameter%fraction = og_hfx_frac

         IF (dft_control%do_admm) THEN
            IF (ASSOCIATED(admm_env%xc_section_primary)) CALL section_vals_release(admm_env%xc_section_primary)
            IF (ASSOCIATED(admm_env%xc_section_aux)) CALL section_vals_release(admm_env%xc_section_aux)
            CALL create_admm_xc_section(qs_env, xc_section, admm_env)

            DO i = 1, SIZE(mos)
               CALL reassign_allocated_mos(mos(i)%mo_set, backup_mos(i)%mo_set)
               CALL deallocate_mo_set(backup_mos(i)%mo_set)
            END DO
            DEALLOCATE (backup_mos)
         END IF
      END IF

      !4,5) Build pseudo-canonical MOs if needed + get related Fock matrix elements
      ALLOCATE (all_struct(nspins), all_coeffs(nspins))
      ALLOCATE (homo_evals(nspins), lumo_evals(nspins))
      CALL dbcsr_get_info(matrix_ks(1)%matrix, row_blk_size=nsgf_blk)
      ALLOCATE (contract_coeffs_backup(nsgf_blk(exat), nspins*ndo_mo))

      DO ispin = 1, nspins
         CALL cp_fm_struct_create(hoho_struct, para_env=para_env, context=blacs_env, &
                                  ncol_global=nhomo(ispin), nrow_global=nhomo(ispin))
         CALL cp_fm_struct_create(lulu_struct, para_env=para_env, context=blacs_env, &
                                  ncol_global=nlumo(ispin), nrow_global=nlumo(ispin))

         CALL cp_fm_create(hoho_work, hoho_struct)
         CALL cp_fm_create(lulu_work, lulu_struct)
         CALL cp_fm_create(homo_work, homo_struct(ispin)%struct)
         CALL cp_fm_create(lumo_work, lumo_struct(ispin)%struct)

         IF (pseudo_canonical) THEN
            !That is where we rotate the MOs to make them pseudo canonical
            !The eigenvalues we get from the diagonalization

            !The Fock matrix in the HOMO subspace
            CALL cp_fm_create(hoho_fock, hoho_struct)
            NULLIFY (homo_evals(ispin)%array)
            ALLOCATE (homo_evals(ispin)%array(nhomo(ispin)))
            CALL cp_dbcsr_sm_fm_multiply(fock_matrix(ispin)%matrix, homo_coeffs(ispin)%matrix, &
                                         homo_work, ncol=nhomo(ispin))
            CALL cp_gemm('T', 'N', nhomo(ispin), nhomo(ispin), nao, 1.0_dp, homo_coeffs(ispin)%matrix, &
                         homo_work, 0.0_dp, hoho_fock)

            !diagonalize and get pseudo-canonical MOs
            CALL choose_eigv_solver(hoho_fock, hoho_work, homo_evals(ispin)%array)
            CALL cp_gemm('N', 'N', nao, nhomo(ispin), nhomo(ispin), 1.0_dp, homo_coeffs(ispin)%matrix, &
                         hoho_work, 0.0_dp, homo_work)
            CALL cp_fm_to_fm(homo_work, homo_coeffs(ispin)%matrix)

            !overwrite the donor_state's contract coeffs with those
            contract_coeffs_backup(:, (ispin - 1)*ndo_mo + 1:ispin*ndo_mo) = &
               donor_state%contract_coeffs(:, (ispin - 1)*ndo_mo + 1:ispin*ndo_mo)
            start_sgf = SUM(nsgf_blk(1:exat - 1)) + 1
            CALL cp_fm_get_submatrix(homo_coeffs(ispin)%matrix, &
                                     donor_state%contract_coeffs(:, (ispin - 1)*ndo_mo + 1:ispin*ndo_mo), &
                                     start_row=start_sgf, start_col=first_domo(ispin), &
                                     n_rows=nsgf_blk(exat), n_cols=ndo_mo)

            !do the same for the pseudo-LUMOs
            CALL cp_fm_create(lulu_fock, lulu_struct)
            NULLIFY (lumo_evals(ispin)%array)
            ALLOCATE (lumo_evals(ispin)%array(nlumo(ispin)))
            CALL cp_dbcsr_sm_fm_multiply(fock_matrix(ispin)%matrix, lumo_coeffs(ispin)%matrix, &
                                         lumo_work, ncol=nlumo(ispin))
            CALL cp_gemm('T', 'N', nlumo(ispin), nlumo(ispin), nao, 1.0_dp, lumo_coeffs(ispin)%matrix, &
                         lumo_work, 0.0_dp, lulu_fock)

            !diagonalize and get pseudo-canonical MOs
            CALL choose_eigv_solver(lulu_fock, lulu_work, lumo_evals(ispin)%array)
            CALL cp_gemm('N', 'N', nao, nlumo(ispin), nlumo(ispin), 1.0_dp, lumo_coeffs(ispin)%matrix, &
                         lulu_work, 0.0_dp, lumo_work)
            CALL cp_fm_to_fm(lumo_work, lumo_coeffs(ispin)%matrix)

            CALL cp_fm_release(lulu_fock)
            CALL cp_fm_release(hoho_fock)

         ELSE !using the generalized Fock matrix diagonal elements

            !Compute their Fock matrix diagonal
            ALLOCATE (homo_evals(ispin)%array(nhomo(ispin)))
            CALL cp_dbcsr_sm_fm_multiply(fock_matrix(ispin)%matrix, homo_coeffs(ispin)%matrix, &
                                         homo_work, ncol=nhomo(ispin))
            CALL cp_gemm('T', 'N', nhomo(ispin), nhomo(ispin), nao, 1.0_dp, homo_coeffs(ispin)%matrix, &
                         homo_work, 0.0_dp, hoho_work)
            CALL cp_fm_get_diag(hoho_work, homo_evals(ispin)%array)

            ALLOCATE (lumo_evals(ispin)%array(nlumo(ispin)))
            CALL cp_dbcsr_sm_fm_multiply(fock_matrix(ispin)%matrix, lumo_coeffs(ispin)%matrix, &
                                         lumo_work, ncol=nlumo(ispin))
            CALL cp_gemm('T', 'N', nlumo(ispin), nlumo(ispin), nao, 1.0_dp, lumo_coeffs(ispin)%matrix, &
                         lumo_work, 0.0_dp, lulu_work)
            CALL cp_fm_get_diag(lulu_work, lumo_evals(ispin)%array)

         END IF
         CALL cp_fm_release(homo_work)
         CALL cp_fm_release(hoho_work)
         CALL cp_fm_struct_release(hoho_struct)
         CALL cp_fm_release(lumo_work)
         CALL cp_fm_release(lulu_work)
         CALL cp_fm_struct_release(lulu_struct)

         !Put back homo and lumo coeffs together, to fit tensor structure
         CALL cp_fm_struct_create(all_struct(ispin)%struct, para_env=para_env, context=blacs_env, &
                                  ncol_global=nhomo(ispin) + nlumo(ispin), nrow_global=nao)
         CALL cp_fm_create(all_coeffs(ispin)%matrix, all_struct(ispin)%struct)
         CALL cp_fm_to_fm(homo_coeffs(ispin)%matrix, all_coeffs(ispin)%matrix, ncol=nhomo(ispin), &
                          source_start=1, target_start=1)
         CALL cp_fm_to_fm(lumo_coeffs(ispin)%matrix, all_coeffs(ispin)%matrix, ncol=nlumo(ispin), &
                          source_start=1, target_start=nhomo(ispin) + 1)

      END DO !ispin

      !get semi-contracted tensor (AOs to MOs, keep RI uncontracted)
      CALL contract_AOs_to_MOs(ja_X, oI_Y, mo_template, all_coeffs, nhomo, nlumo, &
                               donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      !intermediate clean-up
      DO ispin = 1, nspins
         CALL cp_fm_release(all_coeffs(ispin)%matrix)
         CALL cp_fm_release(homo_coeffs(ispin)%matrix)
         CALL cp_fm_release(lumo_coeffs(ispin)%matrix)
         CALL cp_fm_struct_release(all_struct(ispin)%struct)
         CALL cp_fm_struct_release(lumo_struct(ispin)%struct)
         CALL cp_fm_struct_release(homo_struct(ispin)%struct)
      END DO

      !6) GW2X iterations

      IF (nspins == 1) THEN
         !restricted-closed shell: only alpha spin
         CALL GW2X_rcs_iterations(first_domo(1), ja_X(1), oI_Y, mo_template(1, 1), homo_evals(1)%array, &
                                  lumo_evals(1)%array, donor_state, xas_tdp_control, qs_env)
      ELSE
         !open-shell, need both spins
         CALL GW2X_os_iterations(first_domo, ja_X, oI_Y, mo_template, homo_evals, lumo_evals, &
                                 donor_state, xas_tdp_control, qs_env)
      END IF

      !restore proper contract_coeffs
      IF (pseudo_canonical) THEN
         donor_state%contract_coeffs(:, :) = contract_coeffs_backup(:, :)
      END IF

      !Final clean-up
      DO ido_mo = 1, nspins*ndo_mo
         CALL dbcsr_t_destroy(oI_Y(ido_mo))
      END DO
      DO ispin = 1, nspins
         CALL dbcsr_t_destroy(ja_X(ispin))
         DEALLOCATE (homo_evals(ispin)%array)
         DEALLOCATE (lumo_evals(ispin)%array)
         DO jspin = 1, nspins
            CALL dbcsr_t_destroy(mo_template(ispin, jspin))
         END DO
      END DO
      DEALLOCATE (oI_Y, homo_evals, lumo_evals)

      CALL timestop(handle)

   END SUBROUTINE GW2X_shift

! **************************************************************************************************
!> \brief Preforms the GW2X iterations in the restricted-closed shell formalism according to the
!>        Newton-Raphson method
!> \param first_domo index of the first core donor MO to consider
!> \param ja_X semi-contracted tensor with j: occupied MO, a: virtual MO, X: RI basis element
!> \param oI_Y semi-contracted tensors with o: all MOs, I donor core MO, Y: RI basis element
!> \param mo_template tensor template for fully MO contracted tensor
!> \param homo_evals ...
!> \param lumo_evals ...
!> \param donor_state ...
!> \param xas_tdp_control ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE GW2X_rcs_iterations(first_domo, ja_X, oI_Y, mo_template, homo_evals, lumo_evals, &
                                  donor_state, xas_tdp_control, qs_env)

      INTEGER, INTENT(IN)                                :: first_domo
      TYPE(dbcsr_t_type), INTENT(inout)                  :: ja_X
      TYPE(dbcsr_t_type), DIMENSION(:), INTENT(inout)    :: oI_Y
      TYPE(dbcsr_t_type), INTENT(inout)                  :: mo_template
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: homo_evals, lumo_evals
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'GW2X_rcs_iterations'

      INTEGER :: batch_size, bounds_1d(2), bounds_2d(2, 2), handle, i, ibatch, ido_mo, iloop, &
         max_iter, nbatch_occ, nbatch_virt, nblk_occ, nblk_virt, nblks(3), ndo_mo, nhomo, nlumo, &
         occ_bo(2), output_unit, tmp_sum, virt_bo(2)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: mo_blk_size
      REAL(dp)                                           :: c_os, c_ss, dg, diff, ds1, ds2, eps_I, &
                                                            eps_iter, g, omega_k, parts(4), s1, s2
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_t_type)                                 :: aj_Ib, aj_Ib_diff, aj_X, ja_Ik, &
                                                            ja_Ik_diff

      CALL timeset(routineN, handle)

      eps_iter = xas_tdp_control%gw2x_eps
      max_iter = xas_tdp_control%max_gw2x_iter
      c_os = xas_tdp_control%c_os
      c_ss = xas_tdp_control%c_ss
      batch_size = xas_tdp_control%batch_size

      ndo_mo = donor_state%ndo_mo
      output_unit = cp_logger_get_default_io_unit()

      nhomo = SIZE(homo_evals)
      nlumo = SIZE(lumo_evals)

      CALL get_qs_env(qs_env, para_env=para_env)

      !We use the Newton-Raphson method to find the zero of the function:
      !g(omega) = eps_I - omega + mp2 terms, dg(omega) = -1 + d/d_omega (mp2 terms)
      !We simply compute at each iteration: omega_k+1 = omega_k - g(omega_k)/dg(omega_k)

      !need transposed tensor of (ja|X) for optimal contraction scheme (s.t. (aj|X) block is on same
      !processor as (ja|X))
      CALL dbcsr_t_create(ja_X, aj_X)
      CALL dbcsr_t_copy(ja_X, aj_X, order=[2, 1, 3])

      !split the MO blocks into batches for memory friendly batched contraction
      !huge dense tensors never need to be stored
      CALL dbcsr_t_get_info(ja_X, nblks_total=nblks)
      ALLOCATE (mo_blk_size(nblks(1)))
      CALL dbcsr_t_get_info(ja_X, blk_size_1=mo_blk_size)

      tmp_sum = 0
      DO i = 1, nblks(1)
         tmp_sum = tmp_sum + mo_blk_size(i)
         IF (tmp_sum == nhomo) THEN
            nblk_occ = i
            nblk_virt = nblks(1) - i
            EXIT
         END IF
      END DO
      nbatch_occ = MAX(1, nblk_occ/batch_size)
      nbatch_virt = MAX(1, nblk_virt/batch_size)

      !Loop over donor_states
      DO ido_mo = 1, ndo_mo
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T5,A,I2,A,I4,A,/,T5,A)") &
               "- GW2X correction for donor MO with spin ", 1, &
               " and MO index ", donor_state%mo_indices(ido_mo, 1), ":", &
               "                         iteration                convergence (eV)"
            CALL m_flush(output_unit)
         END IF

         !starting values
         eps_I = homo_evals(first_domo + ido_mo - 1)
         omega_k = eps_I
         iloop = 0
         diff = 2.0_dp*eps_iter

         DO WHILE (ABS(diff) > eps_iter)
            iloop = iloop + 1

            !Compute the mp2 terms and their first derivative
            parts = 0.0_dp

            !We do batched contraction for (ja|Ik) and (ja|Ib) to never have to carry the full tensor
            DO ibatch = 1, nbatch_occ

               occ_bo = get_limit(nblk_occ, nbatch_occ, ibatch - 1)
               bounds_1d = [SUM(mo_blk_size(1:occ_bo(1) - 1)) + 1, SUM(mo_blk_size(1:occ_bo(2)))]

               CALL dbcsr_t_create(mo_template, ja_Ik)
               CALL dbcsr_t_contract(alpha=dbcsr_scalar(1.0_dp), tensor_1=ja_X, tensor_2=oI_Y(ido_mo), &
                                     beta=dbcsr_scalar(0.0_dp), tensor_3=ja_Ik, contract_1=[3], &
                                     notcontract_1=[1, 2], contract_2=[2], notcontract_2=[1], &
                                     map_1=[1, 2], map_2=[3], bounds_3=bounds_1d)

               !opposite-spin contribution
               CALL calc_os_oov_contrib(parts(1), parts(2), ja_Ik, homo_evals, lumo_evals, homo_evals, &
                                        omega_k, c_os, nhomo)

               bounds_2d(:, 2) = bounds_1d
               bounds_2d(1, 1) = nhomo + 1
               bounds_2d(2, 1) = nhomo + nlumo

               !same-spin contribution. Contraction only neede if c_ss != 0
               !directly compute the difference (ja|Ik) - (ka|Ij)
               IF (ABS(c_ss) > EPSILON(1.0_dp)) THEN

                  CALL dbcsr_t_create(ja_Ik, ja_Ik_diff, map1_2d=[1], map2_2d=[2, 3])
                  CALL dbcsr_t_copy(ja_Ik, ja_Ik_diff, move_data=.TRUE.)

                  CALL dbcsr_t_contract(alpha=dbcsr_scalar(-1.0_dp), tensor_1=oI_Y(ido_mo), tensor_2=aj_X, &
                                        beta=dbcsr_scalar(1.0_dp), tensor_3=ja_Ik_diff, contract_1=[2], &
                                        notcontract_1=[1], contract_2=[3], notcontract_2=[1, 2], &
                                        map_1=[1], map_2=[2, 3], bounds_2=[1, nhomo], bounds_3=bounds_2d)

                  CALL calc_ss_oov_contrib(parts(1), parts(2), ja_Ik_diff, homo_evals, lumo_evals, omega_k, c_ss)

                  CALL dbcsr_t_destroy(ja_Ik_diff)
               END IF !c_ss != 0

               CALL dbcsr_t_destroy(ja_Ik)
            END DO

            DO ibatch = 1, nbatch_virt

               virt_bo = get_limit(nblk_virt, nbatch_virt, ibatch - 1)
               bounds_1d = [SUM(mo_blk_size(1:nblk_occ + virt_bo(1) - 1)) + 1, &
                            SUM(mo_blk_size(1:nblk_occ + virt_bo(2)))]

               CALL dbcsr_t_create(mo_template, aj_Ib)
               CALL dbcsr_t_contract(alpha=dbcsr_scalar(1.0_dp), tensor_1=aj_X, tensor_2=oI_Y(ido_mo), &
                                     beta=dbcsr_scalar(0.0_dp), tensor_3=aj_Ib, contract_1=[3], &
                                     notcontract_1=[1, 2], contract_2=[2], notcontract_2=[1], &
                                     map_1=[1, 2], map_2=[3], bounds_3=bounds_1d)

               !opposite-spin contribution
               CALL calc_os_ovv_contrib(parts(3), parts(4), aj_Ib, lumo_evals, homo_evals, lumo_evals, &
                                        omega_k, c_os, nhomo, nhomo)

               !same-spin contribution, only if c_ss is not 0
               !directly compute the difference (aj|Ib) - (bj|Ia)
               IF (ABS(c_ss) > EPSILON(1.0_dp)) THEN
                  bounds_2d(1, 1) = 1
                  bounds_2d(2, 1) = nhomo
                  bounds_2d(:, 2) = bounds_1d

                  CALL dbcsr_t_create(aj_Ib, aj_Ib_diff, map1_2d=[1], map2_2d=[2, 3])
                  CALL dbcsr_t_copy(aj_Ib, aj_Ib_diff, move_data=.TRUE.)

                  CALL dbcsr_t_contract(alpha=dbcsr_scalar(-1.0_dp), tensor_1=oI_Y(ido_mo), tensor_2=ja_X, &
                                        beta=dbcsr_scalar(1.0_dp), tensor_3=aj_Ib_diff, contract_1=[2], &
                                        notcontract_1=[1], contract_2=[3], notcontract_2=[1, 2], &
                                        map_1=[1], map_2=[2, 3], &
                                        bounds_2=[nhomo + 1, nhomo + nlumo], bounds_3=bounds_2d)

                  CALL calc_ss_ovv_contrib(parts(3), parts(4), aj_Ib_diff, homo_evals, lumo_evals, omega_k, c_ss)

                  CALL dbcsr_t_destroy(aj_Ib_diff)
               END IF ! c_ss not 0

               CALL dbcsr_t_destroy(aj_Ib)
            END DO

            CALL mp_sum(parts, para_env%group)
            s1 = parts(1); ds1 = parts(2)
            s2 = parts(3); ds2 = parts(4)

            !evaluate g and its derivative
            g = eps_I - omega_k + s1 + s2
            dg = -1.0_dp + ds1 + ds2

            !compute the diff to the new step
            diff = -g/dg

            !and the new omega
            omega_k = omega_k + diff
            diff = diff*evolt

            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(T21,I18,F32.6)") &
                  iloop, diff
               CALL m_flush(output_unit)
            END IF

            IF (iloop > max_iter) THEN
               CPWARN("GW2X iteration not converged.")
               EXIT
            END IF
         END DO !while loop on eps_iter

         !compute the shift and update donor_state
         donor_state%gw2x_evals(ido_mo, 1) = omega_k

         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/T7,A,F11.6,/,T5,A,F11.6)") &
               "Final GW2X shift for this donor MO (eV):", &
               (donor_state%energy_evals(ido_mo, 1) - omega_k)*evolt
         END IF

      END DO !ido_mo

      CALL dbcsr_t_destroy(aj_X)

      CALL timestop(handle)

   END SUBROUTINE GW2X_rcs_iterations

! **************************************************************************************************
!> \brief Preforms the GW2X iterations in the open-shell shell formalism according to the
!>        Newton-Raphson method
!> \param first_domo index of the first core donor MO to consider, for each spin
!> \param ja_X semi-contracted tensors with j: occupied MO, a: virtual MO, X: RI basis element
!> \param oI_Y semi-contracted tensors with o: all MOs, I donor core MO, Y: RI basis element
!> \param mo_template tensor template for fully MO contracted tensor, for each spin combination
!> \param homo_evals ...
!> \param lumo_evals ...
!> \param donor_state ...
!> \param xas_tdp_control ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE GW2X_os_iterations(first_domo, ja_X, oI_Y, mo_template, homo_evals, lumo_evals, &
                                 donor_state, xas_tdp_control, qs_env)

      INTEGER, INTENT(IN)                                :: first_domo(2)
      TYPE(dbcsr_t_type), DIMENSION(:), INTENT(inout)    :: ja_X, oI_Y
      TYPE(dbcsr_t_type), DIMENSION(:, :), INTENT(inout) :: mo_template
      TYPE(cp_1d_r_p_type), DIMENSION(:), INTENT(in)     :: homo_evals, lumo_evals
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'GW2X_os_iterations'

      INTEGER :: batch_size, bounds_1d(2), bounds_2d(2, 2), handle, i, ibatch, ido_mo, iloop, &
         ispin, max_iter, nbatch_occ, nbatch_virt, nblk_occ, nblk_virt, nblks(3), ndo_mo, &
         nhomo(2), nlumo(2), nspins, occ_bo(2), other_spin, output_unit, tmp_sum, virt_bo(2)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: mo_blk_size
      REAL(dp)                                           :: c_os, c_ss, dg, diff, ds1, ds2, eps_I, &
                                                            eps_iter, g, omega_k, parts(4), s1, s2
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_t_type)                                 :: aj_Ib, aj_Ib_diff, aj_X(2), ja_Ik, &
                                                            ja_Ik_diff

      CALL timeset(routineN, handle)

      eps_iter = xas_tdp_control%gw2x_eps
      max_iter = xas_tdp_control%max_gw2x_iter
      c_os = xas_tdp_control%c_os
      c_ss = xas_tdp_control%c_ss
      batch_size = xas_tdp_control%batch_size

      nspins = 2
      ndo_mo = donor_state%ndo_mo
      output_unit = cp_logger_get_default_io_unit()

      DO ispin = 1, nspins
         nhomo(ispin) = SIZE(homo_evals(ispin)%array)
         nlumo(ispin) = SIZE(lumo_evals(ispin)%array)
      END DO

      CALL get_qs_env(qs_env, para_env=para_env)

      !We use the Newton-Raphson method to find the zero of the function:
      !g(omega) = eps_I - omega + mp2 terms, dg(omega) = -1 + d/d_omega (mp2 terms)
      !We simply compute at each iteration: omega_k+1 = omega_k - g(omega_k)/dg(omega_k)

      DO ispin = 1, nspins

         !need transposed tensor of (ja|X) for optimal contraction scheme,
         !s.t. (aj|X) block is on same processor as (ja|X)) and differences can be taken
         CALL dbcsr_t_create(ja_X(ispin), aj_X(ispin))
         CALL dbcsr_t_copy(ja_X(ispin), aj_X(ispin), order=[2, 1, 3])

      END DO ! ispin
      DO ispin = 1, nspins

         other_spin = 3 - ispin

         !split the MO blocks into batches for memory friendly batched contraction
         !huge dense tensors never need to be stored. Split MOs for the current spin
         CALL dbcsr_t_get_info(ja_X(ispin), nblks_total=nblks)
         ALLOCATE (mo_blk_size(nblks(1)))
         CALL dbcsr_t_get_info(ja_X(ispin), blk_size_1=mo_blk_size)

         tmp_sum = 0
         DO i = 1, nblks(1)
            tmp_sum = tmp_sum + mo_blk_size(i)
            IF (tmp_sum == nhomo(ispin)) THEN
               nblk_occ = i
               nblk_virt = nblks(1) - i
               EXIT
            END IF
         END DO
         nbatch_occ = MAX(1, nblk_occ/batch_size)
         nbatch_virt = MAX(1, nblk_virt/batch_size)

         !Loop over donor_states of the current spin
         DO ido_mo = 1, ndo_mo
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/,T5,A,I2,A,I4,A,/,T5,A)") &
                  "- GW2X correction for donor MO with spin ", ispin, &
                  " and MO index ", donor_state%mo_indices(ido_mo, ispin), ":", &
                  "                         iteration                convergence (eV)"
               CALL m_flush(output_unit)
            END IF

            !starting values
            eps_I = homo_evals(ispin)%array(first_domo(ispin) + ido_mo - 1)
            omega_k = eps_I
            iloop = 0
            diff = 2.0_dp*eps_iter

            DO WHILE (ABS(diff) > eps_iter)
               iloop = iloop + 1

               !Compute the mp2 terms and their first derivative
               parts = 0.0_dp

               !We do batched contraction for (ja|Ik) and (ja|Ib) to never have to carry the full tensor
               DO ibatch = 1, nbatch_occ

                  !opposite-spin contribution, i.e. (j_beta a_beta| I_alpha k_alpha) and vice-versa
                  !do the batching along k because same spin as donor MO
                  occ_bo = get_limit(nblk_occ, nbatch_occ, ibatch - 1)
                  bounds_1d = [SUM(mo_blk_size(1:occ_bo(1) - 1)) + 1, SUM(mo_blk_size(1:occ_bo(2)))]

                  CALL dbcsr_t_create(mo_template(other_spin, ispin), ja_Ik)
                  CALL dbcsr_t_contract(alpha=dbcsr_scalar(1.0_dp), tensor_1=ja_X(other_spin), &
                                        tensor_2=oI_Y((ispin - 1)*ndo_mo + ido_mo), &
                                        beta=dbcsr_scalar(0.0_dp), tensor_3=ja_Ik, contract_1=[3], &
                                        notcontract_1=[1, 2], contract_2=[2], notcontract_2=[1], &
                                        map_1=[1, 2], map_2=[3], bounds_3=bounds_1d)

                  CALL calc_os_oov_contrib(parts(1), parts(2), ja_Ik, homo_evals(other_spin)%array, &
                                           lumo_evals(other_spin)%array, homo_evals(ispin)%array, &
                                           omega_k, c_os, nhomo(other_spin))

                  CALL dbcsr_t_destroy(ja_Ik)

                  !same-spin contribution, need to compute (ja|Ik) - (ka|Ij), all with the current spin
                  !skip if c_ss == 0
                  IF (ABS(c_ss) > EPSILON(1.0_dp)) THEN

                     !same batching as opposite spin
                     CALL dbcsr_t_create(mo_template(ispin, ispin), ja_Ik)
                     CALL dbcsr_t_contract(alpha=dbcsr_scalar(1.0_dp), tensor_1=ja_X(ispin), &
                                           tensor_2=oI_Y((ispin - 1)*ndo_mo + ido_mo), &
                                           beta=dbcsr_scalar(0.0_dp), tensor_3=ja_Ik, contract_1=[3], &
                                           notcontract_1=[1, 2], contract_2=[2], notcontract_2=[1], &
                                           map_1=[1, 2], map_2=[3], bounds_3=bounds_1d)

                     bounds_2d(:, 2) = bounds_1d
                     bounds_2d(1, 1) = nhomo(ispin) + 1
                     bounds_2d(2, 1) = nhomo(ispin) + nlumo(ispin)

                     !the tensor difference is directly taken here
                     CALL dbcsr_t_create(ja_Ik, ja_Ik_diff, map1_2d=[1], map2_2d=[2, 3])
                     CALL dbcsr_t_copy(ja_Ik, ja_Ik_diff, move_data=.TRUE.)

                     CALL dbcsr_t_contract(alpha=dbcsr_scalar(-1.0_dp), tensor_1=oI_Y((ispin - 1)*ndo_mo + ido_mo), &
                                           tensor_2=aj_X(ispin), beta=dbcsr_scalar(1.0_dp), tensor_3=ja_Ik_diff, &
                                           contract_1=[2], notcontract_1=[1], contract_2=[3], notcontract_2=[1, 2], &
                                           map_1=[1], map_2=[2, 3], bounds_2=[1, nhomo(ispin)], bounds_3=bounds_2d)

                     CALL calc_ss_oov_contrib(parts(1), parts(2), ja_Ik_diff, homo_evals(ispin)%array, &
                                              lumo_evals(ispin)%array, omega_k, c_ss)

                     CALL dbcsr_t_destroy(ja_Ik_diff)
                     CALL dbcsr_t_destroy(ja_Ik)
                  END IF !c_ss !!= 0

               END DO

               DO ibatch = 1, nbatch_virt

                  !opposite-spin contribution, i.e. (a_beta j_beta| I_alpha b_alpha) and vice-versa
                  !do the batching along b because same spin as donor MO
                  virt_bo = get_limit(nblk_virt, nbatch_virt, ibatch - 1)
                  bounds_1d = [SUM(mo_blk_size(1:nblk_occ + virt_bo(1) - 1)) + 1, &
                               SUM(mo_blk_size(1:nblk_occ + virt_bo(2)))]

                  CALL dbcsr_t_create(mo_template(other_spin, ispin), aj_Ib)
                  CALL dbcsr_t_contract(alpha=dbcsr_scalar(1.0_dp), tensor_1=aj_X(other_spin), &
                                        tensor_2=oI_Y((ispin - 1)*ndo_mo + ido_mo), &
                                        beta=dbcsr_scalar(0.0_dp), tensor_3=aj_Ib, contract_1=[3], &
                                        notcontract_1=[1, 2], contract_2=[2], notcontract_2=[1], &
                                        map_1=[1, 2], map_2=[3], bounds_3=bounds_1d)

                  CALL calc_os_ovv_contrib(parts(3), parts(4), aj_Ib, lumo_evals(other_spin)%array, &
                                           homo_evals(other_spin)%array, lumo_evals(ispin)%array, &
                                           omega_k, c_os, nhomo(other_spin), nhomo(ispin))

                  CALL dbcsr_t_destroy(aj_Ib)

                  !same-spin contribution, need to compute (aj|Ib) - (bj|Ia), all with the current spin
                  !skip if c_ss == 0
                  IF (ABS(c_ss) > EPSILON(1.0_dp)) THEN

                     !same batching as opposite spin
                     CALL dbcsr_t_create(mo_template(ispin, ispin), aj_Ib)
                     CALL dbcsr_t_contract(alpha=dbcsr_scalar(1.0_dp), tensor_1=aj_X(ispin), &
                                           tensor_2=oI_Y((ispin - 1)*ndo_mo + ido_mo), &
                                           beta=dbcsr_scalar(0.0_dp), tensor_3=aj_Ib, contract_1=[3], &
                                           notcontract_1=[1, 2], contract_2=[2], notcontract_2=[1], &
                                           map_1=[1, 2], map_2=[3], bounds_3=bounds_1d)

                     bounds_2d(1, 1) = 1
                     bounds_2d(2, 1) = nhomo(ispin)
                     bounds_2d(:, 2) = bounds_1d

                     CALL dbcsr_t_create(aj_Ib, aj_Ib_diff, map1_2d=[1], map2_2d=[2, 3])
                     CALL dbcsr_t_copy(aj_Ib, aj_Ib_diff, move_data=.TRUE.)

                     CALL dbcsr_t_contract(alpha=dbcsr_scalar(-1.0_dp), tensor_1=oI_Y((ispin - 1)*ndo_mo + ido_mo), &
                                           tensor_2=ja_X(ispin), beta=dbcsr_scalar(1.0_dp), tensor_3=aj_Ib_diff, &
                                           contract_1=[2], notcontract_1=[1], contract_2=[3], &
                                           notcontract_2=[1, 2], map_1=[1], map_2=[2, 3], &
                                           bounds_2=[nhomo(ispin) + 1, nhomo(ispin) + nlumo(ispin)], &
                                           bounds_3=bounds_2d)

                     CALL calc_ss_ovv_contrib(parts(3), parts(4), aj_Ib_diff, homo_evals(ispin)%array, &
                                              lumo_evals(ispin)%array, omega_k, c_ss)

                     CALL dbcsr_t_destroy(aj_Ib_diff)
                     CALL dbcsr_t_destroy(aj_Ib)
                  END IF ! c_ss not 0

               END DO

               CALL mp_sum(parts, para_env%group)
               s1 = parts(1); ds1 = parts(2)
               s2 = parts(3); ds2 = parts(4)

               !evaluate g and its derivative
               g = eps_I - omega_k + s1 + s2
               dg = -1.0_dp + ds1 + ds2

               !compute the diff to the new step
               diff = -g/dg

               !and the new omega
               omega_k = omega_k + diff
               diff = diff*evolt

               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(T21,I18,F32.6)") &
                     iloop, diff
                  CALL m_flush(output_unit)
               END IF

               IF (iloop > max_iter) THEN
                  CPWARN("GW2X iteration not converged.")
                  EXIT
               END IF
            END DO !while loop on eps_iter

            !compute the shift and update donor_state
            donor_state%gw2x_evals(ido_mo, ispin) = omega_k

            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/T7,A,F11.6,/,T5,A,F11.6)") &
                  "Final GW2X shift for this donor MO (eV):", &
                  (donor_state%energy_evals(ido_mo, ispin) - omega_k)*evolt
            END IF

         END DO !ido_mo

         DEALLOCATE (mo_blk_size)
      END DO ! ispin

      DO ispin = 1, nspins
         CALL dbcsr_t_destroy(aj_X(ispin))
      END DO

      CALL timestop(handle)

   END SUBROUTINE GW2X_os_iterations

! **************************************************************************************************
!> \brief Takes the 3-center integrals from the ri_ex_3c tensor and returns a full tensor. Since
!>        ri_ex_3c is only half filled because of symmetry, we have to add the transpose
!>        and scale the diagonal blocks by 0.5
!> \param pq_X the full (desymmetrized) tensor containing the (pq|X) exchange integrals, in a new
!>        3d distribution and optimized block sizes
!> \param exat index of current excited atom
!> \param xas_tdp_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_full_pqX_from_3c_ex(pq_X, exat, xas_tdp_env, qs_env)

      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: pq_X
      INTEGER, INTENT(IN)                                :: exat
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'get_full_pqX_from_3c_ex', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, i, ind(3), natom, nblk_ri, nsgf_x
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: orb_blk_size, proc_dist_1, proc_dist_2, &
                                                            proc_dist_3
      INTEGER, DIMENSION(3)                              :: pdims
      LOGICAL                                            :: found
      REAL(dp)                                           :: factor
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: pblock
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_t_distribution_type)                    :: t_dist
      TYPE(dbcsr_t_iterator_type)                        :: iter
      TYPE(dbcsr_t_pgrid_type)                           :: t_pgrid
      TYPE(dbcsr_t_type)                                 :: pq_X_tmp, work

      NULLIFY (para_env)

      !create work tensor with same 2D dist as pq_X, but only keep excited atom along RI direction
      CALL get_qs_env(qs_env, para_env=para_env, natom=natom)
      CALL dbcsr_t_get_info(xas_tdp_env%ri_3c_ex, pdims=pdims)
      nsgf_x = SIZE(xas_tdp_env%ri_inv_ex, 1)
      nblk_ri = 1

      CALL dbcsr_t_pgrid_create(para_env%group, pdims, t_pgrid)
      ALLOCATE (proc_dist_1(natom), proc_dist_2(natom), orb_blk_size(natom))
      CALL dbcsr_t_get_info(xas_tdp_env%ri_3c_ex, proc_dist_1=proc_dist_1, proc_dist_2=proc_dist_2, &
                            blk_size_1=orb_blk_size)
      CALL dbcsr_t_distribution_new(t_dist, t_pgrid, nd_dist_1=proc_dist_1, nd_dist_2=proc_dist_2, &
                                    nd_dist_3=[(0, i=1, nblk_ri)])

      CALL dbcsr_t_create(work, name="(pq|X)", dist=t_dist, map1_2d=[1], map2_2d=[2, 3], &
                          blk_size_1=orb_blk_size, blk_size_2=orb_blk_size, blk_size_3=[nsgf_x])
      CALL dbcsr_t_distribution_destroy(t_dist)

      !dist of 3c_ex and work match, can simply copy blocks over. Diagonal with factor 0.5
      CALL dbcsr_t_iterator_start(iter, xas_tdp_env%ri_3c_ex)
      DO WHILE (dbcsr_t_iterator_blocks_left(iter))
         CALL dbcsr_t_iterator_next_block(iter, ind, blk)
         CALL dbcsr_t_get_block(xas_tdp_env%ri_3c_ex, ind, pblock, found)

         factor = 1.0_dp
         IF (ind(1) == ind(2)) factor = 0.5_dp
         IF (ind(3) /= exat) CYCLE

         CALL dbcsr_t_put_block(work, [ind(1), ind(2), 1], &
                                [orb_blk_size(ind(1)), orb_blk_size(ind(2)), nsgf_x], &
                                pblock, scale=factor)

         DEALLOCATE (pblock)
      END DO
      CALL dbcsr_t_iterator_stop(iter)
      CALL dbcsr_t_finalize(work)

      !create (pq|X) based on work and copy over
      CALL dbcsr_t_create(work, pq_X_tmp)
      CALL dbcsr_t_copy(work, pq_X_tmp)
      CALL dbcsr_t_copy(work, pq_X_tmp, order=[2, 1, 3], summation=.TRUE., move_data=.TRUE.)

      CALL dbcsr_t_destroy(work)

      !create the pgrid, based on the 2D dbcsr grid
      CALL dbcsr_t_pgrid_destroy(t_pgrid)
      pdims = 0
      CALL dbcsr_t_pgrid_create(para_env%group, pdims, t_pgrid, tensor_dims=[natom, natom, 1])

      !cyclic distribution accross all directions.
      ALLOCATE (proc_dist_3(nblk_ri))
      CALL dbcsr_t_default_distvec(natom, pdims(1), orb_blk_size, proc_dist_1)
      CALL dbcsr_t_default_distvec(natom, pdims(2), orb_blk_size, proc_dist_2)
      CALL dbcsr_t_default_distvec(nblk_ri, pdims(3), [nsgf_x], proc_dist_3)
      CALL dbcsr_t_distribution_new(t_dist, t_pgrid, nd_dist_1=proc_dist_1, nd_dist_2=proc_dist_2, &
                                    nd_dist_3=proc_dist_3)

      CALL dbcsr_t_create(pq_X, name="(pq|X)", dist=t_dist, map1_2d=[2, 3], map2_2d=[1], &
                          blk_size_1=orb_blk_size, blk_size_2=orb_blk_size, blk_size_3=[nsgf_x])
      CALL dbcsr_t_copy(pq_X_tmp, pq_X, move_data=.TRUE.)

      CALL dbcsr_t_distribution_destroy(t_dist)
      CALL dbcsr_t_pgrid_destroy(t_pgrid)
      CALL dbcsr_t_destroy(pq_X_tmp)

   END SUBROUTINE get_full_pqX_from_3c_ex

! **************************************************************************************************
!> \brief Contracts (pq|X) and (rI|Y) from AOs to MOs to (ja|X) and (oI|Y) respectively, where
!>        j is a occupied MO, a is a virtual MO and o is a general MO
!> \param ja_X partial contraction over occupied MOs j, virtual MOs a: (ja|X), for both spins (alpha-alpha or beta-beta)
!> \param oI_Y partial contraction over all MOs o and donor MOs I (can be more than 1 if 2p or open-shell)
!> \param ja_Io_template template to be able to build tensors after calling this routine, for each spin combination
!> \param mo_coeffs ...
!> \param nocc ...
!> \param nvirt ...
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note the multiplication by (X|Y)^-1 is included in the final (oI|Y) tensor. Only integrals with the
!>       same spin on one center are non-zero, i.e. (oI|Y) is non zero only if both o and Y have the same spin
! **************************************************************************************************
   SUBROUTINE contract_AOs_to_MOs(ja_X, oI_Y, ja_Io_template, mo_coeffs, nocc, nvirt, &
                                  donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: ja_X, oI_Y
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: ja_Io_template
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(INOUT)    :: mo_coeffs
      INTEGER, INTENT(IN)                                :: nocc(2), nvirt(2)
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'contract_AOs_to_MOs', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: bo(2), handle, i, ispin, jspin, &
                                                            nblk_aos, nblk_mos(2), nblk_occ(2), &
                                                            nblk_pqX(3), nblk_ri, nblk_virt(2), &
                                                            nspins
      INTEGER, DIMENSION(3)                              :: pdims
      INTEGER, DIMENSION(:), POINTER                     :: ao_blk_size, ao_col_dist, ao_row_dist, &
                                                            mo_dist_3, ri_blk_size, ri_dist_3
      INTEGER, DIMENSION(:, :), POINTER                  :: mat_pgrid
      TYPE(cp_1d_i_p_type), ALLOCATABLE, DIMENSION(:)    :: mo_blk_size, mo_col_dist, mo_row_dist
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type)                      :: mat_dist
      TYPE(dbcsr_distribution_type), POINTER             :: std_mat_dist
      TYPE(dbcsr_t_distribution_type)                    :: t_dist
      TYPE(dbcsr_t_pgrid_type)                           :: t_pgrid
      TYPE(dbcsr_t_type)                                 :: jq_X, pq_X, t_mo_coeffs
      TYPE(dbcsr_type)                                   :: dbcsr_mo_coeffs

      NULLIFY (ao_blk_size, ao_col_dist, ao_row_dist, mo_dist_3, ri_blk_size, ri_dist_3, mat_pgrid, &
               para_env, std_mat_dist)

      CALL timeset(routineN, handle)

      nspins = 1; IF (xas_tdp_control%do_uks .OR. xas_tdp_control%do_roks) nspins = 2

      !There are 2 contractions to do for the first tensor: (pq|X) --> (jq|X) --> (ja|X)
      !Because memory is the main concern, we move_data everytime at the cost of extra copies

      !Some quantities need to be stored for both spins, because they are later combined
      CALL get_qs_env(qs_env, para_env=para_env)
      ALLOCATE (mo_blk_size(nspins), mo_row_dist(nspins), mo_col_dist(nspins))
      ALLOCATE (ja_X(nspins))
      ALLOCATE (oI_Y(nspins*donor_state%ndo_mo))

      DO ispin = 1, nspins

         !First, we need a fully populated pq_X (spin-independent)
         CALL get_full_pqX_from_3c_ex(pq_X, donor_state%at_index, xas_tdp_env, qs_env)

         !Create the tensor pgrid. AOs and RI independent from spin
         IF (ispin == 1) THEN
            CALL dbcsr_t_get_info(pq_X, pdims=pdims, nblks_total=nblk_pqX)
            CALL dbcsr_t_pgrid_create(para_env%group, pdims, t_pgrid)
            nblk_aos = nblk_pqX(1)
            nblk_ri = nblk_pqX(3)
         END IF

         !Define MO block sizes, at worst, take one block per proc
         nblk_occ(ispin) = MAX(pdims(1), nocc(ispin)/16)
         nblk_virt(ispin) = MAX(pdims(2), nvirt(ispin)/16)
         nblk_mos(ispin) = nblk_occ(ispin) + nblk_virt(ispin)
         ALLOCATE (mo_blk_size(ispin)%array(nblk_mos(ispin)))
         DO i = 1, nblk_occ(ispin)
            bo = get_limit(nocc(ispin), nblk_occ(ispin), i - 1)
            mo_blk_size(ispin)%array(i) = bo(2) - bo(1) + 1
         END DO
         DO i = 1, nblk_virt(ispin)
            bo = get_limit(nvirt(ispin), nblk_virt(ispin), i - 1)
            mo_blk_size(ispin)%array(nblk_occ(ispin) + i) = bo(2) - bo(1) + 1
         END DO

         !Convert the fm mo_coeffs into a dbcsr matrix and then a tensor
         CALL get_qs_env(qs_env, dbcsr_dist=std_mat_dist)
         CALL dbcsr_distribution_get(std_mat_dist, pgrid=mat_pgrid)
         ALLOCATE (ao_blk_size(nblk_aos), ri_blk_size(nblk_ri))
         CALL dbcsr_t_get_info(pq_X, blk_size_1=ao_blk_size, blk_size_3=ri_blk_size)

         !we opt for a cyclic dist for the MOs (since they should be rather dense anyways)
         ALLOCATE (ao_row_dist(nblk_aos), mo_col_dist(ispin)%array(nblk_mos(ispin)))
         CALL dbcsr_t_default_distvec(nblk_aos, SIZE(mat_pgrid, 1), ao_blk_size, ao_row_dist)
         CALL dbcsr_t_default_distvec(nblk_mos(ispin), SIZE(mat_pgrid, 2), mo_blk_size(ispin)%array, &
                                      mo_col_dist(ispin)%array)
         CALL dbcsr_distribution_new(mat_dist, group=para_env%group, pgrid=mat_pgrid, &
                                     row_dist=ao_row_dist, col_dist=mo_col_dist(ispin)%array)

         CALL dbcsr_create(dbcsr_mo_coeffs, name="MO coeffs", matrix_type="N", dist=mat_dist, &
                           row_blk_size=ao_blk_size, col_blk_size=mo_blk_size(ispin)%array)
         CALL copy_fm_to_dbcsr(mo_coeffs(ispin)%matrix, dbcsr_mo_coeffs)

         CALL dbcsr_t_create(dbcsr_mo_coeffs, t_mo_coeffs)
         CALL dbcsr_t_copy_matrix_to_tensor(dbcsr_mo_coeffs, t_mo_coeffs)

         !prepare the (jq|X) tensor for the first contraction (over occupied MOs)
         ALLOCATE (mo_row_dist(ispin)%array(nblk_mos(ispin)), ao_col_dist(nblk_aos), ri_dist_3(nblk_ri))
         CALL dbcsr_t_default_distvec(nblk_mos(ispin), pdims(1), mo_blk_size(ispin)%array, mo_row_dist(ispin)%array)
         CALL dbcsr_t_default_distvec(nblk_aos, pdims(2), ao_blk_size, ao_col_dist)
         CALL dbcsr_t_default_distvec(nblk_ri, pdims(3), ri_blk_size, ri_dist_3)
         CALL dbcsr_t_distribution_new(t_dist, t_pgrid, nd_dist_1=mo_row_dist(ispin)%array, &
                                       nd_dist_2=ao_col_dist, nd_dist_3=ri_dist_3)

         CALL dbcsr_t_create(jq_X, name="(jq|X)", dist=t_dist, map1_2d=[1, 3], map2_2d=[2], &
                             blk_size_1=mo_blk_size(ispin)%array, blk_size_2=ao_blk_size, blk_size_3=ri_blk_size)
         CALL dbcsr_t_distribution_destroy(t_dist)

         !contract (pq|X) into (jq|X)
         CALL dbcsr_t_contract(alpha=dbcsr_scalar(1.0_dp), tensor_1=pq_X, tensor_2=t_mo_coeffs, &
                               beta=dbcsr_scalar(0.0_dp), tensor_3=jq_X, contract_1=[1], &
                               notcontract_1=[2, 3], contract_2=[1], notcontract_2=[2], &
                               map_1=[2, 3], map_2=[1], bounds_3=[1, nocc(ispin)], &!only want occupied MOs for j
                               move_data=.TRUE.)

         CALL dbcsr_t_destroy(pq_X)
         CALL dbcsr_t_copy_matrix_to_tensor(dbcsr_mo_coeffs, t_mo_coeffs)

         !prepare (ja|X) tensor for the second contraction (over virtual MOs)
         !only virtual-occupied bit of the first 2 indices is occupied + it should be dense
         !take blk dist such that blocks are evenly distributed
         CALL dbcsr_t_default_distvec(nblk_occ(ispin), pdims(1), mo_blk_size(ispin)%array(1:nblk_occ(ispin)), &
                                      mo_row_dist(ispin)%array(1:nblk_occ(ispin)))
         CALL dbcsr_t_default_distvec(nblk_virt(ispin), pdims(1), &
                                      mo_blk_size(ispin)%array(nblk_occ(ispin) + 1:nblk_mos(ispin)), &
                                      mo_row_dist(ispin)%array(nblk_occ(ispin) + 1:nblk_mos(ispin)))
         CALL dbcsr_t_default_distvec(nblk_occ(ispin), pdims(2), mo_blk_size(ispin)%array(1:nblk_occ(ispin)), &
                                      mo_col_dist(ispin)%array(1:nblk_occ(ispin)))
         CALL dbcsr_t_default_distvec(nblk_virt(ispin), pdims(2), &
                                      mo_blk_size(ispin)%array(nblk_occ(ispin) + 1:nblk_mos(ispin)), &
                                      mo_col_dist(ispin)%array(nblk_occ(ispin) + 1:nblk_mos(ispin)))
         CALL dbcsr_t_distribution_new(t_dist, t_pgrid, nd_dist_1=mo_row_dist(ispin)%array, &
                                       nd_dist_2=mo_col_dist(ispin)%array, nd_dist_3=ri_dist_3)

         CALL dbcsr_t_create(ja_X(ispin), name="(ja|X)", dist=t_dist, map1_2d=[1, 2], map2_2d=[3], &
                             blk_size_1=mo_blk_size(ispin)%array, blk_size_2=mo_blk_size(ispin)%array, &
                             blk_size_3=ri_blk_size)
         CALL dbcsr_t_distribution_destroy(t_dist)

         !contract (jq|X) into (ja|X)
         CALL dbcsr_t_contract(alpha=dbcsr_scalar(1.0_dp), tensor_1=jq_X, tensor_2=t_mo_coeffs, &
                               beta=dbcsr_scalar(0.0_dp), tensor_3=ja_X(ispin), contract_1=[2], &
                               notcontract_1=[1, 3], contract_2=[1], notcontract_2=[2], &
                               map_1=[1, 3], map_2=[2], move_data=.TRUE., &
                               bounds_3=[nocc(ispin) + 1, nocc(ispin) + nvirt(ispin)])

         CALL dbcsr_t_destroy(jq_X)
         CALL dbcsr_t_copy_matrix_to_tensor(dbcsr_mo_coeffs, t_mo_coeffs)

         !Finally, get the oI_Y tensors
         CALL get_oIY_tensors(oI_Y, ispin, ao_blk_size, mo_blk_size(ispin)%array, ri_blk_size, &
                              t_mo_coeffs, donor_state, xas_tdp_env, xas_tdp_control, qs_env)

         !intermediate clen-up
         CALL dbcsr_t_destroy(t_mo_coeffs)
         CALL dbcsr_distribution_release(mat_dist)
         CALL dbcsr_release(dbcsr_mo_coeffs)
         DEALLOCATE (ao_col_dist, ri_dist_3, ri_blk_size, ao_blk_size, ao_row_dist)

      END DO !ispin

      !create a empty tensor template for the fully contracted (ja|Io) MO integrals, for all spin
      !configureations: alpha-alpha|alpha-alpha, alpha-alpha|beta-beta, etc.
      ALLOCATE (ja_Io_template(nspins, nspins))
      DO ispin = 1, nspins
         DO jspin = 1, nspins
            ALLOCATE (mo_dist_3(nblk_mos(jspin)))
            CALL dbcsr_t_default_distvec(nblk_occ(jspin), pdims(3), mo_blk_size(jspin)%array(1:nblk_occ(jspin)), &
                                         mo_dist_3(1:nblk_occ(jspin)))
            CALL dbcsr_t_default_distvec(nblk_virt(jspin), pdims(3), &
                                         mo_blk_size(jspin)%array(nblk_occ(jspin) + 1:nblk_mos(jspin)), &
                                         mo_dist_3(nblk_occ(jspin) + 1:nblk_mos(jspin)))
            CALL dbcsr_t_distribution_new(t_dist, t_pgrid, nd_dist_1=mo_row_dist(ispin)%array, &
                                          nd_dist_2=mo_col_dist(ispin)%array, nd_dist_3=mo_dist_3)

            CALL dbcsr_t_create(ja_Io_template(ispin, jspin), name="(ja|Io)", dist=t_dist, map1_2d=[1, 2], &
                                map2_2d=[3], blk_size_1=mo_blk_size(ispin)%array, &
                                blk_size_2=mo_blk_size(ispin)%array, blk_size_3=mo_blk_size(jspin)%array)
            CALL dbcsr_t_distribution_destroy(t_dist)
            DEALLOCATE (mo_dist_3)
         END DO
      END DO

      !clean-up
      CALL dbcsr_t_pgrid_destroy(t_pgrid)
      DO ispin = 1, nspins
         DEALLOCATE (mo_blk_size(ispin)%array)
         DEALLOCATE (mo_col_dist(ispin)%array)
         DEALLOCATE (mo_row_dist(ispin)%array)
      END DO

      CALL timestop(handle)

   END SUBROUTINE contract_AOs_to_MOs

! **************************************************************************************************
!> \brief Contracts the (oI|Y) tensors, for each donor MO
!> \param oI_Y the contracted tensr. It is assumed to be allocated outside of this routine
!> \param ispin ...
!> \param ao_blk_size ...
!> \param mo_blk_size ...
!> \param ri_blk_size ...
!> \param t_mo_coeffs ...
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_oIY_tensors(oI_Y, ispin, ao_blk_size, mo_blk_size, ri_blk_size, t_mo_coeffs, &
                              donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: oI_Y
      INTEGER, INTENT(IN)                                :: ispin
      INTEGER, DIMENSION(:), POINTER                     :: ao_blk_size, mo_blk_size, ri_blk_size
      TYPE(dbcsr_t_type), INTENT(inout)                  :: t_mo_coeffs
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_oIY_tensors'

      INTEGER                                            :: blk, bo(2), handle, i, ido_mo, ind(2), &
                                                            natom, nblk_aos, nblk_mos, nblk_ri, &
                                                            ndo_mo, pdims_2d(2), proc_id
      INTEGER, DIMENSION(:), POINTER                     :: ao_row_dist, mo_row_dist, ri_col_dist
      INTEGER, DIMENSION(:, :), POINTER                  :: mat_pgrid
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pblock
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type), POINTER             :: std_mat_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: pI_Y
      TYPE(dbcsr_t_distribution_type)                    :: t_dist
      TYPE(dbcsr_t_iterator_type)                        :: iter
      TYPE(dbcsr_t_pgrid_type)                           :: t_pgrid
      TYPE(dbcsr_t_type)                                 :: t_pI_Y, t_work

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, natom=natom, para_env=para_env, dbcsr_dist=std_mat_dist)
      ndo_mo = donor_state%ndo_mo
      nblk_aos = SIZE(ao_blk_size)
      nblk_mos = SIZE(mo_blk_size)
      nblk_ri = SIZE(ri_blk_size)

      !We first contract (pq|X) over q into I using kernel routines (goes over all MOs and spins)
      CALL contract2_AO_to_doMO(pI_Y, "EXCHANGE", donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      !multiply by (X|Y)^-1
      CALL ri_all_blocks_mm(pI_Y, xas_tdp_env%ri_inv_ex)

      !get standaed 2d matrix proc grid
      CALL dbcsr_distribution_get(std_mat_dist, pgrid=mat_pgrid)

      !Loop over donor MOs of this spin
      DO ido_mo = (ispin - 1)*ndo_mo + 1, ispin*ndo_mo

         !cast the matrix into a tensor
         CALL dbcsr_t_create(pI_Y(ido_mo)%matrix, t_work)
         CALL dbcsr_t_copy_matrix_to_tensor(pI_Y(ido_mo)%matrix, t_work)

         !find col proc_id of the only populated column of t_work
         ALLOCATE (ri_col_dist(natom))
         CALL dbcsr_t_get_info(t_work, proc_dist_2=ri_col_dist)
         proc_id = ri_col_dist(donor_state%at_index)
         DEALLOCATE (ri_col_dist)

         !preapre (oI_Y) tensor and (pI|Y) tensor in proper dist and blk sizes
         pdims_2d(1) = SIZE(mat_pgrid, 1); pdims_2d(2) = SIZE(mat_pgrid, 2)
         CALL dbcsr_t_pgrid_create(para_env%group, pdims_2d, t_pgrid)

         ALLOCATE (ri_col_dist(nblk_ri), ao_row_dist(nblk_aos), mo_row_dist(nblk_mos))
         CALL dbcsr_t_get_info(t_work, proc_dist_1=ao_row_dist)
         ri_col_dist = proc_id

         CALL dbcsr_t_distribution_new(t_dist, t_pgrid, nd_dist_1=ao_row_dist, nd_dist_2=ri_col_dist)
         CALL dbcsr_t_create(t_pI_Y, name="(pI|Y)", dist=t_dist, map1_2d=[1], map2_2d=[2], &
                             blk_size_1=ao_blk_size, blk_size_2=ri_blk_size)
         CALL dbcsr_t_distribution_destroy(t_dist)

         !copy block by block, dist match
         CALL dbcsr_t_iterator_start(iter, t_work)
         DO WHILE (dbcsr_t_iterator_blocks_left(iter))
            CALL dbcsr_t_iterator_next_block(iter, ind, blk)
            CALL dbcsr_t_get_block(t_work, ind, pblock, found)

            DO i = 1, nblk_ri
               bo(1) = SUM(ri_blk_size(1:i - 1)) + 1
               bo(2) = bo(1) + ri_blk_size(i) - 1
               CALL dbcsr_t_put_block(t_pI_Y, [ind(1), i], [ao_blk_size(ind(1)), ri_blk_size(i)], &
                                      pblock(:, bo(1):bo(2)))
            END DO

            DEALLOCATE (pblock)
         END DO
         CALL dbcsr_t_iterator_stop(iter)
         CALL dbcsr_t_finalize(t_pI_Y)

         !get optimal pgrid  for (oI|Y)
         CALL dbcsr_t_pgrid_destroy(t_pgrid)
         pdims_2d = 0
         CALL dbcsr_t_pgrid_create(para_env%group, pdims_2d, t_pgrid, tensor_dims=[nblk_mos, nblk_ri])

         CALL dbcsr_t_default_distvec(nblk_aos, pdims_2d(1), ao_blk_size, ao_row_dist)
         CALL dbcsr_t_default_distvec(nblk_mos, pdims_2d(1), mo_blk_size, mo_row_dist)
         CALL dbcsr_t_default_distvec(nblk_ri, pdims_2d(2), ri_blk_size, ri_col_dist)

         !transfer pI_Y to the correct pgrid
         CALL dbcsr_t_destroy(t_work)
         CALL dbcsr_t_distribution_new(t_dist, t_pgrid, nd_dist_1=ao_row_dist, nd_dist_2=ri_col_dist)
         CALL dbcsr_t_create(t_work, name="t_pI_Y", dist=t_dist, map1_2d=[1], map2_2d=[2], &
                             blk_size_1=ao_blk_size, blk_size_2=ri_blk_size)
         CALL dbcsr_t_copy(t_pI_Y, t_work, move_data=.TRUE.)
         CALL dbcsr_t_distribution_destroy(t_dist)

         !create (oI|Y)
         CALL dbcsr_t_distribution_new(t_dist, t_pgrid, nd_dist_1=mo_row_dist, nd_dist_2=ri_col_dist)
         CALL dbcsr_t_create(oI_Y(ido_mo), name="(oI|Y)", dist=t_dist, map1_2d=[1], map2_2d=[2], &
                             blk_size_1=mo_blk_size, blk_size_2=ri_blk_size)
         CALL dbcsr_t_distribution_destroy(t_dist)

         !contract (pI|Y) into (oI|Y)
         CALL dbcsr_t_contract(alpha=dbcsr_scalar(1.0_dp), tensor_1=t_work, tensor_2=t_mo_coeffs, &
                               beta=dbcsr_scalar(0.0_dp), tensor_3=oI_Y(ido_mo), contract_1=[1], &
                               notcontract_1=[2], contract_2=[1], notcontract_2=[2], &
                               map_1=[2], map_2=[1]) !no bound, all MOs needed

         !intermediate clean-up
         CALL dbcsr_t_destroy(t_work)
         CALL dbcsr_t_destroy(t_pI_Y)
         CALL dbcsr_t_pgrid_destroy(t_pgrid)
         DEALLOCATE (ri_col_dist, ao_row_dist, mo_row_dist)

      END DO !ido_mo

      !final clean-up
      CALL dbcsr_deallocate_matrix_set(pI_Y)

      CALL timestop(handle)

   END SUBROUTINE get_oIY_tensors

! **************************************************************************************************
!> \brief Computes the same spin, occupied-occupied-virtual MO contribution to the electron propagator:
!>        0.5 * sum_ajk |<Ia||jk>|^2/(omega + eps_a - epsj_j - eps_k) and its 1st derivative wrt omega:
!>        -0.5 * sum_ajk |<Ia||jk>|^2/(omega + eps_a - epsj_j - eps_k)**2
!> \param contrib ...
!> \param dev the first derivative
!> \param ja_Ik_diff ... contains the (ja|Ik) - (ka|Ij) tensor
!> \param occ_evals ...
!> \param virt_evals ...
!> \param omega ...
!> \param c_ss ...
!> \note since the is same-spin, there is only one possibility for occ_evals and virt_evals
! **************************************************************************************************
   SUBROUTINE calc_ss_oov_contrib(contrib, dev, ja_Ik_diff, occ_evals, virt_evals, omega, c_ss)

      REAL(dp), INTENT(inout)                            :: contrib, dev
      TYPE(dbcsr_t_type), INTENT(inout)                  :: ja_Ik_diff
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: occ_evals, virt_evals
      REAL(dp), INTENT(in)                               :: omega, c_ss

      CHARACTER(len=*), PARAMETER :: routineN = 'calc_ss_oov_contrib', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: a, blk, boff(3), bsize(3), handle, idx1, &
                                                            idx2, idx3, ind(3), j, k, nocc
      LOGICAL                                            :: found
      REAL(dp)                                           :: denom, tmp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: tensor_blk
      TYPE(dbcsr_t_iterator_type)                        :: iter

      CALL timeset(routineN, handle)

      !<Ia||jk> = <Ia|jk> - <Ia|kj> = (Ij|ak) - (Ik|aj) = (ka|Ij) - (ja|Ik)
      !Note: the same spin contribution only involve spib-orbitals that are all of the same spin

      nocc = SIZE(occ_evals, 1)

      !Iterate over the tensors and sum. Both tensors have same dist
      CALL dbcsr_t_iterator_start(iter, ja_Ik_diff)
      DO WHILE (dbcsr_t_iterator_blocks_left(iter))
         CALL dbcsr_t_iterator_next_block(iter, ind, blk, blk_offset=boff, blk_size=bsize)
         CALL dbcsr_t_get_block(ja_Ik_diff, ind, tensor_blk, found)

         IF (found) THEN

            DO idx3 = 1, bsize(3)
               DO idx2 = 1, bsize(2)
                  DO idx1 = 1, bsize(1)

                     !get proper MO indices
                     j = boff(1) + idx1 - 1
                     a = boff(2) + idx2 - 1 - nocc
                     k = boff(3) + idx3 - 1

                     !the denominator
                     denom = omega + virt_evals(a) - occ_evals(j) - occ_evals(k)

                     !the same spin contribution
                     tmp = c_ss*tensor_blk(idx1, idx2, idx3)**2

                     contrib = contrib + 0.5_dp*tmp/denom
                     dev = dev - 0.5_dp*tmp/denom**2

                  END DO
               END DO
            END DO
         END IF
         DEALLOCATE (tensor_blk)
      END DO
      CALL dbcsr_t_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE calc_ss_oov_contrib

! **************************************************************************************************
!> \brief Computes the opposite spin, occupied-occupied-virtual MO contribution to the electron propagator:
!>        0.5 * sum_ajk |<Ia||jk>|^2/(omega + eps_a - epsj_j - eps_k) and its 1st derivative wrt omega:
!>        -0.5 * sum_ajk |<Ia||jk>|^2/(omega + eps_a - epsj_j - eps_k)**2
!> \param contrib ...
!> \param dev the first derivative
!> \param ja_Ik ...
!> \param j_evals ocucpied evals for j MO
!> \param a_evals virtual evals for a MO
!> \param k_evals ocucpied evals for k MO
!> \param omega ...
!> \param c_os ...
!> \param a_offset the number of occupied MOs for the same spin as a MOs
!> \note since this is opposite-spin, evals might be different for different spins
! **************************************************************************************************
   SUBROUTINE calc_os_oov_contrib(contrib, dev, ja_Ik, j_evals, a_evals, k_evals, omega, c_os, a_offset)

      REAL(dp), INTENT(inout)                            :: contrib, dev
      TYPE(dbcsr_t_type), INTENT(inout)                  :: ja_Ik
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: j_evals, a_evals, k_evals
      REAL(dp), INTENT(in)                               :: omega, c_os
      INTEGER, INTENT(IN)                                :: a_offset

      CHARACTER(len=*), PARAMETER :: routineN = 'calc_os_oov_contrib', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: a, blk, boff(3), bsize(3), handle, idx1, &
                                                            idx2, idx3, ind(3), j, k
      LOGICAL                                            :: found
      REAL(dp)                                           :: denom, tmp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: ja_Ik_blk
      TYPE(dbcsr_t_iterator_type)                        :: iter

      CALL timeset(routineN, handle)

      !<Ia||jk> = <Ia|jk> - <Ia|kj> = (Ij|ak) - (Ik|aj) = (ka|Ij) - (ja|Ik)
      !Note: the opposite spin contribution comes in 2 parts, once (ka|Ij) and once (ja|Ik) only,
      !      where both spin-orbitals on one center have the same spin, but it is the opposite of
      !      the spin on the other center. Because it is eventually summed, can consider only one
      !      of the 2 terms, but with a factor 2

      !Iterate over the tensor and sum
      CALL dbcsr_t_iterator_start(iter, ja_Ik)
      DO WHILE (dbcsr_t_iterator_blocks_left(iter))
         CALL dbcsr_t_iterator_next_block(iter, ind, blk, blk_offset=boff, blk_size=bsize)
         CALL dbcsr_t_get_block(ja_Ik, ind, ja_Ik_blk, found)

         IF (found) THEN

            DO idx3 = 1, bsize(3)
               DO idx2 = 1, bsize(2)
                  DO idx1 = 1, bsize(1)

                     !get proper MO indices
                     j = boff(1) + idx1 - 1
                     a = boff(2) + idx2 - 1 - a_offset
                     k = boff(3) + idx3 - 1

                     !the denominator
                     denom = omega + a_evals(a) - j_evals(j) - k_evals(k)

                     !the opposite spin contribution
                     tmp = c_os*ja_Ik_blk(idx1, idx2, idx3)**2

                     !take factor 2 into acocunt (2 x 0.5 = 1)
                     contrib = contrib + tmp/denom
                     dev = dev - tmp/denom**2

                  END DO
               END DO
            END DO
         END IF
         DEALLOCATE (ja_Ik_blk)
      END DO
      CALL dbcsr_t_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE calc_os_oov_contrib

! **************************************************************************************************
!> \brief Computes the same-spin occupied-virtual-virtual MO contribution to the electron propagator:
!>        0.5 * sum_abj |<Ij||ab>|^2/(omega + eps_j - eps_a - eps_b) as well as its first derivative:
!>        -0.5 * sum_abj |<Ij||ab>|^2/(omega + eps_j - eps_a - eps_b)**2
!> \param contrib ...
!> \param dev the first derivative
!> \param aj_Ib_diff contatins the (aj|Ib) - (bj|Ia) tensor
!> \param occ_evals ...
!> \param virt_evals ...
!> \param omega ...
!> \param c_ss ...
!> \note since the is same-spin, there is only one possibility for occ_evals and virt_evals
! **************************************************************************************************
   SUBROUTINE calc_ss_ovv_contrib(contrib, dev, aj_Ib_diff, occ_evals, virt_evals, omega, c_ss)

      REAL(dp), INTENT(inout)                            :: contrib, dev
      TYPE(dbcsr_t_type), INTENT(inout)                  :: aj_Ib_diff
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: occ_evals, virt_evals
      REAL(dp), INTENT(in)                               :: omega, c_ss

      CHARACTER(len=*), PARAMETER :: routineN = 'calc_ss_ovv_contrib', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: a, b, blk, boff(3), bsize(3), handle, &
                                                            idx1, idx2, idx3, ind(3), j, nocc
      LOGICAL                                            :: found
      REAL(dp)                                           :: denom, tmp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: tensor_blk
      TYPE(dbcsr_t_iterator_type)                        :: iter

      CALL timeset(routineN, handle)

      !<Ij||ab> = <Ij|ab> - <Ij|ba> = (Ia|jb) - (Ib|ja) = (jb|Ia) - (ja|Ib)
      !Notes: only non-zero contribution if all MOs have the same spin

      nocc = SIZE(occ_evals, 1)

      !tensors have matching distributions, can do that safely
      CALL dbcsr_t_iterator_start(iter, aj_Ib_diff)
      DO WHILE (dbcsr_t_iterator_blocks_left(iter))
         CALL dbcsr_t_iterator_next_block(iter, ind, blk, blk_offset=boff, blk_size=bsize)
         CALL dbcsr_t_get_block(aj_Ib_diff, ind, tensor_blk, found)

         IF (found) THEN

            DO idx3 = 1, bsize(3)
               DO idx2 = 1, bsize(2)
                  DO idx1 = 1, bsize(1)

                     !get proper MO indices
                     a = boff(1) + idx1 - 1 - nocc
                     j = boff(2) + idx2 - 1
                     b = boff(3) + idx3 - 1 - nocc

                     !the common denominator
                     denom = omega + occ_evals(j) - virt_evals(a) - virt_evals(b)

                     !the same spin contribution
                     tmp = c_ss*tensor_blk(idx1, idx2, idx3)**2

                     contrib = contrib + 0.5_dp*tmp/denom
                     dev = dev - 0.5_dp*tmp/denom**2

                  END DO
               END DO
            END DO
         END IF
         DEALLOCATE (tensor_blk)
      END DO
      CALL dbcsr_t_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE calc_ss_ovv_contrib

! **************************************************************************************************
!> \brief Computes the opposite-spin occupied-virtual-virtual MO contribution to the electron propagator:
!>        0.5 * sum_abj |<Ij||ab>|^2/(omega + eps_j - eps_a - eps_b) as well as its first derivative:
!>        -0.5 * sum_abj |<Ij||ab>|^2/(omega + eps_j - eps_a - eps_b)**2
!> \param contrib ...
!> \param dev the first derivative
!> \param aj_Ib ...
!> \param a_evals virtual evals for a MO
!> \param j_evals occupied evals for j MO
!> \param b_evals virtual evals for b MO
!> \param omega ...
!> \param c_os ...
!> \param a_offset number of occupied MOs for the same spin as a MO
!> \param b_offset number of occupied MOs for the same spin as b MO
!> \note since this is opposite-spin, evals might be different for different spins
! **************************************************************************************************
   SUBROUTINE calc_os_ovv_contrib(contrib, dev, aj_Ib, a_evals, j_evals, b_evals, omega, c_os, &
                                  a_offset, b_offset)

      REAL(dp), INTENT(inout)                            :: contrib, dev
      TYPE(dbcsr_t_type), INTENT(inout)                  :: aj_Ib
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: a_evals, j_evals, b_evals
      REAL(dp), INTENT(in)                               :: omega, c_os
      INTEGER, INTENT(IN)                                :: a_offset, b_offset

      CHARACTER(len=*), PARAMETER :: routineN = 'calc_os_ovv_contrib', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: a, b, blk, boff(3), bsize(3), handle, &
                                                            idx1, idx2, idx3, ind(3), j
      LOGICAL                                            :: found
      REAL(dp)                                           :: denom, tmp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: aj_Ib_blk
      TYPE(dbcsr_t_iterator_type)                        :: iter

      CALL timeset(routineN, handle)

      !<Ij||ab> = <Ij|ab> - <Ij|ba> = (Ia|jb) - (Ib|ja) = (jb|Ia) - (ja|Ib)
      !Notes: only 2 distinct contributions, once from (jb|Ia) and once form (ja|Ib) only, when the 2
      !       MOs on one center have one spin and the 2 MOs on the other center have another spin
      !       In the end, the sum is such that can take one of those with a factor 2

      CALL dbcsr_t_iterator_start(iter, aj_Ib)
      DO WHILE (dbcsr_t_iterator_blocks_left(iter))
         CALL dbcsr_t_iterator_next_block(iter, ind, blk, blk_offset=boff, blk_size=bsize)
         CALL dbcsr_t_get_block(aj_Ib, ind, aj_Ib_blk, found)

         IF (found) THEN

            DO idx3 = 1, bsize(3)
               DO idx2 = 1, bsize(2)
                  DO idx1 = 1, bsize(1)

                     !get proper MO indices
                     a = boff(1) + idx1 - 1 - a_offset
                     j = boff(2) + idx2 - 1
                     b = boff(3) + idx3 - 1 - b_offset

                     !the denominator
                     denom = omega + j_evals(j) - a_evals(a) - b_evals(b)

                     !the opposite-spin contribution. Factor 2 taken into account (2 x 0.5 = 1)
                     tmp = c_os*(aj_Ib_blk(idx1, idx2, idx3))**2

                     contrib = contrib + tmp/denom
                     dev = dev - tmp/denom**2

                  END DO
               END DO
            END DO
         END IF
         DEALLOCATE (aj_Ib_blk)
      END DO
      CALL dbcsr_t_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE calc_os_ovv_contrib

! **************************************************************************************************
!> \brief We try to compute the spin-orbit splitting via perturbation theory. We keep it
!>\        cheap by only inculding the degenerate states (2p, 3d, 3p, etc.).
!> \param soc_shifts the SOC corrected orbital shifts to apply to original energies, for both spins
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_soc_splitting(soc_shifts, donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(out)                                     :: soc_shifts
      TYPE(donor_state_type), POINTER                    :: donor_state
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_soc_splitting'

      COMPLEX(dp), ALLOCATABLE, DIMENSION(:, :)          :: evecs, hami
      INTEGER                                            :: beta_spin, handle, ialpha, ibeta, &
                                                            ido_mo, ispin, nao, ndo_mo, ndo_so, &
                                                            nspins
      REAL(dp)                                           :: alpha_tot_contrib, beta_tot_contrib
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: evals
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: tmp_shifts
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type), POINTER                         :: hami_cfm
      TYPE(cp_fm_struct_type), POINTER                   :: ao_domo_struct, domo_domo_struct, &
                                                            doso_doso_struct
      TYPE(cp_fm_type), POINTER                          :: alpha_gs_coeffs, ao_domo_work, &
                                                            beta_gs_coeffs, domo_domo_work, &
                                                            img_fm, real_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(dbcsr_type), POINTER                          :: orb_soc_x, orb_soc_y, orb_soc_z

      NULLIFY (matrix_ks, para_env, blacs_env, hami_cfm, ao_domo_struct, domo_domo_struct, &
               doso_doso_struct, alpha_gs_coeffs, ao_domo_work, beta_gs_coeffs, domo_domo_work, &
               img_fm, real_fm, orb_soc_x, orb_soc_y, orb_soc_z)

      CALL timeset(routineN, handle)

      ! Idea: we compute the SOC matrix in the space of the degenerate spin-orbitals, add it to
      !       the KS matrix in the same basis, diagonalize the whole thing and get the corrected energies
      !       for SOC

      CALL get_qs_env(qs_env, matrix_ks=matrix_ks, para_env=para_env, blacs_env=blacs_env)

      orb_soc_x => xas_tdp_env%orb_soc(1)%matrix
      orb_soc_y => xas_tdp_env%orb_soc(2)%matrix
      orb_soc_z => xas_tdp_env%orb_soc(3)%matrix

      ! Whether it is open-shell or not, we have 2*ndo_mo spin-orbitals
      nspins = 2
      ndo_mo = donor_state%ndo_mo
      ndo_so = nspins*ndo_mo
      CALL dbcsr_get_info(matrix_ks(1)%matrix, nfullrows_total=nao)

      ! Build the fm infrastructure
      CALL cp_fm_struct_create(ao_domo_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=nao, ncol_global=ndo_mo)
      CALL cp_fm_struct_create(domo_domo_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ndo_mo, ncol_global=ndo_mo)
      CALL cp_fm_struct_create(doso_doso_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=ndo_so, ncol_global=ndo_so)

      CALL cp_fm_create(alpha_gs_coeffs, ao_domo_struct)
      CALL cp_fm_create(beta_gs_coeffs, ao_domo_struct)
      CALL cp_fm_create(ao_domo_work, ao_domo_struct)
      CALL cp_fm_create(domo_domo_work, domo_domo_struct)
      CALL cp_fm_create(real_fm, doso_doso_struct)
      CALL cp_fm_create(img_fm, doso_doso_struct)

      ! Put the gs_coeffs in the correct format.
      IF (xas_tdp_control%do_uks) THEN

         CALL cp_fm_to_fm_submat(msource=donor_state%gs_coeffs, mtarget=alpha_gs_coeffs, nrow=nao, &
                                 ncol=ndo_mo, s_firstrow=1, s_firstcol=1, t_firstrow=1, t_firstcol=1)
         CALL cp_fm_to_fm_submat(msource=donor_state%gs_coeffs, mtarget=beta_gs_coeffs, nrow=nao, &
                                 ncol=ndo_mo, s_firstrow=1, s_firstcol=ndo_mo + 1, t_firstrow=1, t_firstcol=1)

      ELSE

         CALL cp_fm_to_fm(donor_state%gs_coeffs, alpha_gs_coeffs)
         CALL cp_fm_to_fm(donor_state%gs_coeffs, beta_gs_coeffs)
      END IF

      ! Compute the KS matrix in this basis, add it to the real part of the final matrix
      !alpha-alpha block in upper left quadrant
      CALL cp_dbcsr_sm_fm_multiply(matrix_ks(1)%matrix, alpha_gs_coeffs, ao_domo_work, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, alpha_gs_coeffs, ao_domo_work, 0.0_dp, &
                   domo_domo_work)
      CALL cp_fm_to_fm_submat(msource=domo_domo_work, mtarget=real_fm, nrow=ndo_mo, ncol=ndo_mo, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=1, t_firstcol=1)

      !beta-beta block in lower right quadrant
      beta_spin = 1; IF (xas_tdp_control%do_uks .OR. xas_tdp_control%do_roks) beta_spin = 2
      CALL cp_dbcsr_sm_fm_multiply(matrix_ks(beta_spin)%matrix, beta_gs_coeffs, ao_domo_work, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, beta_gs_coeffs, ao_domo_work, 0.0_dp, &
                   domo_domo_work)
      CALL cp_fm_to_fm_submat(msource=domo_domo_work, mtarget=real_fm, nrow=ndo_mo, ncol=ndo_mo, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=ndo_mo + 1, t_firstcol=ndo_mo + 1)

      ! Compute the SOC matrix elements and add them to the real or imaginary part of the matrix
      ! alpha-alpha block, only Hz not zero, purely imaginary, addition
      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, alpha_gs_coeffs, ao_domo_work, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, alpha_gs_coeffs, ao_domo_work, 0.0_dp, &
                   domo_domo_work)
      CALL cp_fm_to_fm_submat(msource=domo_domo_work, mtarget=img_fm, nrow=ndo_mo, ncol=ndo_mo, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=1, t_firstcol=1)

      ! beta-beta block, only Hz not zero, purely imaginary, substraciton
      CALL cp_dbcsr_sm_fm_multiply(orb_soc_z, beta_gs_coeffs, ao_domo_work, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, -1.0_dp, beta_gs_coeffs, ao_domo_work, 0.0_dp, &
                   domo_domo_work)
      CALL cp_fm_to_fm_submat(msource=domo_domo_work, mtarget=img_fm, nrow=ndo_mo, ncol=ndo_mo, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=ndo_mo + 1, t_firstcol=ndo_mo + 1)

      ! alpha-beta block, two non-zero terms in Hx and Hy
      ! Hx term, purely imaginary, addition
      CALL cp_dbcsr_sm_fm_multiply(orb_soc_x, beta_gs_coeffs, ao_domo_work, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, alpha_gs_coeffs, ao_domo_work, 0.0_dp, &
                   domo_domo_work)
      CALL cp_fm_to_fm_submat(msource=domo_domo_work, mtarget=img_fm, nrow=ndo_mo, ncol=ndo_mo, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=1, t_firstcol=ndo_mo + 1)
      ! Hy term, purely real, addition
      CALL cp_dbcsr_sm_fm_multiply(orb_soc_y, beta_gs_coeffs, ao_domo_work, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, alpha_gs_coeffs, ao_domo_work, 0.0_dp, &
                   domo_domo_work)
      CALL cp_fm_to_fm_submat(msource=domo_domo_work, mtarget=real_fm, nrow=ndo_mo, ncol=ndo_mo, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=1, t_firstcol=ndo_mo + 1)

      ! beta-alpha block, two non-zero terms in Hx and Hy
      ! Hx term, purely imaginary, addition
      CALL cp_dbcsr_sm_fm_multiply(orb_soc_x, alpha_gs_coeffs, ao_domo_work, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, 1.0_dp, beta_gs_coeffs, ao_domo_work, 0.0_dp, &
                   domo_domo_work)
      CALL cp_fm_to_fm_submat(msource=domo_domo_work, mtarget=img_fm, nrow=ndo_mo, ncol=ndo_mo, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=ndo_mo + 1, t_firstcol=1)
      ! Hy term, purely real, substraction
      CALL cp_dbcsr_sm_fm_multiply(orb_soc_y, alpha_gs_coeffs, ao_domo_work, ncol=ndo_mo)
      CALL cp_gemm('T', 'N', ndo_mo, ndo_mo, nao, -1.0_dp, beta_gs_coeffs, ao_domo_work, 0.0_dp, &
                   domo_domo_work)
      CALL cp_fm_to_fm_submat(msource=domo_domo_work, mtarget=real_fm, nrow=ndo_mo, ncol=ndo_mo, &
                              s_firstrow=1, s_firstcol=1, t_firstrow=ndo_mo + 1, t_firstcol=1)

      ! Cast everything in complex fm format
      CALL cp_cfm_create(hami_cfm, doso_doso_struct)
      CALL cp_fm_to_cfm(real_fm, img_fm, hami_cfm)

      ! And diagonalize. Since tiny matrix (6x6), diagonalize locally
      ALLOCATE (evals(ndo_so), evecs(ndo_so, ndo_so), hami(ndo_so, ndo_so))
      CALL cp_cfm_get_submatrix(hami_cfm, hami)
      CALL complex_diag(hami, evecs, evals)

      !The SOC corrected KS eigenvalues
      ALLOCATE (tmp_shifts(ndo_mo, 2))

      ialpha = 1; ibeta = 1; 
      DO ido_mo = 1, ndo_so
         !need to find out whether the eigenvalue corresponds to an alpha or beta spin-orbtial
         alpha_tot_contrib = REAL(DOT_PRODUCT(evecs(1:ndo_mo, ido_mo), evecs(1:ndo_mo, ido_mo)))
         beta_tot_contrib = REAL(DOT_PRODUCT(evecs(ndo_mo + 1:ndo_so, ido_mo), evecs(ndo_mo + 1:ndo_so, ido_mo)))

         IF (alpha_tot_contrib > beta_tot_contrib) THEN
            tmp_shifts(ialpha, 1) = evals(ido_mo)
            ialpha = ialpha + 1
         ELSE
            tmp_shifts(ibeta, 2) = evals(ido_mo)
            ibeta = ibeta + 1
         END IF
      END DO

      !compute shift from KS evals
      ALLOCATE (soc_shifts(ndo_mo, SIZE(donor_state%energy_evals, 2)))
      DO ispin = 1, SIZE(donor_state%energy_evals, 2)
         soc_shifts(:, ispin) = tmp_shifts(:, ispin) - donor_state%energy_evals(:, ispin)
      END DO

      ! clean-up
      CALL cp_fm_release(alpha_gs_coeffs)
      CALL cp_fm_release(beta_gs_coeffs)
      CALL cp_fm_release(ao_domo_work)
      CALL cp_fm_release(domo_domo_work)
      CALL cp_fm_release(real_fm)
      CALL cp_fm_release(img_fm)

      CALL cp_cfm_release(hami_cfm)

      CALL cp_fm_struct_release(ao_domo_struct)
      CALL cp_fm_struct_release(domo_domo_struct)
      CALL cp_fm_struct_release(doso_doso_struct)

      CALL timestop(handle)

   END SUBROUTINE get_soc_splitting

END MODULE xas_tdp_correction
