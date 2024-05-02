!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Auxiliary routines needed for RPA-exchange
!>        given blacs_env to another
!> \par History
!>      09.2016 created [Vladimir Rybkin]
!>      03.2019 Renamed [Frederick Stein]
!>      03.2019 Moved Functions from rpa_ri_gpw.F [Frederick Stein]
!>      04.2024 Added open-shell calculations, SOSEX [Frederick Stein]
!> \author Vladimir Rybkin
! **************************************************************************************************
MODULE rpa_exchange
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_p_type,&
                                              cp_fm_struct_release
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat_general,&
                                              cp_fm_type
   USE cp_dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_dot, dbcsr_frobenius_norm, dbcsr_get_info, dbcsr_init_p, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_set, dbcsr_trace, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE hfx_types,                       ONLY: hfx_create,&
                                              hfx_release,&
                                              hfx_type
   USE input_constants,                 ONLY: rpa_exchange_axk,&
                                              rpa_exchange_none,&
                                              rpa_exchange_sosex
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: sqrt2,&
                                              sqrthalf
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type
   USE mp2_types,                       ONLY: mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE rpa_communication,               ONLY: gamma_fm_to_dbcsr
   USE rpa_util,                        ONLY: calc_fm_mat_S_rpa,&
                                              remove_scaling_factor_rpa
   USE scf_control_types,               ONLY: scf_control_type
   USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_exchange'

   PUBLIC :: rpa_exchange_work_type

   TYPE rpa_exchange_env_type
      PRIVATE
      TYPE(qs_environment_type), POINTER :: qs_env => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: mat_hfx => NULL()
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data => NULL()
      TYPE(mp_para_env_type), POINTER                    :: para_env => NULL()
      TYPE(section_vals_type), POINTER                   :: hfx_sections => NULL()
      INTEGER :: n_rep_hf = 0
      INTEGER :: my_group_L_size = 0
      LOGICAL :: my_recalc_hfx_integrals = .FALSE.
      REAL(KIND=dp) :: eps_filter = 0.0_dp
      REAL(KIND=dp) :: e_exchange = 0.0_dp
   CONTAINS
      PROCEDURE, PASS(exchange_env), NON_OVERRIDABLE :: create => hfx_create_subgroup
      PROCEDURE, PASS(exchange_env), NON_OVERRIDABLE :: integrate => integrate_exchange
      PROCEDURE, PASS(exchange_env), NON_OVERRIDABLE :: release => hfx_release_subgroup
   END TYPE

   TYPE rpa_exchange_work_type
      PRIVATE
      INTEGER :: exchange_correction = rpa_exchange_none
      TYPE(rpa_exchange_env_type) :: exchange_env
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: mo_coeff_o
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE :: mo_coeff_v
      INTEGER :: my_group_L_start = 0
      INTEGER :: my_group_L_end = 0
      INTEGER :: my_group_L_size = 0
      INTEGER :: ngroup = 0
      INTEGER :: color_sub = 0
      INTEGER, DIMENSION(:), ALLOCATABLE :: homo, virtual, dimen_ia
      INTEGER :: dimen_RI = 0
      REAL(KIND=dp) :: eps_filter = 0.0_dp
      TYPE(cp_fm_struct_p_type), ALLOCATABLE, DIMENSION(:) :: struct_Gamma
      TYPE(cp_fm_type) :: fm_mat_Q_tmp = cp_fm_type()
      TYPE(cp_fm_type) :: fm_mat_R_half_gemm = cp_fm_type()
      TYPE(cp_fm_type) :: fm_mat_U = cp_fm_type()
      TYPE(mp_para_env_type), POINTER                    :: para_env_RPA => NULL()
      TYPE(mp_para_env_type), POINTER                    :: para_env_sub => NULL()
      INTEGER :: unit_nr = -1
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dbcsr_Gamma_munu_P => NULL()
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: dbcsr_Gamma_inu_P
   CONTAINS
      PROCEDURE, PASS(exchange_work), NON_OVERRIDABLE :: create => rpa_exchange_work_create
      PROCEDURE, PASS(exchange_work), NON_OVERRIDABLE :: compute => rpa_exchange_work_compute
      PROCEDURE, PASS(exchange_work), NON_OVERRIDABLE :: release => rpa_exchange_work_release
   END TYPE

   TYPE dbcsr_matrix_p_set
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:) :: matrix_set
   END TYPE

   LOGICAL, PARAMETER :: use_hfx_implementation = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param exchange_work ...
!> \param mp2_env ...
!> \param para_env_RPA ...
!> \param para_env_sub ...
!> \param dimen_RI ...
!> \param fm_mat_S ...
!> \param fm_mat_Q ...
!> \param fm_mat_Q_gemm ...
!> \param unit_nr ...
!> \param mat_munu ...
! **************************************************************************************************
   SUBROUTINE rpa_exchange_work_create(exchange_work, mp2_env, para_env_RPA, para_env_sub, dimen_RI, &
                                       fm_mat_S, fm_mat_Q, fm_mat_Q_gemm, unit_nr, mat_munu)
      CLASS(rpa_exchange_work_type), INTENT(INOUT) :: exchange_work
      TYPE(mp2_type), INTENT(INOUT) :: mp2_env
      TYPE(mp_para_env_type), TARGET, INTENT(INOUT) :: para_env_RPA, para_env_sub
      INTEGER, INTENT(IN) :: dimen_RI
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN) :: fm_mat_S
      TYPE(cp_fm_type), INTENT(IN) :: fm_mat_Q, fm_mat_Q_gemm
      INTEGER, INTENT(IN) :: unit_nr
      TYPE(dbcsr_type), INTENT(IN)                     :: mat_munu

      INTEGER :: iitmp(2), ispin, nspins

      exchange_work%exchange_correction = mp2_env%ri_rpa%exchange_correction

      IF (exchange_work%exchange_correction == rpa_exchange_none) RETURN

      CALL MOVE_ALLOC(mp2_env%ri_rpa%mo_coeff_o, exchange_work%mo_coeff_o)
      CALL MOVE_ALLOC(mp2_env%ri_rpa%mo_coeff_v, exchange_work%mo_coeff_v)
      exchange_work%eps_filter = mp2_env%mp2_gpw%eps_filter
      exchange_work%para_env_RPA => para_env_RPA
      CALL para_env_RPA%retain()
      exchange_work%para_env_sub => para_env_sub
      CALL para_env_sub%retain()
      exchange_work%unit_nr = unit_nr

      nspins = SIZE(exchange_work%mo_coeff_o)

      ALLOCATE (exchange_work%homo(nspins), exchange_work%virtual(nspins), exchange_work%dimen_ia(nspins))
      DO ispin = 1, nspins
         CALL dbcsr_get_info(exchange_work%mo_coeff_o(ispin), nfullcols_total=exchange_work%homo(ispin))
         CALL dbcsr_get_info(exchange_work%mo_coeff_v(ispin), nfullcols_total=exchange_work%virtual(ispin))
         exchange_work%dimen_ia(ispin) = exchange_work%homo(ispin)*exchange_work%virtual(ispin)
      END DO
      IF (unit_nr > 0) WRITE (unit_nr, *) "dimen_ia", exchange_work%dimen_ia
      exchange_work%dimen_RI = dimen_RI

      ! Get aux sizes
      exchange_work%ngroup = para_env_RPA%num_pe/para_env_sub%num_pe
      exchange_work%color_sub = para_env_RPA%mepos/para_env_sub%num_pe

      iitmp = get_limit(dimen_RI, exchange_work%ngroup, exchange_work%color_sub)
      exchange_work%my_group_L_start = iitmp(1)
      exchange_work%my_group_L_end = iitmp(2)
      exchange_work%my_group_L_size = iitmp(2) - iitmp(1) + 1

      ALLOCATE (exchange_work%struct_Gamma(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_struct_create(exchange_work%struct_Gamma(ispin)%struct, template_fmstruct=fm_mat_S(ispin)%matrix_struct, &
                                  nrow_global=exchange_work%dimen_ia(ispin), ncol_global=exchange_work%dimen_RI)
      END DO

      CALL cp_fm_create(exchange_work%fm_mat_U, fm_mat_Q%matrix_struct, name="fm_mat_U")
      CALL cp_fm_create(exchange_work%fm_mat_Q_tmp, fm_mat_Q%matrix_struct, name="fm_mat_Q_tmp")
      CALL cp_fm_create(exchange_work%fm_mat_R_half_gemm, fm_mat_Q_gemm%matrix_struct)

      ALLOCATE (exchange_work%dbcsr_Gamma_inu_P(nspins))
      CALL dbcsr_allocate_matrix_set(exchange_work%dbcsr_Gamma_munu_P, nspins)
      DO ispin = 1, nspins
         ALLOCATE (exchange_work%dbcsr_Gamma_munu_P(ispin)%matrix)
         CALL dbcsr_create(exchange_work%dbcsr_Gamma_munu_P(ispin)%matrix, template=mat_munu, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(exchange_work%dbcsr_Gamma_munu_P(ispin)%matrix, mat_munu)

         CALL dbcsr_create(exchange_work%dbcsr_Gamma_inu_P(ispin), template=exchange_work%mo_coeff_o(ispin))
         CALL dbcsr_copy(exchange_work%dbcsr_Gamma_inu_P(ispin), exchange_work%mo_coeff_o(ispin))
         CALL dbcsr_set(exchange_work%dbcsr_Gamma_inu_P(ispin), 0.0_dp)
      END DO
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param exchange_work ...
! **************************************************************************************************
   SUBROUTINE rpa_exchange_work_release(exchange_work)
      CLASS(rpa_exchange_work_type), INTENT(INOUT) :: exchange_work

      INTEGER :: ispin

      IF (ALLOCATED(exchange_work%homo)) DEALLOCATE (exchange_work%homo)
      IF (ALLOCATED(exchange_work%virtual)) DEALLOCATE (exchange_work%virtual)
      IF (ALLOCATED(exchange_work%dimen_ia)) DEALLOCATE (exchange_work%dimen_ia)
      IF (ALLOCATED(exchange_work%struct_Gamma)) THEN
      DO ispin = 1, SIZE(exchange_work%struct_Gamma)
         CALL cp_fm_struct_release(exchange_work%struct_Gamma(ispin)%struct)
      END DO
      DEALLOCATE (exchange_work%struct_Gamma)
      END IF
      CALL cp_fm_release(exchange_work%fm_mat_Q_tmp)
      CALL cp_fm_release(exchange_work%fm_mat_U)
      CALL cp_fm_release(exchange_work%fm_mat_R_half_gemm)
      CALL mp_para_env_release(exchange_work%para_env_RPA)
      CALL mp_para_env_release(exchange_work%para_env_sub)

      IF (ALLOCATED(exchange_work%mo_coeff_o)) THEN
         DO ispin = 1, SIZE(exchange_work%mo_coeff_o)
            CALL dbcsr_release(exchange_work%mo_coeff_o(ispin))
         END DO
         DEALLOCATE (exchange_work%mo_coeff_o)
      END IF
      IF (ALLOCATED(exchange_work%mo_coeff_v)) THEN
         DO ispin = 1, SIZE(exchange_work%mo_coeff_v)
            CALL dbcsr_release(exchange_work%mo_coeff_v(ispin))
         END DO
         DEALLOCATE (exchange_work%mo_coeff_v)
      END IF

      IF (ALLOCATED(exchange_work%dbcsr_Gamma_inu_P)) THEN
         DO ispin = 1, SIZE(exchange_work%dbcsr_Gamma_inu_P)
            CALL dbcsr_release(exchange_work%dbcsr_Gamma_inu_P(ispin))
         END DO
         DEALLOCATE (exchange_work%dbcsr_Gamma_inu_P)
      END IF
      IF (ASSOCIATED(exchange_work%dbcsr_Gamma_munu_P)) THEN
         DO ispin = 1, SIZE(exchange_work%dbcsr_Gamma_munu_P)
            CALL dbcsr_release(exchange_work%dbcsr_Gamma_munu_P(ispin)%matrix)
            DEALLOCATE (exchange_work%dbcsr_Gamma_munu_P(ispin)%matrix)
         END DO
         DEALLOCATE (exchange_work%dbcsr_Gamma_munu_P)
      END IF
   END SUBROUTINE

! **************************************************************************************************
!> \brief Main driver for RPA-exchange energies
!> \param exchange_work ...
!> \param qs_env ...
!> \param fm_mat_Q ...
!> \param eig ...
!> \param fm_mat_S ...
!> \param omega ...
!> \param e_exchange_corr ...          exchange energy correctrion for a quadrature point
!> \author Vladimir Rybkin, 07/2016
! **************************************************************************************************
   SUBROUTINE rpa_exchange_work_compute(exchange_work, qs_env, fm_mat_Q, &
                                        eig, fm_mat_S, omega, &
                                        e_exchange_corr)
      CLASS(rpa_exchange_work_type), INTENT(INOUT) :: exchange_work
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_Q
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eig
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_S
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), INTENT(INOUT)                       :: e_exchange_corr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rpa_exchange_work_compute'
      REAL(KIND=dp), PARAMETER                           :: thresh = 0.0000001_dp

      INTEGER :: handle, iib, ispin, L_counter, nspins, P, Q, L_size_B, L_size_Gamma, hom, virt, i, j, a, b
      REAL(KIND=dp) :: e_rpa, sum_a, sum_b, sum_ia, sum_jb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: product_matrix_small
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :) :: mat_Gamma_3, mat_B_iaP
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :) :: product_matrix
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_mat_Gamma_3
      TYPE(dbcsr_matrix_p_set), ALLOCATABLE, &
         DIMENSION(:)                                    :: dbcsr_Gamma_3
      TYPE(mp_comm_type) :: comm_exchange

      CALL timeset(routineN, handle)

      MARK_USED(product_matrix)
      MARK_USED(product_matrix_small)
      MARK_USED(sqrt2)

      IF (exchange_work%exchange_correction == rpa_exchange_none) RETURN

      nspins = SIZE(fm_mat_S)

      ! Eigenvalues
      ALLOCATE (eigenval(exchange_work%dimen_RI))
      eigenval = 0.0_dp

      CALL cp_fm_set_all(matrix=exchange_work%fm_mat_Q_tmp, alpha=0.0_dp)
      CALL cp_fm_set_all(matrix=exchange_work%fm_mat_U, alpha=0.0_dp)

      ! Copy Q to Q_tmp
      CALL cp_fm_to_fm(fm_mat_Q, exchange_work%fm_mat_Q_tmp)
      ! Diagonalize Q
      CALL choose_eigv_solver(exchange_work%fm_mat_Q_tmp, exchange_work%fm_mat_U, eigenval)

      IF (exchange_work%unit_nr > 0) WRITE (exchange_work%unit_nr, *) &
         "DEBUG RPA energy (AXK)", 0.5_dp*SUM(LOG(1.0_dp + eigenval) - eigenval)

      ! Calculate diagonal matrix for R_half

      ! Manipulate eigenvalues to get diagonal matrix
      IF (exchange_work%exchange_correction == rpa_exchange_axk) THEN
         IF (exchange_work%unit_nr > 0) WRITE (exchange_work%unit_nr, *) "Do AXK"
         DO iib = 1, exchange_work%dimen_RI
            IF (ABS(eigenval(iib)) .GE. thresh) THEN
               eigenval(iib) = &
                  SQRT((1.0_dp/(eigenval(iib)**2))*LOG(1.0_dp + eigenval(iib)) &
                       - 1.0_dp/(eigenval(iib)*(eigenval(iib) + 1.0_dp)))
            ELSE
               eigenval(iib) = sqrthalf
            END IF
         END DO
      ELSE IF (exchange_work%exchange_correction == rpa_exchange_sosex) THEN
         IF (exchange_work%unit_nr > 0) WRITE (exchange_work%unit_nr, *) "Do SOSEX"
         DO iib = 1, exchange_work%dimen_RI
            IF (ABS(eigenval(iib)) .GE. thresh) THEN
               eigenval(iib) = &
                  SQRT(-(1.0_dp/(eigenval(iib)**2))*LOG(1.0_dp + eigenval(iib)) &
                       + 1.0_dp/eigenval(iib))
            ELSE
               eigenval(iib) = sqrthalf
            END IF
         END DO
      ELSE
         CPABORT("Unknown RPA exchange correction")
      END IF

      ! fm_mat_U now contains some sqrt of the required matrix-valued function
      CALL cp_fm_column_scale(exchange_work%fm_mat_U, eigenval)

      ! Release memory
      DEALLOCATE (eigenval)

      ! Redistribute fm_mat_U for "rectangular" multiplication: ia*P P*P
      CALL cp_fm_set_all(matrix=exchange_work%fm_mat_R_half_gemm, alpha=0.0_dp)

      CALL cp_fm_to_fm_submat_general(exchange_work%fm_mat_U, exchange_work%fm_mat_R_half_gemm, exchange_work%dimen_RI, &
                                      exchange_work%dimen_RI, 1, 1, 1, 1, &
                                      exchange_work%fm_mat_U%matrix_struct%context)

      ALLOCATE (fm_mat_Gamma_3(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_create(fm_mat_Gamma_3(ispin), exchange_work%struct_Gamma(ispin)%struct)
         CALL cp_fm_set_all(matrix=fm_mat_Gamma_3(ispin), alpha=0.0_dp)

         ! Update G with a new value of Omega: in practice, it is G*S

         ! Scale fm_work_iaP
         CALL calc_fm_mat_S_rpa(fm_mat_S(ispin), .TRUE., exchange_work%virtual(ispin), eig(:, ispin), &
                                exchange_work%homo(ispin), omega, 0.0_dp)

         ! Calculate Gamma_3: Gamma_3 = G*S*R^(1/2) = G*S*R^(1/2)
         CALL parallel_gemm(transa="T", transb="N", m=exchange_work%dimen_ia(ispin), n=exchange_work%dimen_RI, &
                            k=exchange_work%dimen_RI, alpha=1.0_dp, &
                            matrix_a=fm_mat_S(ispin), matrix_b=exchange_work%fm_mat_R_half_gemm, beta=0.0_dp, &
                            matrix_c=fm_mat_Gamma_3(ispin))

         ! Remove extra factor from S after the multiplication (to return to the original matrix)
         CALL remove_scaling_factor_rpa(fm_mat_S(ispin), exchange_work%virtual(ispin), &
                                        eig(:, ispin), exchange_work%homo(ispin), omega)
      END DO

      IF (use_hfx_implementation) THEN
         CALL exchange_work%exchange_env%create(qs_env, exchange_work%dbcsr_gamma_munu_P(1)%matrix, exchange_work%para_env_sub, &
                                                exchange_work%my_group_L_size, exchange_work%eps_filter)
         exchange_work%exchange_env%e_exchange = 0.0_dp

         ALLOCATE (dbcsr_Gamma_3(nspins))
         DO ispin = 1, nspins
            ! Copy Gamma_ia_P^3 to dbcsr matrix set
            CALL gamma_fm_to_dbcsr(fm_mat_Gamma_3(ispin), dbcsr_Gamma_3(ispin)%matrix_set, &
                                   exchange_work%para_env_RPA, exchange_work%para_env_sub, &
                                   exchange_work%homo(ispin), exchange_work%virtual(ispin), exchange_work%mo_coeff_o(ispin), &
                                   exchange_work%ngroup, exchange_work%my_group_L_start, &
                                   exchange_work%my_group_L_end, exchange_work%my_group_L_size)
         END DO
         DEALLOCATE (fm_mat_Gamma_3)

         DO L_counter = 1, exchange_work%my_group_L_size
            DO ispin = 1, nspins
               ! Do dbcsr multiplication: transform the virtual index
               CALL dbcsr_multiply("N", "T", 1.0_dp, exchange_work%mo_coeff_v(ispin), dbcsr_Gamma_3(ispin)%matrix_set(L_counter), &
                                   0.0_dp, exchange_work%dbcsr_Gamma_inu_P(ispin), filter_eps=exchange_work%eps_filter)

               CALL dbcsr_release(dbcsr_Gamma_3(ispin)%matrix_set(L_counter))

               ! Do dbcsr multiplication: transform the occupied index
               CALL dbcsr_multiply("N", "T", 1.0_dp, exchange_work%dbcsr_Gamma_inu_P(ispin), exchange_work%mo_coeff_o(ispin), &
                                   0.0_dp, exchange_work%dbcsr_Gamma_munu_P(ispin)%matrix, filter_eps=exchange_work%eps_filter)
            END DO

            ! Contract DM with exchange integrals
            CALL exchange_work%exchange_env%integrate(exchange_work%dbcsr_Gamma_munu_P)
         END DO

         DO ispin = 1, nspins
            DEALLOCATE (dbcsr_Gamma_3(ispin)%matrix_set)
         END DO
         DEALLOCATE (dbcsr_Gamma_3)

         CALL exchange_work%exchange_env%release(e_exchange_corr)

         IF (nspins == 1) e_exchange_corr = e_exchange_corr
         IF (nspins == 2) e_exchange_corr = e_exchange_corr/sqrthalf

         ! For printing, we need the energy per quadrature point/per integration group, but we obtain the energy per sub group
         e_exchange_corr = e_exchange_corr/REAL(exchange_work%para_env_sub%num_pe, KIND=dp)
         CALL exchange_work%para_env_RPA%sum(e_exchange_corr)
      ELSE
         ! Reorder matrix from (P, i*a) -> (a, i, P) with P being distributed within integration group
         CALL comm_exchange%from_split(exchange_work%para_env_RPA, exchange_work%color_sub)
         e_exchange_corr = 0.0_dp
         e_rpa = 0.0_dp

         DO ispin = 1, nspins
            hom = exchange_work%homo(ispin)
            virt = exchange_work%virtual(ispin)

            CALL redistribute_Gamma_in_subgroups(fm_mat_Gamma_3(ispin), mat_Gamma_3, comm_exchange, &
                                                 hom, virt)
            CALL cp_fm_release(fm_mat_Gamma_3(ispin))

            ! We need only the pure matrix
            CALL remove_scaling_factor_rpa(fm_mat_S(ispin), virt, &
                                           eig(:, ispin), hom, omega)

            ! Reorder matrix from (P, i*a) -> (a, i, P) with P being distributed within subgroups
            CALL redistribute_B_in_subgroups(fm_mat_S(ispin), mat_B_iaP, exchange_work%para_env_sub, &
                                             hom, virt)

            ! Return to the original tensor
            CALL calc_fm_mat_S_rpa(fm_mat_S(ispin), .TRUE., virt, eig(:, ispin), hom, omega, 0.0_dp)

            L_size_Gamma = SIZE(mat_Gamma_3, 3)
            L_size_B = SIZE(mat_B_iaP, 3)

            WRITE (*, *) "sizes", L_size_B, L_size_Gamma, hom, virt
            DO P = 1, L_size_Gamma
            DO Q = 1, L_size_B
               sum_ia = 0.0_dp
               DO i = 1, hom
               DO a = 1, virt
                  sum_ia = sum_ia + mat_Gamma_3(a, i, P)*mat_B_iaP(a, i, Q)
               END DO
               END DO
               sum_jb = 0.0_dp
               DO j = 1, hom
               DO b = 1, virt
                  sum_jb = sum_jb + mat_B_iaP(b, j, Q)*mat_Gamma_3(b, j, P)
               END DO
               END DO
               e_rpa = e_rpa + sum_ia*sum_jb
            END DO
            END DO
            DO P = 1, L_size_Gamma
            DO Q = 1, L_size_B
            DO i = 1, hom
            DO j = 1, hom
               sum_a = 0.0_dp
               DO a = 1, virt
                  sum_a = sum_a + mat_Gamma_3(a, i, P)*mat_B_iaP(a, j, Q)
               END DO
               sum_b = 0.0_dp
               DO b = 1, virt
                  sum_b = sum_b + mat_B_iaP(b, i, Q)*mat_Gamma_3(b, j, P)
               END DO
               e_exchange_corr = e_exchange_corr + sum_a*sum_b
            END DO
            END DO
            END DO
            END DO

            DEALLOCATE (mat_B_iaP, mat_Gamma_3)
            WRITE (*, *) ispin, e_rpa, e_exchange_corr
         END DO
         IF (exchange_work%unit_nr > 0) WRITE (exchange_work%unit_nr, *) "E_rpa/E_ex", e_rpa, e_exchange_corr
         CALL comm_exchange%free()
         DEALLOCATE (fm_mat_Gamma_3)

         e_exchange_corr = e_exchange_corr*2.0_dp
      END IF

      ! Print exchange correlation energy to the file
      IF (exchange_work%unit_nr > 0) WRITE (exchange_work%unit_nr, '(T3,A,T68,F25.14,A4)') &
         'exchange correlation energy for a quadrature point:', &
         e_exchange_corr, ' a.u.'

      ! We return an energy per process as we will have to sum over all RPA integration groups anyways
      e_exchange_corr = e_exchange_corr/REAL(exchange_work%para_env_RPA%num_pe, KIND=dp)

      CALL timestop(handle)

   END SUBROUTINE rpa_exchange_work_compute

! **************************************************************************************************
!> \brief Contract RPA-exchange density matrix with HF exchange integrals and evaluate the correction
!> \param exchange_env ...
!> \param dbcsr_Gamma_munu_P ...  exchange density matrix in AO basis to be contracted
!> \author Vladimir Rybkin, 08/2016
! **************************************************************************************************
   SUBROUTINE integrate_exchange(exchange_env, dbcsr_Gamma_munu_P)
      CLASS(rpa_exchange_env_type), INTENT(INOUT) :: exchange_env
      TYPE(dbcsr_p_type), DIMENSION(:), TARGET          :: dbcsr_Gamma_munu_P

      CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_exchange'

      INTEGER                                            :: handle, ispin
      REAL(KIND=dp)                                      :: e_exchange_P

      CALL timeset(routineN, handle)

      DO ispin = 1, SIZE(exchange_env%mat_hfx)
         CALL dbcsr_set(exchange_env%mat_hfx(ispin)%matrix, 0.0_dp)
      END DO

      CALL tddft_hfx_matrix(exchange_env%mat_hfx, dbcsr_Gamma_munu_P, exchange_env%qs_env, .FALSE., &
                            exchange_env%my_recalc_hfx_integrals, &
                            exchange_env%hfx_sections, exchange_env%x_data, exchange_env%para_env)

      exchange_env%my_recalc_hfx_integrals = .FALSE.
      DO ispin = 1, SIZE(dbcsr_Gamma_munu_P)
         CALL dbcsr_dot(exchange_env%mat_hfx(ispin)%matrix, dbcsr_Gamma_munu_P(ispin)%matrix, e_exchange_p)
         exchange_env%e_exchange = exchange_env%e_exchange - e_exchange_P
      END DO

      CALL timestop(handle)

   END SUBROUTINE integrate_exchange

! **************************************************************************************************
!> \brief ... Initializes x_data on a subgroup
!> \param exchange_env ...
!> \param qs_env ...
!> \param mat_munu ...
!> \param para_env_sub ...
!> \param my_group_L_size ...
!> \param eps_filter ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE hfx_create_subgroup(exchange_env, qs_env, mat_munu, para_env_sub, my_group_L_size, eps_filter)
      CLASS(rpa_exchange_env_type), INTENT(INOUT) :: exchange_env
      TYPE(dbcsr_type), INTENT(IN) :: mat_munu
      TYPE(qs_environment_type), POINTER, INTENT(IN)   :: qs_env
      TYPE(mp_para_env_type), POINTER, INTENT(IN)            :: para_env_sub
      INTEGER, INTENT(IN) :: my_group_L_size
      REAL(KIND=dp), INTENT(IN) :: eps_filter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_create_subgroup'

      INTEGER                                            :: handle, nelectron_total, ispin
      LOGICAL                                            :: do_hfx
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: my_cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      exchange_env%qs_env => qs_env
      exchange_env%para_env => para_env_sub
      exchange_env%eps_filter = eps_filter
      exchange_env%my_group_L_size = my_group_L_size

      NULLIFY (my_cell, atomic_kind_set, particle_set, dft_control, qs_kind_set, scf_control)

      CALL get_qs_env(qs_env, &
                      subsys=subsys, &
                      input=input, &
                      scf_control=scf_control, &
                      nelectron_total=nelectron_total)

      CALL qs_subsys_get(subsys, &
                         cell=my_cell, &
                         atomic_kind_set=atomic_kind_set, &
                         qs_kind_set=qs_kind_set, &
                         particle_set=particle_set)

      do_hfx = .TRUE.
      exchange_env%hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")
      CALL section_vals_get(exchange_env%hfx_sections, explicit=do_hfx, n_repetition=exchange_env%n_rep_hf)
      CALL get_qs_env(qs_env, dft_control=dft_control)

      IF (do_hfx) THEN
         ! Retrieve particle_set and atomic_kind_set
         CALL hfx_create(exchange_env%x_data, para_env_sub, exchange_env%hfx_sections, atomic_kind_set, &
                         qs_kind_set, particle_set, dft_control, my_cell, orb_basis='ORB', &
                         nelectron_total=nelectron_total)
      END IF

      exchange_env%my_recalc_hfx_integrals = .TRUE.

      CALL dbcsr_allocate_matrix_set(exchange_env%mat_hfx, dft_control%nspins)
      DO ispin = 1, dft_control%nspins
         ALLOCATE (exchange_env%mat_hfx(ispin)%matrix)
         CALL dbcsr_init_p(exchange_env%mat_hfx(ispin)%matrix)
         CALL dbcsr_create(exchange_env%mat_hfx(ispin)%matrix, template=mat_munu, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(exchange_env%mat_hfx(ispin)%matrix, mat_munu)
      END DO

      CALL timestop(handle)

   END SUBROUTINE hfx_create_subgroup

! **************************************************************************************************
!> \brief ...
!> \param exchange_env ...
!> \param e_exchange ...
! **************************************************************************************************
   SUBROUTINE hfx_release_subgroup(exchange_env, e_exchange)
      CLASS(rpa_exchange_env_type), INTENT(INOUT) :: exchange_env
      REAL(KIND=dp), INTENT(OUT) :: e_exchange

      INTEGER :: ispin

      NULLIFY (exchange_env%qs_env, exchange_env%para_env, exchange_env%hfx_sections)

      e_exchange = exchange_env%e_exchange
      exchange_env%e_exchange = 0.0_dp

      CALL hfx_release(exchange_env%x_data)
      NULLIFY (exchange_env%x_data)

      DO ispin = 1, SIZE(exchange_env%mat_hfx, 1)
         CALL dbcsr_release(exchange_env%mat_hfx(ispin)%matrix)
         CALL dbcsr_release(exchange_env%mat_hfx(ispin)%matrix)
         DEALLOCATE (exchange_env%mat_hfx(ispin)%matrix)
      END DO
      DEALLOCATE (exchange_env%mat_hfx)
      NULLIFY (exchange_env%mat_hfx)
   END SUBROUTINE hfx_release_subgroup

! **************************************************************************************************
!> \brief ...
!> \param fm_mat ...
!> \param mat ...
!> \param comm ...
!> \param ncols ...
!> \param nrows ...
! **************************************************************************************************
   SUBROUTINE redistribute_Gamma_in_subgroups(fm_mat, mat, comm, ncols, nrows)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: mat

      CLASS(mp_comm_type), INTENT(IN) :: comm
      INTEGER, INTENT(IN) :: ncols, nrows

      INTEGER :: iocc, avirt, ncol_global

      MARK_USED(comm)

      CALL cp_fm_get_info(matrix=fm_mat, &
                          ncol_global=ncol_global)

      ALLOCATE (mat(nrows, ncols, ncol_global))

      DO iocc = 1, ncols
      DO avirt = 1, nrows
         mat(avirt, iocc, :) = fm_mat%local_data((iocc - 1)*nrows + avirt, :)
      END DO
      END DO

      WRITE (*, *) "DEBUG redist", NORM2(mat) - NORM2(fm_mat%local_data)

   END SUBROUTINE redistribute_Gamma_in_subgroups

! **************************************************************************************************
!> \brief ...
!> \param fm_mat ...
!> \param mat ...
!> \param comm ...
!> \param ncols ...
!> \param nrows ...
! **************************************************************************************************
   SUBROUTINE redistribute_B_in_subgroups(fm_mat, mat, comm, ncols, nrows)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: mat

      CLASS(mp_comm_type), INTENT(IN) :: comm
      INTEGER, INTENT(IN) :: ncols, nrows

      INTEGER :: iocc, avirt, nrow_global

      MARK_USED(comm)

      CALL cp_fm_get_info(matrix=fm_mat, &
                          nrow_global=nrow_global)

      ALLOCATE (mat(nrows, ncols, nrow_global))

      DO iocc = 1, ncols
      DO avirt = 1, nrows
         mat(avirt, iocc, :) = fm_mat%local_data(:, (iocc - 1)*nrows + avirt)
      END DO
      END DO

      WRITE (*, *) "DEBUG redist", NORM2(mat) - NORM2(fm_mat%local_data)

   END SUBROUTINE redistribute_B_in_subgroups

END MODULE rpa_exchange
