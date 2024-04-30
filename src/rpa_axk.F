!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Auxiliary routines needed for RPA-AXK
!>        given blacs_env to another
!> \par History
!>      09.2016 created [Vladimir Rybkin]
!>      03.2019 Renamed [Frederick Stein]
!>      03.2019 Moved Functions from rpa_ri_gpw.F [Frederick Stein]
!> \author Vladimir Rybkin
! **************************************************************************************************
MODULE rpa_axk
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_dot, dbcsr_get_info, dbcsr_init_p, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_release, dbcsr_set, dbcsr_trace, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat_general,&
                                              cp_fm_type
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE hfx_types,                       ONLY: hfx_create,&
                                              hfx_release,&
                                              hfx_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: sqrthalf
   USE message_passing,                 ONLY: mp_para_env_type
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_axk'

   PUBLIC :: rpa_axk_work_type

   TYPE rpa_axk_exchange_type
      PRIVATE
      TYPE(qs_environment_type), POINTER :: qs_env => NULL()
      TYPE(dbcsr_type)          :: work_ao
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: mat_hfx => NULL()
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data => NULL()
      TYPE(mp_para_env_type), POINTER                    :: para_env => NULL()
      TYPE(section_vals_type), POINTER                   :: hfx_sections => NULL()
      INTEGER :: n_rep_hf = 0
      INTEGER :: my_group_L_size = 0
      LOGICAL :: my_recalc_hfx_integrals = .FALSE.
      REAL(KIND=dp) :: eps_filter = 0.0_dp
      REAL(KIND=dp) :: e_axk = 0.0_dp
   CONTAINS
      PROCEDURE, PASS(axk_exchange), NON_OVERRIDABLE :: create => hfx_create_subgroup
      PROCEDURE, PASS(axk_exchange), NON_OVERRIDABLE :: integrate => integrate_exchange
      PROCEDURE, PASS(axk_exchange), NON_OVERRIDABLE :: release => hfx_release_subgroup
   END TYPE

   TYPE rpa_axk_work_type
      PRIVATE
      TYPE(rpa_axk_exchange_type) :: axk_exchange
      TYPE(dbcsr_type), DIMENSION(:), POINTER :: mo_coeff_o => NULL()
      TYPE(dbcsr_type), DIMENSION(:), POINTER :: mo_coeff_v => NULL()
      INTEGER :: my_group_L_start = 0
      INTEGER :: my_group_L_end = 0
      INTEGER :: my_group_L_size = 0
      INTEGER :: ngroup = 0
      INTEGER, DIMENSION(:), ALLOCATABLE :: homo, virtual, dimen_ia
      INTEGER :: dimen_RI = 0
      REAL(KIND=dp) :: eps_filter = 0.0_dp
      TYPE(cp_fm_struct_type), POINTER :: struct_Gamma => NULL()
      TYPE(cp_fm_type) :: fm_mat_Q_tmp = cp_fm_type()
      TYPE(cp_fm_type) :: fm_mat_R_half = cp_fm_type()
      TYPE(cp_fm_type) :: fm_mat_R_half_gemm = cp_fm_type()
      TYPE(cp_fm_type) :: fm_mat_U = cp_fm_type()
      TYPE(mp_para_env_type), POINTER                    :: para_env_RPA => NULL()
      TYPE(mp_para_env_type), POINTER                    :: para_env_sub => NULL()
      INTEGER :: unit_nr = -1
   CONTAINS
      PROCEDURE, PASS(axk_work), NON_OVERRIDABLE :: create => rpa_axk_work_create
      PROCEDURE, PASS(axk_work), NON_OVERRIDABLE :: compute => rpa_axk_work_compute
      PROCEDURE, PASS(axk_work), NON_OVERRIDABLE :: release => rpa_axk_work_release
   END TYPE

   TYPE dbcsr_matrix_p_set
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:) :: matrix_set
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param axk_work ...
!> \param mp2_env ...
!> \param para_env_RPA ...
!> \param para_env_sub ...
!> \param dimen_RI ...
!> \param fm_mat_S ...
!> \param fm_mat_Q ...
!> \param fm_mat_Q_gemm ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE rpa_axk_work_create(axk_work, mp2_env, para_env_RPA, para_env_sub, dimen_RI, &
                                  fm_mat_S, fm_mat_Q, fm_mat_Q_gemm, unit_nr)
      CLASS(rpa_axk_work_type), INTENT(INOUT) :: axk_work
      TYPE(mp2_type), INTENT(IN) :: mp2_env
      TYPE(mp_para_env_type), TARGET, INTENT(IN) :: para_env_RPA, para_env_sub
      INTEGER, INTENT(IN) :: dimen_RI
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN) :: fm_mat_S
      TYPE(cp_fm_type), INTENT(IN) :: fm_mat_Q, fm_mat_Q_gemm
      INTEGER, INTENT(IN) :: unit_nr

      INTEGER :: color_sub, iitmp(2), ispin, nspins

      axk_work%mo_coeff_o => mp2_env%ri_rpa%mo_coeff_o
      axk_work%mo_coeff_v => mp2_env%ri_rpa%mo_coeff_v
      axk_work%eps_filter = mp2_env%mp2_gpw%eps_filter
      axk_work%para_env_RPA => para_env_RPA
      axk_work%para_env_sub => para_env_sub
      axk_work%unit_nr = unit_nr

      nspins = SIZE(axk_work%mo_coeff_o)

      ALLOCATE (axk_work%homo(nspins), axk_work%virtual(nspins), axk_work%dimen_ia(nspins))
      DO ispin = 1, nspins
         CALL dbcsr_get_info(axk_work%mo_coeff_o(ispin), nfullcols_total=axk_work%homo(ispin))
         CALL dbcsr_get_info(axk_work%mo_coeff_v(ispin), nfullcols_total=axk_work%virtual(ispin))
         axk_work%dimen_ia(ispin) = axk_work%homo(ispin)*axk_work%virtual(ispin)
      END DO
      axk_work%dimen_RI = dimen_RI

      ! Get aux sizes
      axk_work%ngroup = para_env_RPA%num_pe/para_env_sub%num_pe
      color_sub = para_env_RPA%mepos/para_env_sub%num_pe

      iitmp = get_limit(dimen_RI, axk_work%ngroup, color_sub)
      axk_work%my_group_L_start = iitmp(1)
      axk_work%my_group_L_end = iitmp(2)
      axk_work%my_group_L_size = iitmp(2) - iitmp(1) + 1

      DO ispin = 1, nspins
         NULLIFY (axk_work%struct_Gamma)
         CALL cp_fm_struct_create(axk_work%struct_Gamma, template_fmstruct=fm_mat_S(ispin)%matrix_struct, &
                                  nrow_global=axk_work%dimen_ia(ispin), ncol_global=axk_work%dimen_RI)
      END DO

      CALL cp_fm_create(axk_work%fm_mat_U, fm_mat_Q%matrix_struct, name="fm_mat_U")
      CALL cp_fm_create(axk_work%fm_mat_R_half, fm_mat_Q%matrix_struct, name="fm_mat_R_half")
      CALL cp_fm_create(axk_work%fm_mat_Q_tmp, fm_mat_Q%matrix_struct, name="fm_mat_Q_tmp")
      CALL cp_fm_create(axk_work%fm_mat_R_half_gemm, fm_mat_Q_gemm%matrix_struct)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param axk_work ...
! **************************************************************************************************
   SUBROUTINE rpa_axk_work_release(axk_work)
      CLASS(rpa_axk_work_type), INTENT(INOUT) :: axk_work

      IF (ALLOCATED(axk_work%homo)) DEALLOCATE (axk_work%homo)
      IF (ALLOCATED(axk_work%virtual)) DEALLOCATE (axk_work%virtual)
      IF (ALLOCATED(axk_work%dimen_ia)) DEALLOCATE (axk_work%dimen_ia)
      IF (ASSOCIATED(axk_work%struct_Gamma)) CALL cp_fm_struct_release(axk_work%struct_Gamma)
      CALL cp_fm_release(axk_work%fm_mat_Q_tmp)
      CALL cp_fm_release(axk_work%fm_mat_U)
      CALL cp_fm_release(axk_work%fm_mat_R_half)
      CALL cp_fm_release(axk_work%fm_mat_R_half_gemm)
      NULLIFY (axk_work%para_env_RPA, axk_work%para_env_sub)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Main driver for RPA-AXK energies
!> \param axk_work ...
!> \param qs_env ...
!> \param fm_mat_Q ...
!> \param eig ...
!> \param fm_mat_S ...
!> \param omega ...
!> \param mat_munu ...
!> \param e_axk_corr ...          AXK energy correctrion for a quadrature point
!> \author Vladimir Rybkin, 07/2016
! **************************************************************************************************
   SUBROUTINE rpa_axk_work_compute(axk_work, qs_env, fm_mat_Q, &
                                   eig, fm_mat_S, omega, &
                                   mat_munu, e_axk_corr)
      CLASS(rpa_axk_work_type), INTENT(INOUT) :: axk_work
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_Q
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eig
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_S
      REAL(KIND=dp), INTENT(IN)                          :: omega
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_munu
      REAL(KIND=dp), INTENT(INOUT)                       :: e_axk_corr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rpa_axk_work_compute'
      REAL(KIND=dp), PARAMETER                           :: thresh = 0.0000001_dp

      INTEGER :: handle, iib, ispin, L_counter, nspins
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_mat_Gamma_3
      TYPE(dbcsr_matrix_p_set), ALLOCATABLE, &
         DIMENSION(:)                                    :: dbcsr_Gamma_3
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dbcsr_Gamma_munu_P
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: dbcsr_Gamma_inu_P

      CALL timeset(routineN, handle)

      NULLIFY (dbcsr_Gamma_munu_P)

      nspins = SIZE(fm_mat_S)

      ! Eigenvalues
      ALLOCATE (eigenval(axk_work%dimen_RI))
      eigenval = 0.0_dp

      CALL cp_fm_set_all(matrix=axk_work%fm_mat_Q_tmp, alpha=0.0_dp)
      CALL cp_fm_set_all(matrix=axk_work%fm_mat_U, alpha=0.0_dp)
      CALL cp_fm_set_all(matrix=axk_work%fm_mat_R_half, alpha=0.0_dp)

      ! Copy Q to Q_tmp
      CALL cp_fm_to_fm(fm_mat_Q, axk_work%fm_mat_Q_tmp)

      CALL cp_fm_scale(0.50_dp, axk_work%fm_mat_Q_tmp)
      ! Diagonalize Q
      CALL choose_eigv_solver(axk_work%fm_mat_Q_tmp, axk_work%fm_mat_U, eigenval)

      ! Calculate diagonal matrix for R_half

      ! U*diag stored in U, whereas eigenvectors are in fm_mat_Q_tmp
      CALL cp_fm_to_fm(axk_work%fm_mat_U, axk_work%fm_mat_Q_tmp)

      ! Manipulate eigenvalues to get diagonal matrix
      DO iib = 1, axk_work%dimen_RI
         IF (ABS(eigenval(iib)) .GE. thresh) THEN
            eigenval(iib) = &
               SQRT((1.0_dp/(eigenval(iib)**2))*LOG(1.0_dp + eigenval(iib)) &
                    - 1.0_dp/(eigenval(iib)*(eigenval(iib) + 1.0_dp)))
         ELSE
            eigenval(iib) = sqrthalf
         END IF
      END DO

      CALL cp_fm_column_scale(axk_work%fm_mat_U, eigenval)

      ! Release memory
      DEALLOCATE (eigenval)

      ! Get R_half by multiplication
      CALL parallel_gemm(transa="N", transb="T", m=axk_work%dimen_RI, n=axk_work%dimen_RI, k=axk_work%dimen_RI, alpha=1.0_dp, &
                         matrix_a=axk_work%fm_mat_U, matrix_b=axk_work%fm_mat_Q_tmp, beta=0.0_dp, &
                         matrix_c=axk_work%fm_mat_R_half)

      ! Redistribute fm_mat_R_half for "rectangular" multiplication: ia*P P*P
      CALL cp_fm_set_all(matrix=axk_work%fm_mat_R_half_gemm, alpha=0.0_dp)

      CALL cp_fm_to_fm_submat_general(axk_work%fm_mat_R_half, axk_work%fm_mat_R_half_gemm, axk_work%dimen_RI, &
                                      axk_work%dimen_RI, 1, 1, 1, 1, &
                                      axk_work%fm_mat_R_half%matrix_struct%context)

      ALLOCATE (fm_mat_Gamma_3(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_create(fm_mat_Gamma_3(ispin), axk_work%struct_Gamma)
         CALL cp_fm_set_all(matrix=fm_mat_Gamma_3(ispin), alpha=0.0_dp)

         ! Update G with a new value of Omega: in practice, it is G*S

         ! Scale fm_work_iaP
         CALL calc_fm_mat_S_rpa(fm_mat_S(ispin), .TRUE., axk_work%virtual(ispin), eig(:, ispin), &
                                axk_work%homo(ispin), omega, 0.0_dp)

         ! Calculate Gamma_3: Gamma_3 = G*S*R^(1/2) = G*S*R^(1/2) )
         CALL parallel_gemm(transa="T", transb="N", m=axk_work%dimen_ia(ispin), n=axk_work%dimen_RI, &
                            k=axk_work%dimen_RI, alpha=1.0_dp, &
                            matrix_a=fm_mat_S(ispin), matrix_b=axk_work%fm_mat_R_half_gemm, beta=0.0_dp, &
                            matrix_c=fm_mat_Gamma_3(ispin))

         ! Remove extra factor from S after the multiplication (to return to the original matrix)
         CALL remove_scaling_factor_rpa(fm_mat_S(ispin), axk_work%virtual(ispin), eig(:, ispin), axk_work%homo(ispin), omega)
      END DO

      CALL axk_work%axk_exchange%create(qs_env, mat_munu, axk_work%para_env_sub, axk_work%my_group_L_size, axk_work%eps_filter)
      axk_work%axk_exchange%e_axk = 0.0_dp

      ALLOCATE (dbcsr_Gamma_inu_P(nspins))
      CALL dbcsr_allocate_matrix_set(dbcsr_Gamma_munu_P, nspins)
      DO ispin = 1, nspins
         ALLOCATE (dbcsr_Gamma_munu_P(ispin)%matrix)
         CALL dbcsr_create(dbcsr_Gamma_munu_P(ispin)%matrix, template=mat_munu%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(dbcsr_Gamma_munu_P(ispin)%matrix, mat_munu%matrix)

         CALL dbcsr_create(dbcsr_Gamma_inu_P(ispin), template=axk_work%mo_coeff_o(ispin))
         CALL dbcsr_copy(dbcsr_Gamma_inu_P(ispin), axk_work%mo_coeff_o(ispin))
         CALL dbcsr_set(dbcsr_Gamma_inu_P(ispin), 0.0_dp)
      END DO

      ALLOCATE (dbcsr_Gamma_3(nspins))
      DO ispin = 1, nspins
         ! Copy Gamma_ia_P^3 to dbcsr matrix set
      CALL gamma_fm_to_dbcsr(fm_mat_Gamma_3(ispin), dbcsr_Gamma_3(ispin)%matrix_set, axk_work%para_env_RPA, axk_work%para_env_sub, &
                                axk_work%homo(ispin), axk_work%virtual(ispin), axk_work%mo_coeff_o(ispin), &
                                axk_work%ngroup, axk_work%my_group_L_start, &
                                axk_work%my_group_L_end, axk_work%my_group_L_size)
      END DO

      DO L_counter = 1, axk_work%my_group_L_size
         DO ispin = 1, nspins
            ! Do dbcsr multiplication: transform the virtual index
            CALL dbcsr_multiply("N", "T", 1.0_dp, axk_work%mo_coeff_v(ispin), dbcsr_Gamma_3(ispin)%matrix_set(L_counter), &
                                0.0_dp, dbcsr_Gamma_inu_P(ispin), filter_eps=axk_work%eps_filter)

            CALL dbcsr_release(dbcsr_Gamma_3(ispin)%matrix_set(L_counter))

            !Do dbcsr multiplication: transform the occupied index
            CALL dbcsr_multiply("N", "T", 1.0_dp, dbcsr_Gamma_inu_P(ispin), axk_work%mo_coeff_o(ispin), &
                                0.0_dp, dbcsr_Gamma_munu_P(ispin)%matrix, filter_eps=axk_work%eps_filter)
         END DO

         ! Contract DM with exchange integrals
         CALL axk_work%axk_exchange%integrate(dbcsr_Gamma_munu_P)
      END DO

      DO ispin = 1, nspins
         CALL dbcsr_release(dbcsr_Gamma_inu_P(ispin))
         CALL dbcsr_release(dbcsr_Gamma_munu_P(ispin)%matrix)
         DEALLOCATE (dbcsr_Gamma_3(ispin)%matrix_set, dbcsr_Gamma_munu_P(ispin)%matrix)
      END DO
      DEALLOCATE (dbcsr_Gamma_inu_P, dbcsr_Gamma_munu_P, dbcsr_Gamma_3)

      CALL axk_work%axk_exchange%release(e_axk_corr)

      ! Print AXK correlation energy to the file
      IF (axk_work%unit_nr > 0) WRITE (axk_work%unit_nr, '(T3,A,T68,F25.14,A4)') 'AXK correlation energy for a quadrature point:', &
         e_axk_corr, ' a.u.'

      CALL timestop(handle)

   END SUBROUTINE rpa_axk_work_compute

! **************************************************************************************************
!> \brief Contract RPA-AXK density matrix with HF exchange integrals and evaluate the correction
!> \param axk_exchange ...
!> \param dbcsr_Gamma_munu_P ...  AXK density matrix in AO basis to be contracted
!> \author Vladimir Rybkin, 08/2016
! **************************************************************************************************
   SUBROUTINE integrate_exchange(axk_exchange, dbcsr_Gamma_munu_P)
      CLASS(rpa_axk_exchange_type), INTENT(INOUT) :: axk_exchange
      TYPE(dbcsr_p_type), DIMENSION(:), TARGET          :: dbcsr_Gamma_munu_P

      CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_exchange'

      INTEGER                                            :: handle, ispin
      REAL(KIND=dp)                                      :: e_axk_P

      CALL timeset(routineN, handle)

      CALL dbcsr_set(axk_exchange%mat_hfx(1)%matrix, 0.0_dp)

      CALL tddft_hfx_matrix(axk_exchange%mat_hfx, dbcsr_Gamma_munu_P, axk_exchange%qs_env, .FALSE., &
                            axk_exchange%my_recalc_hfx_integrals, &
                            axk_exchange%hfx_sections, axk_exchange%x_data, axk_exchange%para_env)

      axk_exchange%my_recalc_hfx_integrals = .FALSE.
      DO ispin = 1, SIZE(dbcsr_Gamma_munu_P)
         CALL dbcsr_dot(axk_exchange%mat_hfx(ispin)%matrix, dbcsr_Gamma_munu_P(ispin)%matrix, e_axk_p)
         axk_exchange%e_axk = axk_exchange%e_axk + e_axk_P
      END DO

      CALL timestop(handle)

   END SUBROUTINE integrate_exchange

! **************************************************************************************************
!> \brief ... Initializes x_data on a subgroup
!> \param axk_exchange ...
!> \param qs_env ...
!> \param mat_munu ...
!> \param para_env_sub ...
!> \param my_group_L_size ...
!> \param eps_filter ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE hfx_create_subgroup(axk_exchange, qs_env, mat_munu, para_env_sub, my_group_L_size, eps_filter)
      CLASS(rpa_axk_exchange_type), INTENT(INOUT) :: axk_exchange
      TYPE(dbcsr_p_type), INTENT(IN) :: mat_munu
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

      axk_exchange%qs_env => qs_env
      axk_exchange%para_env => para_env_sub
      axk_exchange%eps_filter = eps_filter
      axk_exchange%my_group_L_size = my_group_L_size

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
      axk_exchange%hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")
      CALL section_vals_get(axk_exchange%hfx_sections, explicit=do_hfx, n_repetition=axk_exchange%n_rep_hf)
      CALL get_qs_env(qs_env, dft_control=dft_control)

      IF (do_hfx) THEN
         ! Retrieve particle_set and atomic_kind_set
         CALL hfx_create(axk_exchange%x_data, para_env_sub, axk_exchange%hfx_sections, atomic_kind_set, &
                         qs_kind_set, particle_set, dft_control, my_cell, orb_basis='ORB', &
                         nelectron_total=nelectron_total)
      END IF

      axk_exchange%my_recalc_hfx_integrals = .TRUE.

      ! create a working matrix
      DO ispin = 1, dft_control%nspins
         CALL dbcsr_create(axk_exchange%work_ao, template=mat_munu%matrix)
      END DO

      CALL dbcsr_allocate_matrix_set(axk_exchange%mat_hfx, dft_control%nspins)
      DO ispin = 1, dft_control%nspins
         ALLOCATE (axk_exchange%mat_hfx(ispin)%matrix)
         CALL dbcsr_init_p(axk_exchange%mat_hfx(ispin)%matrix)
         CALL dbcsr_create(axk_exchange%mat_hfx(ispin)%matrix, template=mat_munu%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(axk_exchange%mat_hfx(ispin)%matrix, mat_munu%matrix)
      END DO

      CALL timestop(handle)

   END SUBROUTINE hfx_create_subgroup

! **************************************************************************************************
!> \brief ...
!> \param axk_exchange ...
!> \param e_axk ...
! **************************************************************************************************
   SUBROUTINE hfx_release_subgroup(axk_exchange, e_axk)
      CLASS(rpa_axk_exchange_type), INTENT(INOUT) :: axk_exchange
      REAL(KIND=dp), INTENT(OUT) :: e_axk

      INTEGER :: ispin

      NULLIFY (axk_exchange%qs_env, axk_exchange%para_env, axk_exchange%hfx_sections)

      e_axk = axk_exchange%e_axk
      axk_exchange%e_axk = 0.0_dp

      CALL hfx_release(axk_exchange%x_data)
      NULLIFY (axk_exchange%x_data)

      DO ispin = 1, SIZE(axk_exchange%mat_hfx, 1)
         CALL dbcsr_release(axk_exchange%mat_hfx(ispin)%matrix)
         CALL dbcsr_release(axk_exchange%mat_hfx(ispin)%matrix)
         DEALLOCATE (axk_exchange%mat_hfx(ispin)%matrix)
      END DO
      DEALLOCATE (axk_exchange%mat_hfx)
      NULLIFY (axk_exchange%mat_hfx)

      CALL dbcsr_release(axk_exchange%work_ao)
   END SUBROUTINE hfx_release_subgroup

END MODULE rpa_axk
