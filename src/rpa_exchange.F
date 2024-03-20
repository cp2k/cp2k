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
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_get_info, dbcsr_init_p, dbcsr_multiply, dbcsr_p_type, &
        dbcsr_release, dbcsr_set, dbcsr_trace, dbcsr_type, dbcsr_type_no_symmetry
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
   USE group_dist_types,                ONLY: create_group_dist,&
                                              get_group_dist,&
                                              group_dist_d1_type,&
                                              group_dist_proc,&
                                              maxsize,&
                                              release_group_dist
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE hfx_types,                       ONLY: hfx_create,&
                                              hfx_release,&
                                              hfx_type
   USE input_constants,                 ONLY: rpa_exchange_axk,&
                                              rpa_exchange_none,&
                                              rpa_exchange_sosex
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE local_gemm_api,                  ONLY: LOCAL_GEMM_PU_GPU
   USE mathconstants,                   ONLY: sqrthalf
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_proc_null
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
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_exchange'

   PUBLIC :: rpa_exchange_work_type, rpa_exchange_needed_mem

   TYPE rpa_exchange_env_type
      PRIVATE
      TYPE(qs_environment_type), POINTER             :: qs_env => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER      :: mat_hfx => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER      :: dbcsr_Gamma_munu_P => NULL()
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)    :: dbcsr_Gamma_inu_P
      ! Workaround GCC 8
      TYPE(dbcsr_type), DIMENSION(:), POINTER :: mo_coeff_o => NULL()
      TYPE(dbcsr_type), DIMENSION(:), POINTER :: mo_coeff_v => NULL()
      TYPE(dbcsr_type)                               :: work_ao
      TYPE(hfx_type), DIMENSION(:, :), POINTER       :: x_data => NULL()
      TYPE(mp_para_env_type), POINTER                :: para_env => NULL()
      TYPE(section_vals_type), POINTER               :: hfx_sections => NULL()
      LOGICAL :: my_recalc_hfx_integrals = .FALSE.
      REAL(KIND=dp) :: eps_filter = 0.0_dp
      TYPE(cp_fm_struct_p_type), DIMENSION(:), ALLOCATABLE :: struct_Gamma
   CONTAINS
      PROCEDURE, PASS(exchange_env), NON_OVERRIDABLE :: create => hfx_create_subgroup
      !PROCEDURE, PASS(exchange_env), NON_OVERRIDABLE :: integrate => integrate_exchange
      PROCEDURE, PASS(exchange_env), NON_OVERRIDABLE :: release => hfx_release_subgroup
   END TYPE

   TYPE dbcsr_matrix_p_set
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:) :: matrix_set
   END TYPE

   TYPE rpa_exchange_work_type
      PRIVATE
      INTEGER :: exchange_correction = rpa_exchange_none
      TYPE(rpa_exchange_env_type) :: exchange_env
      INTEGER, DIMENSION(:), ALLOCATABLE :: homo, virtual, dimen_ia
      TYPE(group_dist_d1_type) :: aux_func_dist = group_dist_d1_type()
      INTEGER, DIMENSION(:), ALLOCATABLE :: aux2send
      INTEGER :: dimen_RI = 0
      INTEGER :: block_size = 0
      INTEGER :: color_sub = 0
      INTEGER :: ngroup = 0
      TYPE(cp_fm_type) :: fm_mat_Q_tmp = cp_fm_type()
      TYPE(cp_fm_type) :: fm_mat_R_half_gemm = cp_fm_type()
      TYPE(cp_fm_type) :: fm_mat_U = cp_fm_type()
      TYPE(mp_para_env_type), POINTER :: para_env_sub => NULL()
   CONTAINS
      PROCEDURE, PUBLIC, PASS(exchange_work), NON_OVERRIDABLE :: create => rpa_exchange_work_create
      PROCEDURE, PUBLIC, PASS(exchange_work), NON_OVERRIDABLE :: compute => rpa_exchange_work_compute
      PROCEDURE, PUBLIC, PASS(exchange_work), NON_OVERRIDABLE :: release => rpa_exchange_work_release
      PROCEDURE, PRIVATE, PASS(exchange_work), NON_OVERRIDABLE :: redistribute_into_subgroups
      PROCEDURE, PRIVATE, PASS(exchange_work), NON_OVERRIDABLE :: compute_fm => rpa_exchange_work_compute_fm
      PROCEDURE, PRIVATE, PASS(exchange_work), NON_OVERRIDABLE :: compute_hfx => rpa_exchange_work_compute_hfx
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param para_env ...
!> \param mem_per_rank ...
!> \param mem_per_repl ...
! **************************************************************************************************
   SUBROUTINE rpa_exchange_needed_mem(mp2_env, homo, virtual, dimen_RI, para_env, mem_per_rank, mem_per_repl)
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      INTEGER, INTENT(IN)                                :: dimen_RI
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp), INTENT(INOUT)                       :: mem_per_rank, mem_per_repl

      INTEGER                                            :: block_size

      ! We need the block size and if it is unknown, an upper bound
      block_size = mp2_env%ri_rpa%exchange_block_size
      IF (block_size <= 0) block_size = MAX(1, (dimen_RI + para_env%num_pe - 1)/para_env%num_pe)

      ! storage of product matrix (upper bound only as it depends on the square of the potential still unknown block size)
      mem_per_rank = mem_per_rank + REAL(MAXVAL(homo), KIND=dp)**2*block_size**2*8.0_dp/(1024_dp**2)

      ! work arrays R (2x) and U, copies of Gamma (2x), communication buffer (as expensive as Gamma)
      mem_per_repl = mem_per_repl + 3.0_dp*dimen_RI*dimen_RI*8.0_dp/(1024_dp**2) &
                     + 3.0_dp*MAXVAL(homo*virtual)*dimen_RI*8.0_dp/(1024_dp**2)
   END SUBROUTINE rpa_exchange_needed_mem

! **************************************************************************************************
!> \brief ...
!> \param exchange_work ...
!> \param qs_env ...
!> \param para_env_sub ...
!> \param mat_munu ...
!> \param dimen_RI ...
!> \param fm_mat_S ...
!> \param fm_mat_Q ...
!> \param fm_mat_Q_gemm ...
!> \param homo ...
!> \param virtual ...
! **************************************************************************************************
   SUBROUTINE rpa_exchange_work_create(exchange_work, qs_env, para_env_sub, mat_munu, dimen_RI, &
                                       fm_mat_S, fm_mat_Q, fm_mat_Q_gemm, homo, virtual)
      CLASS(rpa_exchange_work_type), INTENT(INOUT) :: exchange_work
      TYPE(qs_environment_type), POINTER :: qs_env
      TYPE(mp_para_env_type), POINTER, INTENT(IN) :: para_env_sub
      TYPE(dbcsr_p_type), INTENT(IN) :: mat_munu
      INTEGER, INTENT(IN) :: dimen_RI
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN) :: fm_mat_S
      TYPE(cp_fm_type), INTENT(IN) :: fm_mat_Q, fm_mat_Q_gemm
      INTEGER, DIMENSION(SIZE(fm_mat_S)), INTENT(IN) :: homo, virtual

      INTEGER :: nspins, aux_global, aux_local, my_process_row, proc, ispin
      INTEGER, DIMENSION(:), POINTER :: row_indices, aux_distribution_fm
      TYPE(cp_blacs_env_type), POINTER :: context

      exchange_work%exchange_correction = qs_env%mp2_env%ri_rpa%exchange_correction

      IF (exchange_work%exchange_correction == rpa_exchange_none) RETURN

      ASSOCIATE (para_env => fm_mat_S(1)%matrix_struct%para_env)
         exchange_work%para_env_sub => para_env_sub
         exchange_work%ngroup = para_env%num_pe/para_env_sub%num_pe
         exchange_work%color_sub = para_env%mepos/para_env_sub%num_pe
      END ASSOCIATE

      CALL cp_fm_get_info(fm_mat_S(1), row_indices=row_indices, nrow_locals=aux_distribution_fm, context=context)
      CALL context%get(my_process_row=my_process_row)

      CALL create_group_dist(exchange_work%aux_func_dist, exchange_work%ngroup, dimen_RI)
      ALLOCATE (exchange_work%aux2send(0:exchange_work%ngroup - 1))
      exchange_work%aux2send = 0
      DO aux_local = 1, aux_distribution_fm(my_process_row)
         aux_global = row_indices(aux_local)
         proc = group_dist_proc(exchange_work%aux_func_dist, aux_global)
         exchange_work%aux2send(proc) = exchange_work%aux2send(proc) + 1
      END DO

      nspins = SIZE(fm_mat_S)

      ALLOCATE (exchange_work%homo(nspins), exchange_work%virtual(nspins), exchange_work%dimen_ia(nspins))
      exchange_work%homo(:) = homo
      exchange_work%virtual(:) = virtual
      exchange_work%dimen_ia(:) = homo*virtual
      exchange_work%dimen_RI = dimen_RI

      exchange_work%block_size = qs_env%mp2_env%ri_rpa%exchange_block_size
      IF (exchange_work%block_size <= 0) exchange_work%block_size = dimen_RI

      CALL cp_fm_create(exchange_work%fm_mat_U, fm_mat_Q%matrix_struct, name="fm_mat_U")
      CALL cp_fm_create(exchange_work%fm_mat_Q_tmp, fm_mat_Q%matrix_struct, name="fm_mat_Q_tmp")
      CALL cp_fm_create(exchange_work%fm_mat_R_half_gemm, fm_mat_Q_gemm%matrix_struct)

      IF (qs_env%mp2_env%ri_rpa%use_hfx_implementation) THEN
         CALL exchange_work%exchange_env%create(qs_env, mat_munu%matrix, para_env_sub, fm_mat_S)
      END IF

      IF (ASSOCIATED(qs_env%mp2_env%ri_rpa%mo_coeff_o)) THEN
         DO ispin = 1, SIZE(qs_env%mp2_env%ri_rpa%mo_coeff_o)
            CALL dbcsr_release(qs_env%mp2_env%ri_rpa%mo_coeff_o(ispin))
         END DO
         DEALLOCATE (qs_env%mp2_env%ri_rpa%mo_coeff_o)
      END IF

      IF (ASSOCIATED(qs_env%mp2_env%ri_rpa%mo_coeff_v)) THEN
         DO ispin = 1, SIZE(qs_env%mp2_env%ri_rpa%mo_coeff_v)
            CALL dbcsr_release(qs_env%mp2_env%ri_rpa%mo_coeff_v(ispin))
         END DO
         DEALLOCATE (qs_env%mp2_env%ri_rpa%mo_coeff_v)
      END IF
   END SUBROUTINE

! **************************************************************************************************
!> \brief ... Initializes x_data on a subgroup
!> \param exchange_env ...
!> \param qs_env ...
!> \param mat_munu ...
!> \param para_env_sub ...
!> \param fm_mat_S ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE hfx_create_subgroup(exchange_env, qs_env, mat_munu, para_env_sub, fm_mat_S)
      CLASS(rpa_exchange_env_type), INTENT(INOUT) :: exchange_env
      TYPE(dbcsr_type), INTENT(IN) :: mat_munu
      TYPE(qs_environment_type), POINTER   :: qs_env
      TYPE(mp_para_env_type), POINTER, INTENT(IN)            :: para_env_sub
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN) :: fm_mat_S

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_create_subgroup'

      INTEGER                                            :: handle, nelectron_total, ispin, &
                                                            number_of_aos, nspins, dimen_RI, dimen_ia
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: my_cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      exchange_env%mo_coeff_o => qs_env%mp2_env%ri_rpa%mo_coeff_o
      exchange_env%mo_coeff_v => qs_env%mp2_env%ri_rpa%mo_coeff_v
      NULLIFY (qs_env%mp2_env%ri_rpa%mo_coeff_o, qs_env%mp2_env%ri_rpa%mo_coeff_v)

      nspins = SIZE(exchange_env%mo_coeff_o)

      exchange_env%qs_env => qs_env
      exchange_env%para_env => para_env_sub
      exchange_env%eps_filter = qs_env%mp2_env%mp2_gpw%eps_filter

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

      exchange_env%hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")
      CALL get_qs_env(qs_env, dft_control=dft_control)

      ! Retrieve particle_set and atomic_kind_set
      CALL hfx_create(exchange_env%x_data, para_env_sub, exchange_env%hfx_sections, atomic_kind_set, &
                      qs_kind_set, particle_set, dft_control, my_cell, orb_basis='ORB', &
                      nelectron_total=nelectron_total)

      exchange_env%my_recalc_hfx_integrals = .TRUE.

      CALL dbcsr_allocate_matrix_set(exchange_env%mat_hfx, nspins)
      DO ispin = 1, nspins
         ALLOCATE (exchange_env%mat_hfx(ispin)%matrix)
         CALL dbcsr_init_p(exchange_env%mat_hfx(ispin)%matrix)
         CALL dbcsr_create(exchange_env%mat_hfx(ispin)%matrix, template=mat_munu, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(exchange_env%mat_hfx(ispin)%matrix, mat_munu)
      END DO

      CALL dbcsr_get_info(mat_munu, nfullcols_total=number_of_aos)

      CALL dbcsr_create(exchange_env%work_ao, template=mat_munu, &
                        matrix_type=dbcsr_type_no_symmetry)

      ALLOCATE (exchange_env%dbcsr_Gamma_inu_P(nspins))
      CALL dbcsr_allocate_matrix_set(exchange_env%dbcsr_Gamma_munu_P, nspins)
      DO ispin = 1, nspins
         ALLOCATE (exchange_env%dbcsr_Gamma_munu_P(ispin)%matrix)
         CALL dbcsr_create(exchange_env%dbcsr_Gamma_munu_P(ispin)%matrix, template=mat_munu, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(exchange_env%dbcsr_Gamma_munu_P(ispin)%matrix, mat_munu)
         CALL dbcsr_set(exchange_env%dbcsr_Gamma_munu_P(ispin)%matrix, 0.0_dp)

         CALL dbcsr_create(exchange_env%dbcsr_Gamma_inu_P(ispin), template=exchange_env%mo_coeff_o(ispin))
         CALL dbcsr_copy(exchange_env%dbcsr_Gamma_inu_P(ispin), exchange_env%mo_coeff_o(ispin))
         CALL dbcsr_set(exchange_env%dbcsr_Gamma_inu_P(ispin), 0.0_dp)
      END DO

      ALLOCATE (exchange_env%struct_Gamma(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_get_info(fm_mat_S(ispin), nrow_global=dimen_RI, ncol_global=dimen_ia)
         CALL cp_fm_struct_create(exchange_env%struct_Gamma(ispin)%struct, template_fmstruct=fm_mat_S(ispin)%matrix_struct, &
                                  nrow_global=dimen_ia, ncol_global=dimen_RI)
      END DO

      CALL timestop(handle)

   END SUBROUTINE hfx_create_subgroup

! **************************************************************************************************
!> \brief ...
!> \param exchange_work ...
! **************************************************************************************************
   SUBROUTINE rpa_exchange_work_release(exchange_work)
      CLASS(rpa_exchange_work_type), INTENT(INOUT) :: exchange_work

      IF (ALLOCATED(exchange_work%homo)) DEALLOCATE (exchange_work%homo)
      IF (ALLOCATED(exchange_work%virtual)) DEALLOCATE (exchange_work%virtual)
      IF (ALLOCATED(exchange_work%dimen_ia)) DEALLOCATE (exchange_work%dimen_ia)
      NULLIFY (exchange_work%para_env_sub)
      CALL release_group_dist(exchange_work%aux_func_dist)
      IF (ALLOCATED(exchange_work%aux2send)) DEALLOCATE (exchange_work%aux2send)
      CALL cp_fm_release(exchange_work%fm_mat_Q_tmp)
      CALL cp_fm_release(exchange_work%fm_mat_U)
      CALL cp_fm_release(exchange_work%fm_mat_R_half_gemm)

      CALL exchange_work%exchange_env%release()
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param exchange_env ...
! **************************************************************************************************
   SUBROUTINE hfx_release_subgroup(exchange_env)
      CLASS(rpa_exchange_env_type), INTENT(INOUT) :: exchange_env

      INTEGER :: ispin

      NULLIFY (exchange_env%para_env, exchange_env%hfx_sections)

      IF (ASSOCIATED(exchange_env%x_data)) THEN
         CALL hfx_release(exchange_env%x_data)
         NULLIFY (exchange_env%x_data)
      END IF

      CALL dbcsr_release(exchange_env%work_ao)

      IF (ASSOCIATED(exchange_env%dbcsr_Gamma_munu_P)) THEN
         DO ispin = 1, SIZE(exchange_env%mat_hfx, 1)
            CALL dbcsr_release(exchange_env%dbcsr_Gamma_munu_P(ispin)%matrix)
            CALL dbcsr_release(exchange_env%mat_hfx(ispin)%matrix)
            CALL dbcsr_release(exchange_env%dbcsr_Gamma_inu_P(ispin))
            CALL dbcsr_release(exchange_env%mo_coeff_o(ispin))
            CALL dbcsr_release(exchange_env%mo_coeff_v(ispin))
            DEALLOCATE (exchange_env%mat_hfx(ispin)%matrix)
            DEALLOCATE (exchange_env%dbcsr_Gamma_munu_P(ispin)%matrix)
         END DO
         DEALLOCATE (exchange_env%mat_hfx, exchange_env%dbcsr_Gamma_munu_P)
         DEALLOCATE (exchange_env%dbcsr_Gamma_inu_P, exchange_env%mo_coeff_o, exchange_env%mo_coeff_v)
         NULLIFY (exchange_env%mat_hfx, exchange_env%dbcsr_Gamma_munu_P)
      END IF
      IF (ALLOCATED(exchange_env%struct_Gamma)) THEN
      DO ispin = 1, SIZE(exchange_env%struct_Gamma)
         CALL cp_fm_struct_release(exchange_env%struct_Gamma(ispin)%struct)
      END DO
      DEALLOCATE (exchange_env%struct_Gamma)
      END IF
   END SUBROUTINE hfx_release_subgroup

! **************************************************************************************************
!> \brief Main driver for RPA-exchange energies
!> \param exchange_work ...
!> \param fm_mat_Q ...
!> \param eig ...
!> \param fm_mat_S ...
!> \param omega ...
!> \param e_exchange_corr exchange energy correction for a quadrature point
!> \param mp2_env ...
!> \author Vladimir Rybkin, 07/2016
! **************************************************************************************************
   SUBROUTINE rpa_exchange_work_compute(exchange_work, fm_mat_Q, eig, fm_mat_S, omega, &
                                        e_exchange_corr, mp2_env)
      CLASS(rpa_exchange_work_type), INTENT(INOUT) :: exchange_work
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_Q
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eig
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT)         :: fm_mat_S
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), INTENT(INOUT)                       :: e_exchange_corr
      TYPE(mp2_type), INTENT(INOUT) :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rpa_exchange_work_compute'
      REAL(KIND=dp), PARAMETER                           :: thresh = 0.0000001_dp

      INTEGER :: handle, nspins, dimen_RI, iiB
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval

      IF (exchange_work%exchange_correction == rpa_exchange_none) RETURN

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(fm_mat_Q, ncol_global=dimen_RI)

      nspins = SIZE(fm_mat_S)

      ! Eigenvalues
      ALLOCATE (eigenval(dimen_RI))
      eigenval = 0.0_dp

      CALL cp_fm_set_all(matrix=exchange_work%fm_mat_Q_tmp, alpha=0.0_dp)
      CALL cp_fm_set_all(matrix=exchange_work%fm_mat_U, alpha=0.0_dp)

      ! Copy Q to Q_tmp
      CALL cp_fm_to_fm(fm_mat_Q, exchange_work%fm_mat_Q_tmp)
      ! Diagonalize Q
      CALL choose_eigv_solver(exchange_work%fm_mat_Q_tmp, exchange_work%fm_mat_U, eigenval)

      ! Calculate diagonal matrix for R_half

      ! Manipulate eigenvalues to get diagonal matrix
      IF (exchange_work%exchange_correction == rpa_exchange_axk) THEN
         DO iib = 1, dimen_RI
            IF (ABS(eigenval(iib)) .GE. thresh) THEN
               eigenval(iib) = &
                  SQRT((1.0_dp/(eigenval(iib)**2))*LOG(1.0_dp + eigenval(iib)) &
                       - 1.0_dp/(eigenval(iib)*(eigenval(iib) + 1.0_dp)))
            ELSE
               eigenval(iib) = sqrthalf
            END IF
         END DO
      ELSE IF (exchange_work%exchange_correction == rpa_exchange_sosex) THEN
         DO iib = 1, dimen_RI
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

      CALL cp_fm_to_fm_submat_general(exchange_work%fm_mat_U, exchange_work%fm_mat_R_half_gemm, dimen_RI, &
                                      dimen_RI, 1, 1, 1, 1, exchange_work%fm_mat_U%matrix_struct%context)

      IF (mp2_env%ri_rpa%use_hfx_implementation) THEN
         CALL exchange_work%compute_hfx(fm_mat_S, eig, omega, e_exchange_corr)
      ELSE
         CALL exchange_work%compute_fm(fm_mat_S, eig, omega, e_exchange_corr, mp2_env)
      END IF

      CALL timestop(handle)

   END SUBROUTINE rpa_exchange_work_compute

! **************************************************************************************************
!> \brief Main driver for RPA-exchange energies
!> \param exchange_work ...
!> \param fm_mat_S ...
!> \param eig ...
!> \param omega ...
!> \param e_exchange_corr exchange energy correction for a quadrature point
!> \param mp2_env ...
!> \author Frederick Stein, May-June 2024
! **************************************************************************************************
   SUBROUTINE rpa_exchange_work_compute_fm(exchange_work, fm_mat_S, eig, omega, &
                                           e_exchange_corr, mp2_env)
      CLASS(rpa_exchange_work_type), INTENT(INOUT) :: exchange_work
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_S
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eig
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), INTENT(INOUT)                       :: e_exchange_corr
      TYPE(mp2_type), INTENT(INOUT) :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rpa_exchange_work_compute_fm'

      INTEGER :: handle, ispin, nspins, P, Q, L_size_Gamma, hom, virt, i, &
                 send_proc, recv_proc, recv_size, max_aux_size, proc_shift, dimen_ia, &
                 block_size, P_start, P_end, P_size, Q_start, Q_size, Q_end, handle2, my_aux_size, my_virt
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), TARGET :: mat_Gamma_3_3D
      REAL(KIND=dp), POINTER, DIMENSION(:), CONTIGUOUS :: mat_Gamma_3_1D
      REAL(KIND=dp), POINTER, DIMENSION(:, :), CONTIGUOUS :: mat_Gamma_3_2D
      REAL(KIND=dp), ALLOCATABLE, TARGET, DIMENSION(:) :: recv_buffer_1D
      REAL(KIND=dp), POINTER, DIMENSION(:, :), CONTIGUOUS :: recv_buffer_2D
      REAL(KIND=dp), POINTER, DIMENSION(:, :, :), CONTIGUOUS :: recv_buffer_3D
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :) :: mat_B_iaP
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET :: product_matrix_1D
      REAL(KIND=dp), POINTER, DIMENSION(:, :), CONTIGUOUS :: product_matrix_2D
      REAL(KIND=dp), POINTER, DIMENSION(:, :, :, :), CONTIGUOUS :: product_matrix_4D
      TYPE(cp_fm_type)        :: fm_mat_Gamma_3
      TYPE(mp_para_env_type), POINTER :: para_env
      TYPE(group_dist_d1_type)                           :: virt_dist

      CALL timeset(routineN, handle)

      nspins = SIZE(fm_mat_S)

      CALL get_group_dist(exchange_work%aux_func_dist, exchange_work%color_sub, sizes=my_aux_size)

      e_exchange_corr = 0.0_dp
      max_aux_size = maxsize(exchange_work%aux_func_dist)

      ! local_gemm_ctx has a very large footprint the first time this routine is
      ! called.
      CALL mp2_env%local_gemm_ctx%create(LOCAL_GEMM_PU_GPU)
      CALL mp2_env%local_gemm_ctx%set_op_threshold_gpu(128*128*128*2)

      DO ispin = 1, nspins
         hom = exchange_work%homo(ispin)
         virt = exchange_work%virtual(ispin)
         dimen_ia = hom*virt
         IF (hom < 1 .OR. virt < 1) CYCLE

         CALL cp_fm_get_info(fm_mat_S(ispin), para_env=para_env)

         CALL cp_fm_create(fm_mat_Gamma_3, fm_mat_S(ispin)%matrix_struct)
         CALL cp_fm_set_all(matrix=fm_mat_Gamma_3, alpha=0.0_dp)

         ! Update G with a new value of Omega: in practice, it is G*S

         ! Scale fm_work_iaP
         CALL calc_fm_mat_S_rpa(fm_mat_S(ispin), .TRUE., virt, eig(:, ispin), &
                                omega, 0.0_dp)

         ! Calculate Gamma_3: Gamma_3 = G*S*R^(1/2) = G*S*R^(1/2)
         CALL parallel_gemm(transa="T", transb="N", m=exchange_work%dimen_RI, n=dimen_ia, k=exchange_work%dimen_RI, alpha=1.0_dp, &
                            matrix_a=exchange_work%fm_mat_R_half_gemm, matrix_b=fm_mat_S(ispin), beta=0.0_dp, &
                            matrix_c=fm_mat_Gamma_3)

         CALL create_group_dist(virt_dist, exchange_work%para_env_sub%num_pe, virt)

         ! Remove extra factor from S after the multiplication (to return to the original matrix)
         CALL remove_scaling_factor_rpa(fm_mat_S(ispin), virt, eig(:, ispin), hom, omega)

         CALL exchange_work%redistribute_into_subgroups(fm_mat_Gamma_3, mat_Gamma_3_3D, ispin, virt_dist)
         CALL cp_fm_release(fm_mat_Gamma_3)

         ! We need only the pure matrix
         CALL remove_scaling_factor_rpa(fm_mat_S(ispin), virt, eig(:, ispin), hom, omega)

         ! Reorder matrix from (P, i*a) -> (a, i, P) with P being distributed within subgroups
         CALL exchange_work%redistribute_into_subgroups(fm_mat_S(ispin), mat_B_iaP, ispin, virt_dist)

         ! Return to the original tensor
         CALL calc_fm_mat_S_rpa(fm_mat_S(ispin), .TRUE., virt, eig(:, ispin), omega, 0.0_dp)

         L_size_Gamma = SIZE(mat_Gamma_3_3D, 3)
         my_virt = SIZE(mat_Gamma_3_3D, 1)
         block_size = exchange_work%block_size

         mat_Gamma_3_1D(1:INT(my_virt, KIND=int_8)*hom*my_aux_size) => mat_Gamma_3_3D(:, :, 1:my_aux_size)
         mat_Gamma_3_2D(1:my_virt, 1:hom*my_aux_size) => mat_Gamma_3_1D(1:INT(my_virt, KIND=int_8)*hom*my_aux_size)

         ALLOCATE (product_matrix_1D(INT(hom*MIN(block_size, L_size_gamma), KIND=int_8)* &
                                     INT(hom*MIN(block_size, max_aux_size), KIND=int_8)))
         ALLOCATE (recv_buffer_1D(INT(virt, KIND=int_8)*hom*max_aux_size))
         recv_buffer_2D(1:my_virt, 1:hom*max_aux_size) => recv_buffer_1D(1:INT(virt, KIND=int_8)*hom*max_aux_size)
         recv_buffer_3D(1:my_virt, 1:hom, 1:max_aux_size) => recv_buffer_1D(1:INT(virt, KIND=int_8)*hom*max_aux_size)
         DO proc_shift = 0, para_env%num_pe - 1, exchange_work%para_env_sub%num_pe
            send_proc = MODULO(para_env%mepos + proc_shift, para_env%num_pe)
            recv_proc = MODULO(para_env%mepos - proc_shift, para_env%num_pe)

            CALL get_group_dist(exchange_work%aux_func_dist, recv_proc/exchange_work%para_env_sub%num_pe, sizes=recv_size)

            IF (recv_size == 0) recv_proc = mp_proc_null

            CALL para_env%sendrecv(mat_B_iaP, send_proc, recv_buffer_3D(:, :, 1:recv_size), recv_proc)

            IF (recv_size == 0) CYCLE

            DO P_start = 1, L_size_Gamma, block_size
               P_end = MIN(L_size_Gamma, P_start + block_size - 1)
               P_size = P_end - P_start + 1
               DO Q_start = 1, recv_size, block_size
                  Q_end = MIN(recv_size, Q_start + block_size - 1)
                  Q_size = Q_end - Q_start + 1

                  ! Reassign product_matrix pointers to enforce contiguity of target array
                  product_matrix_2D(1:hom*P_size, 1:hom*Q_size) => &
                     product_matrix_1D(1:INT(hom*P_size, KIND=int_8)*INT(hom*Q_size, KIND=int_8))
                  product_matrix_4D(1:hom, 1:P_size, 1:hom, 1:Q_size) => &
                     product_matrix_1D(1:INT(hom*P_size, KIND=int_8)*INT(hom*Q_size, KIND=int_8))

                  CALL timeset(routineN//"_gemm", handle2)
                  CALL mp2_env%local_gemm_ctx%gemm("T", "N", hom*P_size, hom*Q_size, my_virt, 1.0_dp, &
                                                   mat_Gamma_3_2D(:, hom*(P_start - 1) + 1:hom*P_end), my_virt, &
                                                   recv_buffer_2D(:, hom*(Q_start - 1) + 1:hom*Q_end), my_virt, &
                                                   0.0_dp, product_matrix_2D, hom*P_size)
                  CALL timestop(handle2)

                  CALL timeset(routineN//"_energy", handle2)
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(P_size, Q_size, hom, product_matrix_4D) &
!$OMP             COLLAPSE(3) REDUCTION(+: e_exchange_corr) PRIVATE(P, Q, i)
                  DO P = 1, P_size
                  DO Q = 1, Q_size
                  DO i = 1, hom
                     e_exchange_corr = e_exchange_corr + DOT_PRODUCT(product_matrix_4D(i, P, :, Q), product_matrix_4D(:, P, i, Q))
                  END DO
                  END DO
                  END DO
                  CALL timestop(handle2)
               END DO
            END DO
         END DO

         CALL release_group_dist(virt_dist)
         IF (ALLOCATED(mat_B_iaP)) DEALLOCATE (mat_B_iaP)
         IF (ALLOCATED(mat_Gamma_3_3D)) DEALLOCATE (mat_Gamma_3_3D)
         IF (ALLOCATED(product_matrix_1D)) DEALLOCATE (product_matrix_1D)
         IF (ALLOCATED(recv_buffer_1D)) DEALLOCATE (recv_buffer_1D)
      END DO

      CALL mp2_env%local_gemm_ctx%destroy()

      IF (nspins == 2) e_exchange_corr = e_exchange_corr*2.0_dp
      IF (nspins == 1) e_exchange_corr = e_exchange_corr*4.0_dp

      CALL timestop(handle)

   END SUBROUTINE rpa_exchange_work_compute_fm

! **************************************************************************************************
!> \brief Contract RPA-exchange density matrix with HF exchange integrals and evaluate the correction
!> \param exchange_work ...
!> \param fm_mat_S ...
!> \param eig ...
!> \param omega ...
!> \param e_exchange_corr ...
!> \author Vladimir Rybkin, 08/2016
! **************************************************************************************************
   SUBROUTINE rpa_exchange_work_compute_hfx(exchange_work, fm_mat_S, eig, omega, e_exchange_corr)
      CLASS(rpa_exchange_work_type), INTENT(INOUT) :: exchange_work
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT) :: fm_mat_S
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eig
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), INTENT(OUT) :: e_exchange_corr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rpa_exchange_work_compute_hfx'

      INTEGER                                            :: handle, ispin, my_aux_start, my_aux_end, &
                                                            my_aux_size, nspins, L_counter, dimen_ia, hom, virt
      REAL(KIND=dp)                                      :: e_exchange_P
      TYPE(dbcsr_matrix_p_set), DIMENSION(:), ALLOCATABLE          :: dbcsr_Gamma_3
      TYPE(cp_fm_type) :: fm_mat_Gamma_3
      TYPE(mp_para_env_type), POINTER :: para_env

      CALL timeset(routineN, handle)

      e_exchange_corr = 0.0_dp

      nspins = SIZE(fm_mat_S)

      CALL get_group_dist(exchange_work%aux_func_dist, exchange_work%color_sub, my_aux_start, my_aux_end, my_aux_size)

      ALLOCATE (dbcsr_Gamma_3(nspins))
      DO ispin = 1, nspins
         hom = exchange_work%homo(ispin)
         virt = exchange_work%virtual(ispin)
         dimen_ia = hom*virt
         IF (hom < 1 .OR. virt < 1) CYCLE

         CALL cp_fm_get_info(fm_mat_S(ispin), para_env=para_env)

         CALL cp_fm_create(fm_mat_Gamma_3, exchange_work%exchange_env%struct_Gamma(ispin)%struct)
         CALL cp_fm_set_all(matrix=fm_mat_Gamma_3, alpha=0.0_dp)

         ! Update G with a new value of Omega: in practice, it is G*S

         ! Scale fm_work_iaP
         CALL calc_fm_mat_S_rpa(fm_mat_S(ispin), .TRUE., virt, eig(:, ispin), &
                                omega, 0.0_dp)

         ! Calculate Gamma_3: Gamma_3 = G*S*R^(1/2) = G*S*R^(1/2)
         CALL parallel_gemm(transa="T", transb="N", m=dimen_ia, n=exchange_work%dimen_RI, &
                            k=exchange_work%dimen_RI, alpha=1.0_dp, &
                            matrix_a=fm_mat_S(ispin), matrix_b=exchange_work%fm_mat_R_half_gemm, beta=0.0_dp, &
                            matrix_c=fm_mat_Gamma_3)

         ! Remove extra factor from S after the multiplication (to return to the original matrix)
         CALL remove_scaling_factor_rpa(fm_mat_S(ispin), virt, eig(:, ispin), hom, omega)

         ! Copy Gamma_ia_P^3 to dbcsr matrix set
         CALL gamma_fm_to_dbcsr(fm_mat_Gamma_3, dbcsr_Gamma_3(ispin)%matrix_set, &
                                para_env, exchange_work%para_env_sub, hom, virt, &
                                exchange_work%exchange_env%mo_coeff_o(ispin), &
                                exchange_work%ngroup, my_aux_start, my_aux_end, my_aux_size)
      END DO

      DO L_counter = 1, my_aux_size
         DO ispin = 1, nspins
            ! Do dbcsr multiplication: transform the virtual index
            CALL dbcsr_multiply("N", "T", 1.0_dp, exchange_work%exchange_env%mo_coeff_v(ispin), &
                                dbcsr_Gamma_3(ispin)%matrix_set(L_counter), &
                                0.0_dp, exchange_work%exchange_env%dbcsr_Gamma_inu_P(ispin), &
                                filter_eps=exchange_work%exchange_env%eps_filter)

            CALL dbcsr_release(dbcsr_Gamma_3(ispin)%matrix_set(L_counter))

            ! Do dbcsr multiplication: transform the occupied index
            CALL dbcsr_multiply("N", "T", 0.5_dp, exchange_work%exchange_env%dbcsr_Gamma_inu_P(ispin), &
                                exchange_work%exchange_env%mo_coeff_o(ispin), &
                                0.0_dp, exchange_work%exchange_env%dbcsr_Gamma_munu_P(ispin)%matrix, &
                                filter_eps=exchange_work%exchange_env%eps_filter)
            CALL dbcsr_multiply("N", "T", 0.5_dp, exchange_work%exchange_env%mo_coeff_o(ispin), &
                                exchange_work%exchange_env%dbcsr_Gamma_inu_P(ispin), &
                                1.0_dp, exchange_work%exchange_env%dbcsr_Gamma_munu_P(ispin)%matrix, &
                                filter_eps=exchange_work%exchange_env%eps_filter)

            CALL dbcsr_set(exchange_work%exchange_env%mat_hfx(ispin)%matrix, 0.0_dp)
         END DO

         CALL tddft_hfx_matrix(exchange_work%exchange_env%mat_hfx, exchange_work%exchange_env%dbcsr_Gamma_munu_P, &
                               exchange_work%exchange_env%qs_env, .FALSE., &
                               exchange_work%exchange_env%my_recalc_hfx_integrals, &
                               exchange_work%exchange_env%hfx_sections, exchange_work%exchange_env%x_data, &
                               exchange_work%exchange_env%para_env)

         exchange_work%exchange_env%my_recalc_hfx_integrals = .FALSE.
         DO ispin = 1, nspins
            CALL dbcsr_multiply("N", "T", 1.0_dp, exchange_work%exchange_env%mat_hfx(ispin)%matrix, &
                                exchange_work%exchange_env%dbcsr_Gamma_munu_P(ispin)%matrix, &
                                0.0_dp, exchange_work%exchange_env%work_ao, filter_eps=exchange_work%exchange_env%eps_filter)
            CALL dbcsr_trace(exchange_work%exchange_env%work_ao, e_exchange_P)
            e_exchange_corr = e_exchange_corr - e_exchange_P
         END DO
      END DO

      IF (nspins == 2) e_exchange_corr = e_exchange_corr
      IF (nspins == 1) e_exchange_corr = e_exchange_corr*4.0_dp

      CALL timestop(handle)

   END SUBROUTINE rpa_exchange_work_compute_hfx

! **************************************************************************************************
!> \brief ...
!> \param exchange_work ...
!> \param fm_mat ...
!> \param mat ...
!> \param ispin ...
!> \param virt_dist ...
! **************************************************************************************************
   SUBROUTINE redistribute_into_subgroups(exchange_work, fm_mat, mat, ispin, virt_dist)
      CLASS(rpa_exchange_work_type), INTENT(IN) :: exchange_work
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: mat
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(group_dist_d1_type), INTENT(IN)               :: virt_dist

      CHARACTER(LEN=*), PARAMETER :: routineN = 'redistribute_into_subgroups'

      INTEGER :: aux_counter, aux_global, aux_local, aux_proc, avirt, dimen_RI, handle, handle2, &
                 ia_global, ia_local, iocc, max_number_recv, max_number_send, my_aux_end, my_aux_size, &
                 my_aux_start, my_process_column, my_process_row, my_virt_end, my_virt_size, &
                 my_virt_start, proc, proc_shift, recv_proc, send_proc, virt_counter, virt_proc, group_size
      INTEGER, ALLOCATABLE, DIMENSION(:) :: data2send, recv_col_indices, &
                                            recv_row_indices, send_aux_indices, send_virt_indices, virt2send
      INTEGER, DIMENSION(2)                              :: recv_shape
      INTEGER, DIMENSION(:), POINTER                     :: aux_distribution_fm, col_indices, &
                                                            ia_distribution_fm, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: mpi2blacs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET   :: recv_buffer, send_buffer
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: recv_ptr, send_ptr
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm_mat, &
                          nrow_locals=aux_distribution_fm, &
                          col_indices=col_indices, &
                          row_indices=row_indices, &
                          ncol_locals=ia_distribution_fm, &
                          context=context, &
                          nrow_global=dimen_RI, &
                          para_env=para_env)

      IF (exchange_work%homo(ispin) <= 0 .OR. exchange_work%virtual(ispin) <= 0) THEN
         CALL get_group_dist(virt_dist, exchange_work%para_env_sub%mepos, my_virt_start, my_virt_end, my_virt_size)
         ALLOCATE (mat(exchange_work%homo(ispin), my_virt_size, dimen_RI))
         CALL timestop(handle)
         RETURN
      END IF

      group_size = exchange_work%para_env_sub%num_pe

      CALL timeset(routineN//"_prep", handle2)
      CALL get_group_dist(exchange_work%aux_func_dist, exchange_work%color_sub, my_aux_start, my_aux_end, my_aux_size)
      CALL get_group_dist(virt_dist, exchange_work%para_env_sub%mepos, my_virt_start, my_virt_end, my_virt_size)
      CALL context%get(my_process_column=my_process_column, my_process_row=my_process_row, mpi2blacs=mpi2blacs)

      ! Determine the number of columns to send
      ALLOCATE (send_aux_indices(MAXVAL(exchange_work%aux2send)))
      ALLOCATE (virt2send(0:group_size - 1))
      virt2send = 0
      DO ia_local = 1, ia_distribution_fm(my_process_column)
         ia_global = col_indices(ia_local)
         avirt = MOD(ia_global - 1, exchange_work%virtual(ispin)) + 1
         proc = group_dist_proc(virt_dist, avirt)
         virt2send(proc) = virt2send(proc) + 1
      END DO

      ALLOCATE (data2send(0:para_env%num_pe - 1))
      DO aux_proc = 0, exchange_work%ngroup - 1
      DO virt_proc = 0, group_size - 1
         data2send(aux_proc*group_size + virt_proc) = exchange_work%aux2send(aux_proc)*virt2send(virt_proc)
      END DO
      END DO

      ALLOCATE (send_virt_indices(MAXVAL(virt2send)))
      max_number_send = MAXVAL(data2send)

      ALLOCATE (send_buffer(INT(max_number_send, KIND=int_8)*exchange_work%homo(ispin)))
      max_number_recv = max_number_send
      CALL para_env%max(max_number_recv)
      ALLOCATE (recv_buffer(max_number_recv))

      ALLOCATE (mat(my_virt_size, exchange_work%homo(ispin), my_aux_size))

      CALL timestop(handle2)

      CALL timeset(routineN//"_own", handle2)
      ! Start with own data
      DO aux_local = 1, aux_distribution_fm(my_process_row)
         aux_global = row_indices(aux_local)
         IF (aux_global < my_aux_start .OR. aux_global > my_aux_end) CYCLE
         DO ia_local = 1, ia_distribution_fm(my_process_column)
            ia_global = fm_mat%matrix_struct%col_indices(ia_local)

            iocc = (ia_global - 1)/exchange_work%virtual(ispin) + 1
            avirt = MOD(ia_global - 1, exchange_work%virtual(ispin)) + 1

            IF (my_virt_start > avirt .OR. my_virt_end < avirt) CYCLE

            mat(avirt - my_virt_start + 1, iocc, aux_global - my_aux_start + 1) = fm_mat%local_data(aux_local, ia_local)
         END DO
      END DO
      CALL timestop(handle2)

      DO proc_shift = 1, para_env%num_pe - 1
         send_proc = MODULO(para_env%mepos + proc_shift, para_env%num_pe)
         recv_proc = MODULO(para_env%mepos - proc_shift, para_env%num_pe)

         CALL timeset(routineN//"_pack_buffer", handle2)
         send_ptr(1:virt2send(MOD(send_proc, group_size)), &
                  1:exchange_work%aux2send(send_proc/group_size)) => &
            send_buffer(1:INT(virt2send(MOD(send_proc, group_size)), KIND=int_8)* &
                        exchange_work%aux2send(send_proc/group_size))
! Pack send buffer
         aux_counter = 0
         DO aux_local = 1, aux_distribution_fm(my_process_row)
            aux_global = row_indices(aux_local)
            proc = group_dist_proc(exchange_work%aux_func_dist, aux_global)
            IF (proc /= send_proc/group_size) CYCLE
            aux_counter = aux_counter + 1
            virt_counter = 0
            DO ia_local = 1, ia_distribution_fm(my_process_column)
               ia_global = col_indices(ia_local)
               avirt = MOD(ia_global - 1, exchange_work%virtual(ispin)) + 1

               proc = group_dist_proc(virt_dist, avirt)
               IF (proc /= MOD(send_proc, group_size)) CYCLE
               virt_counter = virt_counter + 1
               send_ptr(virt_counter, aux_counter) = fm_mat%local_data(aux_local, ia_local)
               send_virt_indices(virt_counter) = ia_global
            END DO
            send_aux_indices(aux_counter) = aux_global
         END DO
         CALL timestop(handle2)

         CALL timeset(routineN//"_ex_size", handle2)
         recv_shape = [1, 1]
         CALL para_env%sendrecv(SHAPE(send_ptr), send_proc, recv_shape, recv_proc)
         CALL timestop(handle2)

         IF (SIZE(send_ptr) == 0) send_proc = mp_proc_null
         IF (PRODUCT(recv_shape) == 0) recv_proc = mp_proc_null

         CALL timeset(routineN//"_ex_idx", handle2)
         ALLOCATE (recv_row_indices(recv_shape(1)), recv_col_indices(recv_shape(2)))
         CALL para_env%sendrecv(send_virt_indices(1:virt_counter), send_proc, recv_row_indices, recv_proc)
         CALL para_env%sendrecv(send_aux_indices(1:aux_counter), send_proc, recv_col_indices, recv_proc)
         CALL timestop(handle2)

         ! Prepare pointer to recv buffer (consider transposition while packing the send buffer)
         recv_ptr(1:recv_shape(1), 1:MAX(1, recv_shape(2))) => recv_buffer(1:recv_shape(1)*MAX(1, recv_shape(2)))

         CALL timeset(routineN//"_sendrecv", handle2)
! Perform communication
         CALL para_env%sendrecv(send_ptr, send_proc, recv_ptr, recv_proc)
         CALL timestop(handle2)

         IF (recv_proc == mp_proc_null) THEN
            DEALLOCATE (recv_row_indices, recv_col_indices)
            CYCLE
         END IF

         CALL timeset(routineN//"_unpack", handle2)
! Unpack receive buffer
         DO aux_local = 1, SIZE(recv_col_indices)
            aux_global = recv_col_indices(aux_local)

            DO ia_local = 1, SIZE(recv_row_indices)
               ia_global = recv_row_indices(ia_local)

               iocc = (ia_global - 1)/exchange_work%virtual(ispin) + 1
               avirt = MOD(ia_global - 1, exchange_work%virtual(ispin)) + 1

               mat(avirt - my_virt_start + 1, iocc, aux_global - my_aux_start + 1) = recv_ptr(ia_local, aux_local)
            END DO
         END DO
         CALL timestop(handle2)

         IF (ALLOCATED(recv_row_indices)) DEALLOCATE (recv_row_indices)
         IF (ALLOCATED(recv_col_indices)) DEALLOCATE (recv_col_indices)
      END DO

      DEALLOCATE (send_aux_indices, send_virt_indices)

      CALL timestop(handle)

   END SUBROUTINE redistribute_into_subgroups

END MODULE rpa_exchange
