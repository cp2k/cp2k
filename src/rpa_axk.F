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
        dbcsr_copy, dbcsr_create, dbcsr_init_p, dbcsr_multiply, dbcsr_p_type, dbcsr_release, &
        dbcsr_set, dbcsr_trace, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat_general,&
                                              cp_fm_type
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_ri,                          ONLY: hfx_ri_update_ks
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

   PUBLIC :: compute_axk_ener

   TYPE rpa_axk_exchange_type
      TYPE(qs_environment_type), POINTER :: qs_env => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_work_ao => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_2d => NULL()
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data => NULL()
      TYPE(mp_para_env_type), POINTER                    :: para_env_sub => NULL()
      TYPE(section_vals_type), POINTER                   :: hfx_sections => NULL()
      INTEGER :: n_rep_hf = 0
      INTEGER :: my_group_L_size = 0
      LOGICAL :: my_recalc_hfx_integrals = .FALSE.
      REAL(KIND=dp) :: eps_filter = 0.0_dp
   CONTAINS
      PROCEDURE, PASS(axk_exchange), NON_OVERRIDABLE :: create => hfx_create_subgroup
      PROCEDURE, PASS(axk_exchange), NON_OVERRIDABLE :: integrate => integrate_exchange
      PROCEDURE, PASS(axk_exchange), NON_OVERRIDABLE :: release => hfx_release_subgroup
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief Main driver for RPA-AXK energies
!> \param qs_env ...
!> \param fm_mat_Q ...
!> \param fm_mat_Q_gemm ...
!> \param dimen_RI ...
!> \param dimen_ia ...
!> \param para_env_sub ...
!> \param para_env_RPA ...
!> \param eig ...
!> \param fm_mat_S ...
!> \param homo ...
!> \param virtual ...
!> \param omega ...
!> \param mp2_env ...
!> \param mat_munu ...
!> \param unit_nr ...
!> \param e_axk_corr ...          AXK energy correctrion for a quadrature point
!> \author Vladimir Rybkin, 07/2016
! **************************************************************************************************
   SUBROUTINE compute_axk_ener(qs_env, fm_mat_Q, fm_mat_Q_gemm, dimen_RI, dimen_ia, &
                               para_env_sub, para_env_RPA, &
                               eig, fm_mat_S, homo, virtual, omega, &
                               mp2_env, mat_munu, unit_nr, e_axk_corr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_Q, fm_mat_Q_gemm
      INTEGER, INTENT(IN)                                :: dimen_RI
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dimen_ia
      TYPE(mp_para_env_type), POINTER                    :: para_env_sub, para_env_RPA
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eig
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_S
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      REAL(KIND=dp), INTENT(IN)                          :: omega
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_munu
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), INTENT(INOUT)                       :: e_axk_corr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_axk_ener'
      REAL(KIND=dp), PARAMETER                           :: thresh = 0.0000001_dp

      INTEGER :: color_sub, handle, iib, iitmp(2), ispin, kkb, L_counter, my_group_L_end, &
         my_group_L_size, my_group_L_start, ngroup, nspins
      REAL(KIND=dp)                                      :: eps_filter, trace_corr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_mat_Q_tmp, fm_mat_R_half, &
                                                            fm_mat_R_half_gemm, fm_mat_U
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_mat_Gamma_3
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dbcsr_Gamma_3, dbcsr_Gamma_munu_P
      TYPE(dbcsr_type)                                   :: dbcsr_Gamma_inu_P
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_o, mo_coeff_v
      TYPE(rpa_axk_exchange_type)                        :: axk_exchange

      CALL timeset(routineN, handle)

      nspins = SIZE(dimen_ia)

      ! Eigenvalues
      ALLOCATE (eigenval(dimen_RI))
      eigenval = 0.0_dp
      ! create the R_half and U matrices with a different blacs env similar to Q
      ! and a tmp_Q needed for diagonalization

      NULLIFY (fm_struct)

      CALL cp_fm_get_info(matrix=fm_mat_Q, &
                          matrix_struct=fm_struct)
      CALL cp_fm_create(fm_mat_U, fm_struct, name="fm_mat_U")
      CALL cp_fm_create(fm_mat_R_half, fm_struct, name="fm_mat_R_half")
      CALL cp_fm_create(fm_mat_Q_tmp, fm_struct, name="fm_mat_Q_tmp")
      CALL cp_fm_set_all(matrix=fm_mat_Q_tmp, alpha=0.0_dp)
      CALL cp_fm_set_all(matrix=fm_mat_U, alpha=0.0_dp)
      CALL cp_fm_set_all(matrix=fm_mat_R_half, alpha=0.0_dp)

      ! Copy Q to Q_tmp
      CALL cp_fm_to_fm(fm_mat_Q, fm_mat_Q_tmp)

      CALL cp_fm_scale(0.50_dp, fm_mat_Q_tmp)
      ! Diagonalize Q
      CALL choose_eigv_solver(fm_mat_Q_tmp, fm_mat_U, eigenval)

      !Calculate diagonal matrix for R_half

      ! U*diag stored in U, whereas eigenvectors are in fm_mat_Q_tmp
      CALL cp_fm_to_fm(fm_mat_U, fm_mat_Q_tmp)

      ! Manipulate eigenvalues to get diagonal matrix
      DO iib = 1, dimen_RI
         IF (ABS(eigenval(iib)) .GE. thresh) THEN
            eigenval(iib) = &
               SQRT((1.0_dp/(eigenval(iib)**2))*LOG(1.0_dp + eigenval(iib)) &
                    - 1.0_dp/(eigenval(iib)*(eigenval(iib) + 1.0_dp)))
         ELSE
            eigenval(iib) = sqrthalf
         END IF
      END DO

      CALL cp_fm_column_scale(fm_mat_U, eigenval)

      ! Release memory
      DEALLOCATE (eigenval)

      ! Get R_half by multiplication
      CALL parallel_gemm(transa="N", transb="T", m=dimen_RI, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                         matrix_a=fm_mat_U, matrix_b=fm_mat_Q_tmp, beta=0.0_dp, &
                         matrix_c=fm_mat_R_half)

      ! Release full matrix stuff
      CALL cp_fm_release(fm_mat_Q_tmp)
      CALL cp_fm_release(fm_mat_U)

      ! Redistribute fm_mat_R_half for "rectangular" multiplication: ia*P P*P
      CALL cp_fm_create(fm_mat_R_half_gemm, fm_mat_Q_gemm%matrix_struct)
      CALL cp_fm_set_all(matrix=fm_mat_R_half_gemm, alpha=0.0_dp)

      CALL cp_fm_to_fm_submat_general(fm_mat_R_half, fm_mat_R_half_gemm, dimen_RI, dimen_RI, 1, 1, 1, 1, &
                                      fm_mat_R_half%matrix_struct%context)

      CALL cp_fm_release(fm_mat_R_half)

      ALLOCATE (fm_mat_Gamma_3(nspins))
      DO ispin = 1, nspins
         ! get info of fm_mat_S and initialize Gamma_3
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, template_fmstruct=fm_mat_S(ispin)%matrix_struct, &
                                  nrow_global=dimen_ia(ispin), ncol_global=dimen_RI)
         CALL cp_fm_create(fm_mat_Gamma_3(ispin), fm_struct)
         CALL cp_fm_struct_release(fm_struct)
         CALL cp_fm_set_all(matrix=fm_mat_Gamma_3(ispin), alpha=0.0_dp)

         ! Update G with a new value of Omega: in practice, it is G*S

         ! Scale fm_work_iaP
         CALL calc_fm_mat_S_rpa(fm_mat_S(ispin), .TRUE., virtual(ispin), eig(:, ispin), &
                                homo(ispin), omega, 0.0_dp)

         ! Calculate Gamma_3: Gamma_3 = G*S*R^(1/2) = G*S*R^(1/2) )
         CALL parallel_gemm(transa="T", transb="N", m=dimen_ia(ispin), n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                            matrix_a=fm_mat_S(ispin), matrix_b=fm_mat_R_half_gemm, beta=0.0_dp, &
                            matrix_c=fm_mat_Gamma_3(ispin))

         ! Remove extra factor from S after the multiplication (to return to the original matrix)
         CALL remove_scaling_factor_rpa(fm_mat_S(ispin), virtual(ispin), eig(:, ispin), homo(ispin), omega)
      END DO

      ! Release full matrix stuff
      CALL cp_fm_release(fm_mat_R_half_gemm)

      ! Get aux sizes
      ngroup = para_env_RPA%num_pe/para_env_sub%num_pe
      color_sub = para_env_RPA%mepos/para_env_sub%num_pe

      iitmp = get_limit(dimen_RI, ngroup, color_sub)
      my_group_L_start = iitmp(1)
      my_group_L_end = iitmp(2)
      my_group_L_size = iitmp(2) - iitmp(1) + 1

      DO ispin = 1, nspins
         ! Retrieve mo coefficients in dbcsr format
         NULLIFY (mo_coeff_o, mo_coeff_v)
         mo_coeff_o => mp2_env%ri_rpa%mo_coeff_o
         mo_coeff_v => mp2_env%ri_rpa%mo_coeff_v

         ! Copy Gamma_ia_P^3 to dbcsr matrix set
         CALL gamma_fm_to_dbcsr(fm_mat_Gamma_3(ispin), dbcsr_Gamma_3, para_env_RPA, para_env_sub, &
                                homo(ispin), virtual(ispin), mo_coeff_o, ngroup, my_group_L_start, &
                                my_group_L_end, my_group_L_size)

         ! Create more dbcsr matrices

         NULLIFY (dbcsr_Gamma_munu_P)
         CALL dbcsr_allocate_matrix_set(dbcsr_Gamma_munu_P, my_group_L_size)
         eps_filter = mp2_env%mp2_gpw%eps_filter

         CALL dbcsr_create(dbcsr_Gamma_inu_P, template=mo_coeff_o)
         CALL dbcsr_copy(dbcsr_Gamma_inu_P, mo_coeff_o)
         CALL dbcsr_set(dbcsr_Gamma_inu_P, 0.0_dp)

         L_counter = 0
         DO kkb = my_group_L_start, my_group_L_end
            L_counter = L_counter + 1
            ! Init Gamma_3 in AO basis
            ALLOCATE (dbcsr_Gamma_munu_P(L_counter)%matrix)
            CALL dbcsr_init_p(dbcsr_Gamma_munu_P(L_counter)%matrix)
            CALL dbcsr_create(dbcsr_Gamma_munu_P(L_counter)%matrix, template=mat_munu%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_copy(dbcsr_Gamma_munu_P(L_counter)%matrix, mat_munu%matrix)
            CALL dbcsr_set(dbcsr_Gamma_munu_P(L_counter)%matrix, 0.0_dp)
         END DO

         CALL axk_exchange%create(qs_env, mat_munu, para_env_sub, my_group_L_size, eps_filter)

      !! Loup over auxiliary basis functions: multiplication
         L_counter = 0
         DO kkb = my_group_L_start, my_group_L_end
            L_counter = L_counter + 1
            ! Do dbcsr multiplication: transform the virtual index
            CALL dbcsr_multiply("N", "T", 1.0_dp, mo_coeff_v, dbcsr_Gamma_3(L_counter)%matrix, &
                                0.0_dp, dbcsr_Gamma_inu_P, filter_eps=eps_filter)

            CALL dbcsr_release(dbcsr_Gamma_3(L_counter)%matrix)
            DEALLOCATE (dbcsr_Gamma_3(L_counter)%matrix)

            !Do dbcsr multiplication: transform the occupied index
            CALL dbcsr_multiply("N", "T", 1.0_dp, dbcsr_Gamma_inu_P, mo_coeff_o, &
                                0.0_dp, dbcsr_Gamma_munu_P(L_counter)%matrix, filter_eps=eps_filter)

            CALL dbcsr_trace(dbcsr_Gamma_munu_P(L_counter)%matrix, trace_corr)
         END DO
         DEALLOCATE (dbcsr_Gamma_3)
         CALL dbcsr_release(dbcsr_Gamma_inu_P)

         ! Contract DM with exchange integrals
         CALL axk_exchange%integrate(dbcsr_Gamma_munu_P, e_axk_corr)

         L_counter = 0
         DO kkb = my_group_L_start, my_group_L_end
            L_counter = L_counter + 1
            CALL dbcsr_release(dbcsr_Gamma_munu_P(L_counter)%matrix)
            DEALLOCATE (dbcsr_Gamma_munu_P(L_counter)%matrix)
         END DO
         DEALLOCATE (dbcsr_Gamma_munu_P)
      END DO

      CALL axk_exchange%release()

      ! Print AXK correlation energy to the file
      IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T68,F25.14,A4)') 'AXK correlation energy for a quadrature point:', &
         e_axk_corr, ' a.u.'

      CALL timestop(handle)

   END SUBROUTINE compute_axk_ener

! **************************************************************************************************
!> \brief Contract RPA-AXK density matrix with HF exchange integrals and evaluate the correction
!> \param axk_exchange ...
!> \param dbcsr_Gamma_munu_P ...  AXK density matrix in AO basis to be contracted
!> \param axk_corr ...  The AXK energy correction
!> \author Vladimir Rybkin, 08/2016
! **************************************************************************************************
   SUBROUTINE integrate_exchange(axk_exchange, dbcsr_Gamma_munu_P, axk_corr)
      CLASS(rpa_axk_exchange_type), INTENT(INOUT) :: axk_exchange
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dbcsr_Gamma_munu_P
      REAL(KIND=dp), INTENT(OUT)                         :: axk_corr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_exchange'

      INTEGER                                            :: aux, handle, irep, ns
      REAL(KIND=dp)                                      :: e_axk_P, ehfx
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_2d

      CALL timeset(routineN, handle)

      ns = SIZE(axk_exchange%rho_work_ao)

      ! The loop over auxiliary basis functions
      axk_corr = 0.0_dp
      DO aux = 1, axk_exchange%my_group_L_size

         CALL dbcsr_copy(axk_exchange%rho_work_ao(1)%matrix, dbcsr_Gamma_munu_P(aux)%matrix)

         DO irep = 1, axk_exchange%n_rep_hf
            rho_ao_2d(1:ns, 1:1) => axk_exchange%rho_work_ao(1:ns)

            CALL dbcsr_set(axk_exchange%mat_2d(1, 1)%matrix, 0.0_dp)

            IF (axk_exchange%x_data(irep, 1)%do_hfx_ri) THEN
               CALL hfx_ri_update_ks(axk_exchange%qs_env, axk_exchange%x_data(irep, 1)%ri_data, axk_exchange%mat_2d, ehfx, &
                                     rho_ao=rho_ao_2d, geometry_did_change=axk_exchange%my_recalc_hfx_integrals, &
                                     nspins=ns, hf_fraction=axk_exchange%x_data(irep, 1)%general_parameter%fraction)

            ELSE
               CALL integrate_four_center(axk_exchange%qs_env, axk_exchange%x_data, axk_exchange%mat_2d, ehfx, &
                                          rho_ao_2d, axk_exchange%hfx_sections, &
                                          axk_exchange%para_env_sub, axk_exchange%my_recalc_hfx_integrals, irep, .TRUE., &
                                          ispin=1)
            END IF
         END DO

         axk_exchange%my_recalc_hfx_integrals = .FALSE.
         ! One more dbcsr multiplication and trace
         CALL dbcsr_multiply("T", "N", 1.0_dp, axk_exchange%mat_2d(1, 1)%matrix, axk_exchange%rho_work_ao(1)%matrix, &
                             0.0_dp, dbcsr_Gamma_munu_P(aux)%matrix, filter_eps=axk_exchange%eps_filter)
         CALL dbcsr_trace(dbcsr_Gamma_munu_P(aux)%matrix, e_axk_p)
         axk_corr = axk_corr + e_axk_P
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
      axk_exchange%para_env_sub => para_env_sub
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

      ! create a working rho environment
      CALL dbcsr_allocate_matrix_set(axk_exchange%rho_work_ao, dft_control%nspins)
      DO ispin = 1, dft_control%nspins
         ALLOCATE (axk_exchange%rho_work_ao(ispin)%matrix)
         CALL dbcsr_init_p(axk_exchange%rho_work_ao(ispin)%matrix)
         CALL dbcsr_create(axk_exchange%rho_work_ao(ispin)%matrix, template=mat_munu%matrix)
      END DO

      CALL dbcsr_allocate_matrix_set(axk_exchange%mat_2d, dft_control%nspins, 1)
      DO ispin = 1, dft_control%nspins
         ALLOCATE (axk_exchange%mat_2d(ispin, 1)%matrix)
         CALL dbcsr_init_p(axk_exchange%mat_2d(ispin, 1)%matrix)
         CALL dbcsr_create(axk_exchange%mat_2d(ispin, 1)%matrix, template=mat_munu%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(axk_exchange%mat_2d(ispin, 1)%matrix, mat_munu%matrix)
      END DO

      CALL timestop(handle)

   END SUBROUTINE hfx_create_subgroup

! **************************************************************************************************
!> \brief ...
!> \param axk_exchange ...
! **************************************************************************************************
   SUBROUTINE hfx_release_subgroup(axk_exchange)
      CLASS(rpa_axk_exchange_type), INTENT(INOUT) :: axk_exchange

      INTEGER :: ispin

      NULLIFY (axk_exchange%qs_env, axk_exchange%para_env_sub, axk_exchange%hfx_sections)

      CALL hfx_release(axk_exchange%x_data)
      NULLIFY (axk_exchange%x_data)

      DO ispin = 1, SIZE(axk_exchange%mat_2d, 1)
         CALL dbcsr_release(axk_exchange%mat_2d(ispin, 1)%matrix)
         CALL dbcsr_release(axk_exchange%mat_2d(ispin, 1)%matrix)
         DEALLOCATE (axk_exchange%mat_2d(ispin, 1)%matrix)
      END DO
      DEALLOCATE (axk_exchange%mat_2d)
      NULLIFY (axk_exchange%mat_2d)

      DO ispin = 1, SIZE(axk_exchange%rho_work_ao, 1)
         CALL dbcsr_release(axk_exchange%rho_work_ao(ispin)%matrix)
         DEALLOCATE (axk_exchange%rho_work_ao(ispin)%matrix)
      END DO
      DEALLOCATE (axk_exchange%rho_work_ao)
      NULLIFY (axk_exchange%rho_work_ao)
   END SUBROUTINE hfx_release_subgroup

END MODULE rpa_axk
