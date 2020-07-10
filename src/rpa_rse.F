!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to compute singles correction to RPA (RSE)
!> \par History
!>      08.2019 created [Vladimir Rybkin]
!> \author Vladimir Rybkin
! **************************************************************************************************
MODULE rpa_rse

   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set,&
                                              dbcsr_type_symmetric
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
   USE mp2_types,                       ONLY: mp2_type
   USE pw_types,                        ONLY: pw_p_type,&
                                              pw_release
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_ks_utils,                     ONLY: compute_matrix_vxc
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_rse'

   PUBLIC :: rse_energy

CONTAINS

! **************************************************************************************************
!> \brief Single excitations energy corrections for RPA
!> \param qs_env ...
!> \param mp2_env ...
!> \param para_env ...
!> \param dft_control ...
!> \param mo_coeff ...
!> \param nmo ...
!> \param homo ...
!> \param Eigenval ...
!> \author Vladimir Rybkin, 08/2019
! **************************************************************************************************
   SUBROUTINE rse_energy(qs_env, mp2_env, para_env, dft_control, &
                         mo_coeff, nmo, homo, Eigenval)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(mp2_type), INTENT(IN), POINTER                :: mp2_env
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(dft_control_type), INTENT(IN), POINTER        :: dft_control
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: mo_coeff
      INTEGER, INTENT(IN)                                :: nmo
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rse_energy'

      INTEGER                                            :: dimen, handle, i_global, iiB, ispin, &
                                                            j_global, jjB, n_rep_hf, ncol_local, &
                                                            nrow_local, nspins
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: beta, do_hfx, hfx_treat_lsd_in_core
      REAL(KIND=dp)                                      :: coeff, corr, rse_corr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diag_diff
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: fm_P_mu_nu, fm_X_mo, fm_XC_mo
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: fm_ao, fm_ao_mo
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_mu_nu, matrix_s, rho_ao
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input

      CALL timeset(routineN, handle)

      nspins = dft_control%nspins
      beta = (nspins == 2)

      ! Pick the diagonal terms
      CALL cp_fm_get_info(matrix=mo_coeff(1)%matrix, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ! start collecting stuff
      dimen = nmo
      NULLIFY (input, matrix_s, blacs_env, rho, energy, sab_orb)
      CALL get_qs_env(qs_env, &
                      input=input, &
                      matrix_s=matrix_s, &
                      blacs_env=blacs_env, &
                      rho=rho, &
                      energy=energy, &
                      sab_orb=sab_orb)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      ! hfx section
      NULLIFY (hfx_sections)
      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")
      CALL section_vals_get(hfx_sections, explicit=do_hfx, n_repetition=n_rep_hf)
      IF (do_hfx) THEN
         CALL section_vals_val_get(hfx_sections, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core, &
                                   i_rep_section=1)
      END IF

      ! create work array
      NULLIFY (mat_mu_nu)
      CALL dbcsr_allocate_matrix_set(mat_mu_nu, nspins)
      DO ispin = 1, nspins
         ALLOCATE (mat_mu_nu(ispin)%matrix)
         CALL dbcsr_create(matrix=mat_mu_nu(ispin)%matrix, template=matrix_s(1)%matrix, name="T_mu_nu", &
                           matrix_type=dbcsr_type_symmetric, nze=0)
         CALL cp_dbcsr_alloc_block_from_nbl(mat_mu_nu(ispin)%matrix, sab_orb)
         CALL dbcsr_set(mat_mu_nu(ispin)%matrix, 0.0_dp)
      END DO

      ! Dense (full) matrices
      ALLOCATE (fm_P_mu_nu(nspins))
      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                               nrow_global=dimen, ncol_global=dimen)
      DO ispin = 1, nspins
         NULLIFY (fm_P_mu_nu(ispin)%matrix)
         CALL cp_fm_create(fm_P_mu_nu(ispin)%matrix, fm_struct_tmp, name="P_mu_nu")
         CALL cp_fm_set_all(fm_P_mu_nu(ispin)%matrix, 0.0_dp)
      END DO
      CALL cp_fm_struct_release(fm_struct_tmp)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                               nrow_global=dimen, ncol_global=dimen)
      ALLOCATE (fm_X_mo(nspins), fm_XC_mo(nspins))
      DO ispin = 1, nspins
         NULLIFY (fm_X_mo(ispin)%matrix, fm_XC_mo(ispin)%matrix)
         CALL cp_fm_create(fm_X_mo(ispin)%matrix, fm_struct_tmp, name="f_X_mo")
         CALL cp_fm_create(fm_XC_mo(ispin)%matrix, fm_struct_tmp, name="f_XC_mo")
         CALL cp_fm_set_all(fm_X_mo(ispin)%matrix, 0.0_dp)
         CALL cp_fm_set_all(fm_XC_mo(ispin)%matrix, 0.0_dp)
      END DO
      CALL cp_fm_struct_release(fm_struct_tmp)

      NULLIFY (fm_ao)
      CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                               nrow_global=dimen, ncol_global=dimen)
      CALL cp_fm_create(fm_ao, fm_struct_tmp, name="f_ao")
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_fm_set_all(fm_ao, 0.0_dp)
      NULLIFY (fm_ao_mo)
      CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                               nrow_global=dimen, ncol_global=dimen)
      CALL cp_fm_create(fm_ao_mo, fm_struct_tmp, name="f_ao_mo")
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_fm_set_all(fm_ao_mo, 0.0_dp)

      !
      !     Ready with preparations, do the real staff
      !

      ! Obtain density matrix like quantity

      coeff = 1.0_dp
      IF (beta) coeff = 0.5_dp
      DO ispin = 1, nspins
         CALL cp_gemm(transa='N', transb='T', m=dimen, n=dimen, k=dimen, alpha=coeff, &
                      matrix_a=mo_coeff(ispin)%matrix, matrix_b=mo_coeff(ispin)%matrix, &
                      beta=0.0_dp, matrix_c=fm_P_mu_nu(ispin)%matrix)
      END DO

      ! Calculate exact exchange contribution
      CALL exchange_contribution(qs_env, para_env, dimen, mo_coeff, &
                                 hfx_sections, n_rep_hf, &
                                 rho, mat_mu_nu, fm_P_mu_nu, &
                                 fm_ao, fm_X_mo, fm_ao_mo, &
                                 mp2_env%ri_mp2%free_hfx_buffer)

      ! Calculate DFT exchange-correlation contribution
      CALL xc_contribution(qs_env, fm_ao, fm_ao_mo, fm_XC_mo, mo_coeff, dimen)

      ALLOCATE (diag_diff(dimen))
      rse_corr = 0.0_dp

      DO ispin = 1, nspins
         ! Compute the correction matrix: it is stored in fm_X_mo
         CALL cp_fm_scale_and_add(1.0_dp, fm_X_mo(ispin)%matrix, -1.0_dp, fm_XC_mo(ispin)%matrix)

         ! Pick the diagonal terms
         CALL cp_fm_get_diag(fm_X_mo(ispin)%matrix, diag_diff)

         ! Compute the correction
         CALL cp_fm_get_info(matrix=fm_X_mo(ispin)%matrix, &
                             nrow_local=nrow_local, &
                             ncol_local=ncol_local, &
                             row_indices=row_indices, &
                             col_indices=col_indices)

         corr = 0.0_dp

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,i_global,j_global) &
!$OMP             REDUCTION(+: corr) &
!$OMP             SHARED(ncol_local,nrow_local,col_indices,row_indices,diag_diff,eigenval,fm_X_mo,homo,ispin)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)
               IF ((i_global .LE. homo(ispin)) .AND. (j_global .GT. homo(ispin))) THEN
                  corr = corr + fm_X_mo(ispin)%matrix%local_data(iib, jjb)**2.0_dp/ &
                         (eigenval(i_global, ispin) - eigenval(j_global, ispin) - diag_diff(i_global) + diag_diff(j_global))
               ENDIF
            END DO
         END DO
!$OMP END PARALLEL DO

         rse_corr = rse_corr + corr
      END DO

      CALL mp_sum(rse_corr, para_env%group)

      rse_corr = rse_corr/nspins

      mp2_env%ri_rpa%rse_corr_diag = rse_corr

      CALL non_diag_rse(fm_X_mo, eigenval, dimen, homo, para_env, blacs_env, rse_corr)

      mp2_env%ri_rpa%rse_corr = rse_corr

      ! Release staff
      DEALLOCATE (diag_diff)
      CALL cp_fm_release(fm_ao)
      CALL cp_fm_release(fm_ao_mo)
      DO ispin = 1, nspins
         CALL cp_fm_release(fm_P_mu_nu(ispin)%matrix)
         CALL cp_fm_release(fm_X_mo(ispin)%matrix)
         CALL cp_fm_release(fm_XC_mo(ispin)%matrix)
         CALL dbcsr_release(mat_mu_nu(ispin)%matrix)
         DEALLOCATE (mat_mu_nu(ispin)%matrix)
      ENDDO
      DEALLOCATE (fm_P_mu_nu, fm_X_mo, fm_XC_mo, mat_mu_nu)

      CALL timestop(handle)

   END SUBROUTINE rse_energy

! **************************************************************************************************
!> \brief HF exchange occupied-virtual matrix
!> \param qs_env ...
!> \param para_env ...
!> \param dimen ...
!> \param mo_coeff ...
!> \param hfx_sections ...
!> \param n_rep_hf ...
!> \param rho_work ...
!> \param mat_mu_nu ...
!> \param fm_P_mu_nu ...
!> \param fm_X_ao ...
!> \param fm_X_mo ...
!> \param fm_X_ao_mo ...
!> \param recalc_hfx_integrals ...
! **************************************************************************************************
   SUBROUTINE exchange_contribution(qs_env, para_env, dimen, mo_coeff, &
                                    hfx_sections, n_rep_hf, &
                                    rho_work, mat_mu_nu, fm_P_mu_nu, &
                                    fm_X_ao, fm_X_mo, fm_X_ao_mo, &
                                    recalc_hfx_integrals)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, INTENT(IN)                                :: dimen
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: mo_coeff
      TYPE(section_vals_type), INTENT(IN), POINTER       :: hfx_sections
      INTEGER, INTENT(IN)                                :: n_rep_hf
      TYPE(qs_rho_type), INTENT(IN), POINTER             :: rho_work
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: mat_mu_nu
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_P_mu_nu
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: fm_X_ao
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_X_mo
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: fm_X_ao_mo
      LOGICAL, INTENT(IN)                                :: recalc_hfx_integrals

      CHARACTER(LEN=*), PARAMETER :: routineN = 'exchange_contribution'

      INTEGER                                            :: handle, irep, is, ns
      LOGICAL                                            :: my_recalc_hfx_integrals
      REAL(KIND=dp)                                      :: ehfx
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: P_mu_nu, rho_work_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_2d, rho_ao_2d

      CALL timeset(routineN, handle)

      my_recalc_hfx_integrals = recalc_hfx_integrals

      CALL qs_rho_get(rho_work, rho_ao=rho_work_ao)
      ns = SIZE(rho_work_ao)
      NULLIFY (P_mu_nu)
      CALL dbcsr_allocate_matrix_set(P_mu_nu, ns)
      DO is = 1, ns
         CALL dbcsr_init_p(P_mu_nu(is)%matrix)
         CALL dbcsr_create(P_mu_nu(is)%matrix, template=rho_work_ao(1)%matrix)
         CALL dbcsr_copy(P_mu_nu(is)%matrix, rho_work_ao(1)%matrix)
         CALL dbcsr_set(P_mu_nu(is)%matrix, 0.0_dp)
      ENDDO

      DO is = 1, ns
         CALL copy_fm_to_dbcsr(fm_P_mu_nu(is)%matrix, P_mu_nu(1)%matrix, keep_sparsity=.TRUE.)

         CALL dbcsr_set(mat_mu_nu(1)%matrix, 0.0_dp)

         DO irep = 1, n_rep_hf
            rho_ao_2d(1:ns, 1:1) => P_mu_nu(1:ns)
            mat_2d(1:ns, 1:1) => mat_mu_nu(1:ns)
            CALL integrate_four_center(qs_env, qs_env%mp2_env%ri_rpa%x_data, mat_2d, ehfx, rho_ao_2d, hfx_sections, &
                                       para_env, my_recalc_hfx_integrals, irep, .TRUE., &
                                       ispin=1)

            my_recalc_hfx_integrals = .FALSE.
         END DO

         ! copy back to fm
         CALL cp_fm_set_all(fm_X_ao, 0.0_dp)
         CALL copy_dbcsr_to_fm(matrix=mat_mu_nu(1)%matrix, fm=fm_X_ao)
         CALL cp_fm_set_all(fm_X_mo(is)%matrix, 0.0_dp)

         ! First index
         CALL cp_gemm('T', 'N', dimen, dimen, dimen, 1.0_dp, &
                      mo_coeff(is)%matrix, fm_X_ao, 0.0_dp, fm_X_ao_mo)

         ! Second index
         CALL cp_gemm('N', 'N', dimen, dimen, dimen, 1.0_dp, &
                      fm_X_ao_mo, mo_coeff(is)%matrix, 1.0_dp, fm_X_mo(is)%matrix)

      END DO

      ! Release dbcsr objects
      DO is = 1, SIZE(P_mu_nu)
         CALL dbcsr_release(P_mu_nu(is)%matrix)
         DEALLOCATE (P_mu_nu(is)%matrix)
      ENDDO
      DEALLOCATE (P_mu_nu)

      CALL timestop(handle)

   END SUBROUTINE exchange_contribution

! **************************************************************************************************
!> \brief Exchange-correlation occupied-virtual matrix
!> \param qs_env ...
!> \param fm_XC_ao ...
!> \param fm_XC_ao_mo ...
!> \param fm_XC_mo ...
!> \param mo_coeff ...
!> \param dimen ...
! **************************************************************************************************
   SUBROUTINE xc_contribution(qs_env, fm_XC_ao, fm_XC_ao_mo, fm_XC_mo, mo_coeff, dimen)
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: fm_XC_ao, fm_XC_ao_mo
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_XC_mo, mo_coeff
      INTEGER, INTENT(IN)                                :: dimen

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'xc_contribution'

      INTEGER                                            :: handle, i
      REAL(KIND=dp)                                      :: exc
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_vxc
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: tau_rspace, v_rspace
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: input, xc_section

      CALL timeset(routineN, handle)

      NULLIFY (matrix_vxc, v_rspace, tau_rspace, input, xc_section, ks_env, &
               rho)
      CALL get_qs_env(qs_env, matrix_vxc=matrix_vxc, input=input, ks_env=ks_env, rho=rho)
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")

      ! Compute XC matrix in AO basis
      CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=xc_section, &
                         vxc_rho=v_rspace, vxc_tau=tau_rspace, exc=exc)

      CALL compute_matrix_vxc(qs_env=qs_env, v_rspace=v_rspace, matrix_vxc=matrix_vxc)

      DO i = 1, SIZE(v_rspace)
         CALL pw_release(v_rspace(i)%pw)
      ENDDO
      DEALLOCATE (v_rspace)

      DO i = 1, SIZE(matrix_vxc)
         CALL cp_fm_set_all(fm_XC_ao, 0.0_dp)
         CALL copy_dbcsr_to_fm(matrix=matrix_vxc(i)%matrix, fm=fm_XC_ao)
         CALL cp_fm_set_all(fm_XC_mo(i)%matrix, 0.0_dp)

         ! First index
         CALL cp_gemm('T', 'N', dimen, dimen, dimen, 1.0_dp, &
                      mo_coeff(i)%matrix, fm_XC_ao, 0.0_dp, fm_XC_ao_mo)

         ! Second index
         CALL cp_gemm('N', 'N', dimen, dimen, dimen, 1.0_dp, &
                      fm_XC_ao_mo, mo_coeff(i)%matrix, 1.0_dp, fm_XC_mo(i)%matrix)

      ENDDO

      DO i = 1, SIZE(matrix_vxc)
         CALL dbcsr_release(matrix_vxc(i)%matrix)
         DEALLOCATE (matrix_vxc(i)%matrix)
      ENDDO
      DEALLOCATE (matrix_vxc)

      CALL timestop(handle)

   END SUBROUTINE xc_contribution

! **************************************************************************************************
!> \brief ...
!> \param fm_F_mo ...
!> \param eigenval ...
!> \param dimen ...
!> \param homo ...
!> \param para_env ...
!> \param blacs_env ...
!> \param rse_corr ...
! **************************************************************************************************
   SUBROUTINE non_diag_rse(fm_F_mo, eigenval, dimen, homo, para_env, &
                           blacs_env, rse_corr)
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_F_mo
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      INTEGER, INTENT(IN)                                :: dimen
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(cp_blacs_env_type), INTENT(IN), POINTER       :: blacs_env
      REAL(KIND=dp), INTENT(OUT)                         :: rse_corr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'non_diag_rse'

      INTEGER                                            :: handle, i_global, iiB, ispin, j_global, &
                                                            jjB, ncol_local, nrow_local, nspins, &
                                                            virtual
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: corr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eig_o, eig_semi_can, eig_v
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: fm_F_oo, fm_F_ov, fm_F_vv, fm_O, fm_tmp, &
                                                            fm_U

      CALL timeset(routineN, handle)

      nspins = SIZE(fm_f_mo)

      DO ispin = 1, nspins
         ! Add eigenvalues on the diagonal
         CALL cp_fm_get_info(matrix=fm_F_mo(ispin)%matrix, &
                             nrow_local=nrow_local, &
                             ncol_local=ncol_local, &
                             row_indices=row_indices, &
                             col_indices=col_indices)

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,i_global,j_global) &
!$OMP             SHARED(ncol_local,nrow_local,col_indices,row_indices,fm_F_mo,eigenval,ispin)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)
               IF (i_global .EQ. j_global) fm_F_mo(ispin)%matrix%local_data(iib, jjb) = &
                  fm_F_mo(ispin)%matrix%local_data(iib, jjb) + eigenval(i_global, ispin)
            END DO
         END DO
!$OMP END PARALLEL DO
      ENDDO

      rse_corr = 0.0_dp

      DO ispin = 1, nspins
         ! Create the occupied-occupied and virtual-virtual blocks, eigenvectors
         NULLIFY (fm_F_oo, fm_O, fm_struct_tmp)
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                                  nrow_global=homo(ispin), ncol_global=homo(ispin))
         CALL cp_fm_create(fm_F_oo, fm_struct_tmp, name="F_oo")
         CALL cp_fm_create(fm_O, fm_struct_tmp, name="O")
         CALL cp_fm_set_all(fm_F_oo, 0.0_dp)
         CALL cp_fm_set_all(fm_O, 0.0_dp)
         CALL cp_fm_struct_release(fm_struct_tmp)

         CALL cp_fm_to_fm_submat(msource=fm_F_mo(ispin)%matrix, mtarget=fm_F_oo, &
                                 nrow=homo(ispin), ncol=homo(ispin), &
                                 s_firstrow=1, s_firstcol=1, &
                                 t_firstrow=1, t_firstcol=1)
         virtual = dimen - homo(ispin)
         NULLIFY (fm_F_vv, fm_U, fm_struct_tmp)
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                                  nrow_global=virtual, ncol_global=virtual)
         CALL cp_fm_create(fm_F_vv, fm_struct_tmp, name="F_vv")
         CALL cp_fm_create(fm_U, fm_struct_tmp, name="U")
         CALL cp_fm_set_all(fm_F_vv, 0.0_dp)
         CALL cp_fm_set_all(fm_U, 0.0_dp)
         CALL cp_fm_struct_release(fm_struct_tmp)

         CALL cp_fm_to_fm_submat(msource=fm_F_mo(ispin)%matrix, mtarget=fm_F_vv, &
                                 nrow=virtual, ncol=virtual, &
                                 s_firstrow=homo(ispin) + 1, s_firstcol=homo(ispin) + 1, &
                                 t_firstrow=1, t_firstcol=1)

         ! Diagonalize occupied-occupied and virtual-virtual matrices
         ALLOCATE (eig_o(homo(ispin)))
         ALLOCATE (eig_v(virtual))
         eig_v = 0.0_dp
         eig_o = 0.0_dp
         CALL choose_eigv_solver(fm_F_oo, fm_O, eig_o)
         CALL choose_eigv_solver(fm_F_vv, fm_U, eig_v)

         ! Collect the eigenvalues to one array
         ALLOCATE (eig_semi_can(dimen))
         eig_semi_can = 0.0_dp
         eig_semi_can(1:homo(ispin)) = eig_o(:)
         eig_semi_can(homo(ispin) + 1:dimen) = eig_v(:)

         ! Create occupied-virtual block
         NULLIFY (fm_F_ov, fm_tmp, fm_struct_tmp)
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                                  nrow_global=homo(ispin), ncol_global=virtual)
         CALL cp_fm_create(fm_F_ov, fm_struct_tmp, name="F_ov")
         CALL cp_fm_create(fm_tmp, fm_struct_tmp, name="tmp")
         CALL cp_fm_set_all(fm_F_ov, 0.0_dp)
         CALL cp_fm_set_all(fm_tmp, 0.0_dp)
         CALL cp_fm_struct_release(fm_struct_tmp)

         CALL cp_fm_to_fm_submat(msource=fm_F_mo(ispin)%matrix, mtarget=fm_F_ov, &
                                 nrow=homo(ispin), ncol=virtual, &
                                 s_firstrow=1, s_firstcol=homo(ispin) + 1, &
                                 t_firstrow=1, t_firstcol=1)

         CALL cp_gemm(transa='N', transb='N', m=homo(ispin), n=virtual, k=homo(ispin), alpha=1.0_dp, &
                      matrix_a=fm_O, matrix_b=fm_F_ov, beta=0.0_dp, matrix_c=fm_tmp)

         CALL cp_gemm(transa='N', transb='N', m=homo(ispin), n=virtual, k=virtual, alpha=1.0_dp, &
                      matrix_a=fm_tmp, matrix_b=fm_U, beta=0.0_dp, matrix_c=fm_F_ov)

         ! Compute the correction
         CALL cp_fm_get_info(matrix=fm_F_ov, &
                             nrow_local=nrow_local, &
                             ncol_local=ncol_local, &
                             row_indices=row_indices, &
                             col_indices=col_indices)
         corr = 0.0_dp
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,i_global,j_global) &
!$OMP             REDUCTION(+:corr) &
!$OMP                SHARED(ncol_local,nrow_local,col_indices,row_indices,fm_F_ov,eig_semi_can,homo,ispin)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)
               corr = corr + fm_F_ov%local_data(iib, jjb)**2.0_dp/ &
                      (eig_semi_can(i_global) - eig_semi_can(j_global + homo(ispin)))
            END DO
         END DO
!$OMP    END PARALLEL DO

         rse_corr = rse_corr + corr

         ! Release
         DEALLOCATE (eig_semi_can)
         DEALLOCATE (eig_o)
         DEALLOCATE (eig_v)

         CALL cp_fm_release(fm_F_ov)
         CALL cp_fm_release(fm_F_oo)
         CALL cp_fm_release(fm_F_vv)
         CALL cp_fm_release(fm_U)
         CALL cp_fm_release(fm_O)
         CALL cp_fm_release(fm_tmp)

      ENDDO

      CALL mp_sum(rse_corr, para_env%group)
      rse_corr = rse_corr/nspins

      CALL timestop(handle)

   END SUBROUTINE non_diag_rse

END MODULE rpa_rse
