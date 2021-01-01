!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to apply a delta pulse for RTP and EMD
! **************************************************************************************************

MODULE rt_delta_pulse
   USE cell_types,                      ONLY: cell_type
   USE commutator_rpnl,                 ONLY: build_com_mom_nl
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                              cp_cfm_gemm
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                              cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert,&
                                              cp_fm_cholesky_reduce,&
                                              cp_fm_cholesky_restore
   USE cp_fm_diag,                      ONLY: cp_fm_syevd
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_filter, dbcsr_init_p, &
        dbcsr_p_type, dbcsr_set, dbcsr_type, dbcsr_type_antisymmetric
   USE input_section_types,             ONLY: section_get_lval,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_moments,                      ONLY: build_berry_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE rt_propagation_types,            ONLY: get_rtp,&
                                              rt_prop_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_delta_pulse'

   PUBLIC :: apply_delta_pulse_periodic, &
             apply_delta_pulse

CONTAINS

! **************************************************************************************************
!> \brief uses perturbation theory to get the proper initial conditions
!> \param qs_env ...
!> \param mos_old ...
!> \param mos_new ...
!> \author Joost & Martin (2011)
! **************************************************************************************************

   SUBROUTINE apply_delta_pulse_periodic(qs_env, mos_old, mos_new)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_old, mos_new

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_delta_pulse_periodic'

      COMPLEX(KIND=dp), DIMENSION(:), POINTER            :: eigenvalues_sqrt
      INTEGER                                            :: handle, icol, idir, irow, ispin, nao, &
                                                            ncol_local, nmo, nrow_global, &
                                                            nrow_local, nvirt
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: com_nl
      REAL(KIND=dp)                                      :: eps_ppnl, factor
      REAL(KIND=dp), DIMENSION(3)                        :: kvec
      REAL(kind=dp), DIMENSION(:), POINTER               :: eigenvalues
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: local_data
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_cfm_type), POINTER                         :: oo_c, oo_v, oo_vt
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: eigenvectors, mat_ks, mat_tmp, momentum, &
                                                            oo_1, oo_2, S_chol, S_inv_fm, tmpS, &
                                                            virtuals
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_rv, matrix_s
      TYPE(dbcsr_type), POINTER                          :: S_inv
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      NULLIFY (cell, mos, rtp, matrix_s, matrix_ks, input, dft_control, particle_set)
      ! we need the overlap and ks matrix for a full diagionalization
      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      mos=mos, &
                      rtp=rtp, &
                      matrix_s=matrix_s, &
                      matrix_ks=matrix_ks, &
                      dft_control=dft_control, &
                      input=input, &
                      particle_set=particle_set)
      com_nl = section_get_lval(section_vals=input, keyword_name="DFT%REAL_TIME_PROPAGATION%COM_NL")
      CALL get_rtp(rtp=rtp, S_inv=S_inv)
      CALL cp_fm_create(S_chol, matrix_struct=rtp%ao_ao_fmstruct, name="S_chol")
      CALL cp_fm_create(S_inv_fm, matrix_struct=rtp%ao_ao_fmstruct, name="S_inv_fm")
      CALL cp_fm_create(tmpS, matrix_struct=rtp%ao_ao_fmstruct)
      CALL copy_dbcsr_to_fm(S_inv, S_inv_fm)
      CALL cp_fm_upper_to_full(S_inv_fm, tmpS)
      CALL cp_fm_get_info(S_inv_fm, nrow_global=nrow_global)
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, S_chol)
      CALL cp_fm_cholesky_decompose(S_chol)
      NULLIFY (mat_ks, eigenvectors, mat_tmp)
      CALL cp_fm_create(mat_ks, matrix_struct=S_inv_fm%matrix_struct, name="mat_ks")
      CALL cp_fm_create(eigenvectors, matrix_struct=S_inv_fm%matrix_struct, name="eigenvectors")

      ! calculate non-local commutator
      IF (com_nl) THEN
         NULLIFY (qs_kind_set, sab_orb, sap_ppnl)
         CALL get_qs_env(qs_env, &
                         sap_ppnl=sap_ppnl, &
                         sab_orb=sab_orb, &
                         qs_kind_set=qs_kind_set)
         eps_ppnl = dft_control%qs_control%eps_ppnl

         NULLIFY (matrix_rv)
         CALL dbcsr_allocate_matrix_set(matrix_rv, 3)
         DO idir = 1, 3
            CALL dbcsr_init_p(matrix_rv(idir)%matrix)
            CALL dbcsr_create(matrix_rv(idir)%matrix, template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(matrix_rv(idir)%matrix, sab_orb)
            CALL dbcsr_set(matrix_rv(idir)%matrix, 0._dp)
         ENDDO
         CALL build_com_mom_nl(qs_kind_set, sab_orb, sap_ppnl, eps_ppnl, particle_set, matrix_rv=matrix_rv)
      ENDIF

      DO ispin = 1, SIZE(matrix_ks)
         ALLOCATE (eigenvalues(nrow_global))
         CALL cp_fm_create(mat_tmp, matrix_struct=S_inv_fm%matrix_struct, name="mat_tmp")

         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, mat_ks)
         CALL cp_fm_cholesky_reduce(mat_ks, S_chol)
         CALL cp_fm_syevd(mat_ks, mat_tmp, eigenvalues)
         CALL cp_fm_cholesky_restore(mat_tmp, nrow_global, S_chol, eigenvectors, "SOLVE")

         ! virtuals
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, nao=nao, nmo=nmo)
         nvirt = nao - nmo
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=S_inv_fm%matrix_struct%para_env, context=S_inv_fm%matrix_struct%context, &
                                  nrow_global=nrow_global, ncol_global=nvirt)
         CALL cp_fm_create(virtuals, matrix_struct=fm_struct_tmp, name="virtuals")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_to_fm(eigenvectors, virtuals, nvirt, nmo + 1, 1)

         ! occupied
         CALL cp_fm_to_fm(eigenvectors, mos_old(2*ispin - 1)%matrix, nmo, 1, 1)

         CALL cp_fm_struct_create(fm_struct_tmp, para_env=S_inv_fm%matrix_struct%para_env, context=S_inv_fm%matrix_struct%context, &
                                  nrow_global=nvirt, ncol_global=nmo)
         CALL cp_fm_create(momentum, matrix_struct=fm_struct_tmp, name="momentum")
         CALL cp_fm_struct_release(fm_struct_tmp)

         ! the momentum operator (in a given direction)
         CALL cp_fm_set_all(mos_new(2*ispin - 1)%matrix, 0.0_dp)

         ! the prefactor (strength of the electric field)
         kvec(:) = cell%h_inv(1, :)*dft_control%rtp_control%delta_pulse_direction(1) + &
                   cell%h_inv(2, :)*dft_control%rtp_control%delta_pulse_direction(2) + &
                   cell%h_inv(3, :)*dft_control%rtp_control%delta_pulse_direction(3)
         kvec = -kvec*twopi*dft_control%rtp_control%delta_pulse_scale

         DO idir = 1, 3
            factor = kvec(idir)
            IF (factor .NE. 0.0_dp) THEN
               CALL cp_dbcsr_sm_fm_multiply(matrix_s(idir + 1)%matrix, mos_old(2*ispin - 1)%matrix, &
                                            mos_old(2*ispin)%matrix, ncol=nmo)
               CALL cp_fm_scale_and_add(1.0_dp, mos_new(2*ispin - 1)%matrix, factor, mos_old(2*ispin)%matrix)
               IF (com_nl) THEN
                  CALL cp_fm_set_all(mos_old(2*ispin)%matrix, 0.0_dp)
                  CALL cp_dbcsr_sm_fm_multiply(matrix_rv(idir)%matrix, mos_old(2*ispin - 1)%matrix, &
                                               mos_old(2*ispin)%matrix, ncol=nmo)
                  CALL cp_fm_scale_and_add(1.0_dp, mos_new(2*ispin - 1)%matrix, factor, mos_old(2*ispin)%matrix)
               ENDIF
            ENDIF
         ENDDO

         CALL cp_gemm('T', 'N', nvirt, nmo, nao, 1.0_dp, virtuals, mos_new(2*ispin - 1)%matrix, 0.0_dp, momentum)

         ! the tricky bit ... rescale by the eigenvalue difference
         CALL cp_fm_get_info(momentum, nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices, local_data=local_data)
         DO icol = 1, ncol_local
            DO irow = 1, nrow_local
               factor = 1/(eigenvalues(col_indices(icol)) - eigenvalues(nmo + row_indices(irow)))
               local_data(irow, icol) = factor*local_data(irow, icol)
            ENDDO
         ENDDO
         CALL cp_fm_release(mat_tmp)
         DEALLOCATE (eigenvalues)

         ! now obtain the initial condition in mos_old
         CALL cp_fm_to_fm(eigenvectors, mos_old(2*ispin - 1)%matrix, nmo, 1, 1)
         CALL cp_gemm("N", "N", nao, nmo, nvirt, 1.0_dp, virtuals, momentum, 0.0_dp, mos_old(2*ispin)%matrix)

         CALL cp_fm_release(virtuals)
         CALL cp_fm_release(momentum)

         ! orthonormalize afterwards
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=S_inv_fm%matrix_struct%para_env, context=S_inv_fm%matrix_struct%context, &
                                  nrow_global=nmo, ncol_global=nmo)
         CALL cp_fm_create(oo_1, matrix_struct=fm_struct_tmp, name="oo_1")
         CALL cp_fm_create(oo_2, matrix_struct=fm_struct_tmp, name="oo_2")
         CALL cp_fm_struct_release(fm_struct_tmp)

         CALL cp_fm_create(mat_tmp, matrix_struct=mos_old(2*ispin - 1)%matrix%matrix_struct, name="tmp_mat")
         ! get the complex overlap matrix
         ! x^T S x + y^T S y + i (-y^TS x+x^T S y)
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, mos_old(2*ispin - 1)%matrix, &
                                      mat_tmp, ncol=nmo)

         CALL cp_gemm("T", "N", nmo, nmo, nao, 1.0_dp, mos_old(2*ispin - 1)%matrix, mat_tmp, 0.0_dp, oo_1)
         CALL cp_gemm("T", "N", nmo, nmo, nao, -1.0_dp, mos_old(2*ispin)%matrix, mat_tmp, 0.0_dp, oo_2)

         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, mos_old(2*ispin)%matrix, &
                                      mat_tmp, ncol=nmo)
         CALL cp_gemm("T", "N", nmo, nmo, nao, 1.0_dp, mos_old(2*ispin)%matrix, mat_tmp, 1.0_dp, oo_1)
         CALL cp_gemm("T", "N", nmo, nmo, nao, 1.0_dp, mos_old(2*ispin - 1)%matrix, mat_tmp, 1.0_dp, oo_2)
         CALL cp_fm_release(mat_tmp)

         CALL cp_cfm_create(oo_c, oo_1%matrix_struct)
         CALL cp_cfm_create(oo_v, oo_1%matrix_struct)
         CALL cp_cfm_create(oo_vt, oo_1%matrix_struct)
         oo_c%local_data = CMPLX(oo_1%local_data, oo_2%local_data, KIND=dp)

         ! compute inv(sqrt(overlap))
         ALLOCATE (eigenvalues(nmo))
         ALLOCATE (eigenvalues_sqrt(nmo))
         CALL cp_cfm_heevd(oo_c, oo_v, eigenvalues)
         eigenvalues_sqrt = CMPLX(1.0_dp/SQRT(eigenvalues), 0.0_dp, dp)
         CALL cp_cfm_to_cfm(oo_v, oo_vt)
         CALL cp_cfm_column_scale(oo_v, eigenvalues_sqrt)
         DEALLOCATE (eigenvalues)
         DEALLOCATE (eigenvalues_sqrt)
         CALL cp_cfm_gemm('N', 'C', nmo, nmo, nmo, (1.0_dp, 0.0_dp), &
                          oo_v, oo_vt, (0.0_dp, 0.0_dp), oo_c)
         oo_1%local_data = REAL(oo_c%local_data, KIND=dp)
         oo_2%local_data = AIMAG(oo_c%local_data)
         CALL cp_cfm_release(oo_c)
         CALL cp_cfm_release(oo_v)
         CALL cp_cfm_release(oo_vt)

         ! use this to compute the orthonormal vectors
         CALL cp_gemm("N", "N", nao, nmo, nmo, 1.0_dp, mos_old(2*ispin - 1)%matrix, oo_1, 0.0_dp, mos_new(2*ispin - 1)%matrix)
         CALL cp_gemm("N", "N", nao, nmo, nmo, 1.0_dp, mos_old(2*ispin - 1)%matrix, oo_2, 0.0_dp, mos_new(2*ispin)%matrix)

         CALL cp_gemm("N", "N", nao, nmo, nmo, -1.0_dp, mos_old(2*ispin)%matrix, oo_2, 0.0_dp, mos_old(2*ispin - 1)%matrix)
         CALL cp_fm_scale_and_add(1.0_dp, mos_old(2*ispin - 1)%matrix, 1.0_dp, mos_new(2*ispin - 1)%matrix)

         CALL cp_gemm("N", "N", nao, nmo, nmo, 1.0_dp, mos_old(2*ispin)%matrix, oo_1, 1.0_dp, mos_new(2*ispin)%matrix)
         CALL cp_fm_to_fm(mos_new(2*ispin)%matrix, mos_old(2*ispin)%matrix)

         CALL cp_fm_release(oo_1)
         CALL cp_fm_release(oo_2)
      END DO

      CALL cp_fm_release(S_chol)
      CALL cp_fm_release(mat_ks)
      CALL cp_fm_release(eigenvectors)

      IF (com_nl) CALL dbcsr_deallocate_matrix_set(matrix_rv)

!***************************************************************
!remove later
      CALL cp_fm_release(S_inv_fm)
      CALL cp_fm_release(tmpS)
!**************************************************************
      CALL timestop(handle)

   END SUBROUTINE apply_delta_pulse_periodic

! **************************************************************************************************
!> \brief applies exp(ikr) to the wavefunction.... stored in mos_old...
!> \param qs_env ...
!> \param mos_old ...
!> \param mos_new ...
!> \author Joost & Martin (2011)
! **************************************************************************************************

   SUBROUTINE apply_delta_pulse(qs_env, mos_old, mos_new)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_old, mos_new

      CHARACTER(len=*), PARAMETER                        :: routineN = 'apply_delta_pulse'

      COMPLEX(KIND=dp), DIMENSION(:), POINTER            :: eigenvalues_sqrt
      INTEGER                                            :: handle, i, nao, nmo
      REAL(KIND=dp), DIMENSION(3)                        :: kvec
      REAL(kind=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_cfm_type), POINTER                         :: oo_c, oo_v, oo_vt
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: mat_S, oo_1, oo_2, S_inv_fm, tmp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: cosmat, S_inv, sinmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(rt_prop_type), POINTER                        :: rtp

      NULLIFY (dft_control)

      CALL timeset(routineN, handle)

      ! we need the inverse overlap

      CALL get_qs_env(qs_env, &
                      mos=mos, &
                      rtp=rtp, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control)
      CALL get_rtp(rtp=rtp, S_inv=S_inv)

      CALL cp_fm_create(S_inv_fm, matrix_struct=rtp%ao_ao_fmstruct, name="tmp_mat")

      CALL cp_fm_create(tmp, matrix_struct=rtp%ao_ao_fmstruct, name="tmp_mat")

      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, S_inv_fm)
      CALL cp_fm_cholesky_decompose(S_inv_fm)
      CALL cp_fm_cholesky_invert(S_inv_fm)
      CALL cp_fm_upper_to_full(S_inv_fm, tmp)

      CALL cp_fm_create(mat_S, matrix_struct=S_inv_fm%matrix_struct, name="mat_S")
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, mat_S)
      CALL cp_fm_upper_to_full(mat_S, tmp)

      CALL cp_fm_release(tmp)

      ! we need the berry matrix
      CALL get_qs_env(qs_env, cell=cell)

      ! direction ... unscaled, this will yield a exp(ikr) that is periodic with the cell
      kvec(:) = cell%h_inv(1, :)*dft_control%rtp_control%delta_pulse_direction(1) + &
                cell%h_inv(2, :)*dft_control%rtp_control%delta_pulse_direction(2) + &
                cell%h_inv(3, :)*dft_control%rtp_control%delta_pulse_direction(3)
      kvec = -kvec*twopi
      ! scaling will make the things not periodic with the cell, which would only be good for gas phase systems ?
      kvec(:) = dft_control%rtp_control%delta_pulse_scale*kvec

      ALLOCATE (cosmat, sinmat)
      CALL dbcsr_copy(cosmat, matrix_s(1)%matrix, 'COS MOM')
      CALL dbcsr_copy(sinmat, matrix_s(1)%matrix, 'SIN MOM')
      CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec)

      ! apply inv(S)*operator to C
      DO i = 1, SIZE(mos)
         CALL get_mo_set(mos(i)%mo_set, nao=nao, nmo=nmo)
         CALL cp_dbcsr_sm_fm_multiply(cosmat, mos(i)%mo_set%mo_coeff, mos_new(2*i - 1)%matrix, ncol=nmo)
         CALL cp_dbcsr_sm_fm_multiply(sinmat, mos(i)%mo_set%mo_coeff, mos_new(2*i)%matrix, ncol=nmo)

         CALL cp_gemm("N", "N", nao, nmo, nao, 1.0_dp, S_inv_fm, mos_new(2*i - 1)%matrix, 0.0_dp, mos_old(2*i - 1)%matrix)
         CALL cp_gemm("N", "N", nao, nmo, nao, 1.0_dp, S_inv_fm, mos_new(2*i)%matrix, 0.0_dp, mos_old(2*i)%matrix)

         ! in a finite basis, unfortunately, inv(S)*operator is not unitary, so orthonormalize afterwards
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=S_inv_fm%matrix_struct%para_env, context=S_inv_fm%matrix_struct%context, &
                                  nrow_global=nmo, ncol_global=nmo)
         CALL cp_fm_create(oo_1, matrix_struct=fm_struct_tmp, name="oo_1")
         CALL cp_fm_create(oo_2, matrix_struct=fm_struct_tmp, name="oo_2")
         CALL cp_fm_struct_release(fm_struct_tmp)

         CALL cp_fm_create(tmp, matrix_struct=mos_old(2*i - 1)%matrix%matrix_struct, name="tmp_mat")
         ! get the complex overlap matrix
         ! x^T S x + y^T S y + i (-y^TS x+x^T S y)
         CALL cp_gemm("N", "N", nao, nmo, nao, 1.0_dp, mat_S, mos_old(2*i - 1)%matrix, 0.0_dp, tmp)
         CALL cp_gemm("T", "N", nmo, nmo, nao, 1.0_dp, mos_old(2*i - 1)%matrix, tmp, 0.0_dp, oo_1)
         CALL cp_gemm("T", "N", nmo, nmo, nao, -1.0_dp, mos_old(2*i)%matrix, tmp, 0.0_dp, oo_2)

         CALL cp_gemm("N", "N", nao, nmo, nao, 1.0_dp, mat_S, mos_old(2*i)%matrix, 0.0_dp, tmp)
         CALL cp_gemm("T", "N", nmo, nmo, nao, 1.0_dp, mos_old(2*i)%matrix, tmp, 1.0_dp, oo_1)
         CALL cp_gemm("T", "N", nmo, nmo, nao, 1.0_dp, mos_old(2*i - 1)%matrix, tmp, 1.0_dp, oo_2)
         CALL cp_fm_release(tmp)

         CALL cp_cfm_create(oo_c, oo_1%matrix_struct)
         CALL cp_cfm_create(oo_v, oo_1%matrix_struct)
         CALL cp_cfm_create(oo_vt, oo_1%matrix_struct)
         oo_c%local_data = CMPLX(oo_1%local_data, oo_2%local_data, KIND=dp)

         ! compute inv(sqrt(overlap))
         ALLOCATE (eigenvalues(nmo))
         ALLOCATE (eigenvalues_sqrt(nmo))
         CALL cp_cfm_heevd(oo_c, oo_v, eigenvalues)
         eigenvalues_sqrt = CMPLX(1.0_dp/SQRT(eigenvalues), 0.0_dp, dp)
         CALL cp_cfm_to_cfm(oo_v, oo_vt)
         CALL cp_cfm_column_scale(oo_v, eigenvalues_sqrt)
         DEALLOCATE (eigenvalues)
         DEALLOCATE (eigenvalues_sqrt)
         CALL cp_cfm_gemm('N', 'C', nmo, nmo, nmo, (1.0_dp, 0.0_dp), &
                          oo_v, oo_vt, (0.0_dp, 0.0_dp), oo_c)
         oo_1%local_data = REAL(oo_c%local_data, KIND=dp)
         oo_2%local_data = AIMAG(oo_c%local_data)
         CALL cp_cfm_release(oo_c)
         CALL cp_cfm_release(oo_v)
         CALL cp_cfm_release(oo_vt)

         ! use this to compute the orthonormal vectors
         CALL cp_gemm("N", "N", nao, nmo, nmo, 1.0_dp, mos_old(2*i - 1)%matrix, oo_1, 0.0_dp, mos_new(2*i - 1)%matrix)
         CALL cp_gemm("N", "N", nao, nmo, nmo, 1.0_dp, mos_old(2*i - 1)%matrix, oo_2, 0.0_dp, mos_new(2*i)%matrix)

         CALL cp_gemm("N", "N", nao, nmo, nmo, -1.0_dp, mos_old(2*i)%matrix, oo_2, 0.0_dp, mos_old(2*i - 1)%matrix)
         CALL cp_fm_scale_and_add(1.0_dp, mos_old(2*i - 1)%matrix, 1.0_dp, mos_new(2*i - 1)%matrix)

         CALL cp_gemm("N", "N", nao, nmo, nmo, 1.0_dp, mos_old(2*i)%matrix, oo_1, 1.0_dp, mos_new(2*i)%matrix)
         CALL cp_fm_to_fm(mos_new(2*i)%matrix, mos_old(2*i)%matrix)

         CALL cp_fm_release(oo_1)
         CALL cp_fm_release(oo_2)
      END DO

      CALL cp_fm_release(mat_S)

      CALL dbcsr_deallocate_matrix(cosmat)
      CALL dbcsr_deallocate_matrix(sinmat)

!***************************************************************
!remove later
      CALL copy_fm_to_dbcsr(S_inv_fm, S_inv)
      CALL dbcsr_filter(S_inv, rtp%filter_eps)
      CALL cp_fm_release(S_inv_fm)
!**************************************************************

      CALL timestop(handle)

   END SUBROUTINE apply_delta_pulse

END MODULE rt_delta_pulse
