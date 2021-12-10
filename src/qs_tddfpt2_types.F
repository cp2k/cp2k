!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_types
   USE admm_types,                      ONLY: admm_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type,&
                                              get_blacs_info
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pool_create,&
                                              fm_pool_create_fm,&
                                              fm_pool_release
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_p_type,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_distribution_type,&
                                              dbcsr_get_info,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_release_p,&
                                              dbcsr_type
   USE ewald_environment_types,         ONLY: ewald_env_release,&
                                              ewald_environment_type
   USE ewald_pw_types,                  ONLY: ewald_pw_release,&
                                              ewald_pw_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_dbcsr_create_by_dist,&
                                              tddfpt_subgroup_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_types'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   ! number of first derivative components (3: d/dx, d/dy, d/dz)
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC :: tddfpt_ground_state_mos, tddfpt_work_matrices
   PUBLIC :: tddfpt_create_work_matrices, stda_create_work_matrices, tddfpt_release_work_matrices

! **************************************************************************************************
!> \brief Ground state molecular orbitals.
!> \par History
!>   * 06.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_ground_state_mos
      !> occupied MOs stored in a matrix form [nao x nmo_occ]
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      !> virtual MOs stored in a matrix form [nao x nmo_virt]
      TYPE(cp_fm_type), POINTER                          :: mos_virt
      !> negated occupied orbital energy matrix [nmo_occ x nmo_occ]: - mos_occ^T * KS * mos_occ .
      !> Allocated when orbital energy correction is in use, otherwise it is just a diagonal
      !> matrix with 'evals_occ' on its diagonal
      TYPE(cp_fm_type), POINTER                          :: evals_occ_matrix
      !> (non-corrected) occupied orbital energies
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals_occ
      !> (non-corrected) virtual orbital energies
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals_virt
      !> phase of occupied MOs; +1.0 -- positive, -1.0 -- negative;
      !> it is mainly needed to make the restart file transferable
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: phases_occ
   END TYPE tddfpt_ground_state_mos

! **************************************************************************************************
!> \brief Set of temporary ("work") matrices.
!> \par History
!>   * 01.2017 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE tddfpt_work_matrices
      !
      ! *** globally distributed dense matrices ***
      !
      !> pool of dense [nao x nmo_occ(spin)] matrices;
      !> used mainly to dynamically expand the list of trial vectors
      TYPE(cp_fm_pool_p_type), ALLOCATABLE, DIMENSION(:) :: fm_pool_ao_mo_occ
      !> S * mos_occ(spin)
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: S_C0
      !> S * \rho_0(spin)
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: S_C0_C0T
      !
      ! *** dense matrices distributed across parallel (sub)groups ***
      !
      !> evects_sub(1:nspins, 1:nstates): a copy of the last 'nstates' trial vectors distributed
      !> across parallel (sub)groups. Here 'nstates' is the number of requested excited states which
      !> is typically much smaller than the total number of Krylov's vectors. Allocated only if
      !> the number of parallel groups > 1, otherwise we use the original globally distributed vectors.
      !> evects_sub(spin, state) == null() means that the trial vector is assigned to a different (sub)group
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: evects_sub
      !> action of TDDFPT operator on trial vectors distributed across parallel (sub)groups
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: Aop_evects_sub
      !> electron density expressed in terms of atomic orbitals using primary basis set
      TYPE(cp_fm_type), POINTER                          :: rho_ao_orb_fm_sub
      !
      ! NOTE: we do not need the next 2 matrices in case of a sparse matrix 'tddfpt_subgroup_env_type%admm_A'
      !
      !> electron density expressed in terms of atomic orbitals using auxiliary basis set;
      !> can be seen as a group-specific version of the matrix 'admm_type%work_aux_aux'
      TYPE(cp_fm_type), POINTER                          :: rho_ao_aux_fit_fm_sub
      !> group-specific version of the matrix 'admm_type%work_aux_orb' with shape [nao_aux x nao]
      TYPE(cp_fm_type), POINTER                          :: wfm_aux_orb_sub
      !
      ! *** sparse matrices distributed across parallel (sub)groups ***
      !
      !> sparse matrix with shape [nao x nao] distributed across subgroups;
      !> Aop_evects_sub(spin,:) = A_ia_munu_sub(spin) * mos_occ(spin)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: A_ia_munu_sub
      !
      ! *** structures to store electron densities distributed across parallel (sub)groups ***
      !
      !> electron density in terms of primary basis set
      TYPE(qs_rho_type), POINTER                         :: rho_orb_struct_sub
      !> electron density in terms of auxiliary basis set
      TYPE(qs_rho_type), POINTER                         :: rho_aux_fit_struct_sub
      !> group-specific copy of a Coulomb/xc-potential on a real-space grid
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: A_ia_rspace_sub
      !> group-specific copy of a reciprocal-space grid
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: wpw_gspace_sub
      !> group-specific copy of a real-space grid
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: wpw_rspace_sub
      !
      ! *** globally distributed matrices required to compute exact exchange terms ***
      !
      !> globally distributed version of the matrix 'rho_ao_orb_fm_sub' to store the electron density
      TYPE(cp_fm_type), POINTER                          :: hfx_fm_ao_ao
      !> sparse matrix to store the electron density in terms of auxiliary (ADMM calculation)
      !> or primary (regular calculation) basis set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: hfx_rho_ao
      !> exact exchange expressed in terms of auxiliary or primary basis set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: hfx_hmat
      !
      ! *** matrices required for sTDA kernel, all matrices are within subgroups
      !
      ! Short-range gamma exchange matrix
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gamma_exchange
      !Lowdin MO coefficients: NAO*NOCC
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: ctransformed
      !S^1/2
      TYPE(dbcsr_type), POINTER                          :: shalf
      !Eigenvalues/eigenvectors of the overlap matrix, used in sTDA forces (Lowdin derivatives)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: S_eigenvalues
      TYPE(cp_fm_type), POINTER                          :: S_eigenvectors
      TYPE(cp_fm_type), POINTER                          :: slambda
      !Ewald environments
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
   END TYPE tddfpt_work_matrices

CONTAINS

! **************************************************************************************************
!> \brief Allocate work matrices for full kernel
!> \param work_matrices  work matrices (allocated on exit)
!> \param gs_mos         occupied and virtual molecular orbitals optimised for the ground state
!> \param nstates        number of excited states to converge
!> \param do_hfx         flag that requested to allocate work matrices required for computation
!>                       of exact-exchange terms
!> \param do_admm ...
!> \param qs_env         Quickstep environment
!> \param sub_env        parallel group environment
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_create_work_matrices(work_matrices, gs_mos, nstates, do_hfx, do_admm, &
                                          qs_env, sub_env)
      TYPE(tddfpt_work_matrices), INTENT(out)            :: work_matrices
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      INTEGER, INTENT(in)                                :: nstates
      LOGICAL, INTENT(in)                                :: do_hfx, do_admm
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_create_work_matrices'

      INTEGER                                            :: handle, igroup, ispin, istate, nao, &
                                                            nao_aux, ngroups, nspins
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ, nmo_virt
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_p_type), DIMENSION(maxspins)     :: fm_struct_evects
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_s_aux_fit, rho_ia_ao
      TYPE(dbcsr_type), POINTER                          :: dbcsr_template_hfx
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_hfx
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      ! sTDA
      NULLIFY (work_matrices%shalf)
      NULLIFY (work_matrices%ewald_env)
      NULLIFY (work_matrices%ewald_pw)
      NULLIFY (work_matrices%gamma_exchange)
      NULLIFY (work_matrices%ctransformed)
      NULLIFY (work_matrices%S_eigenvalues)
      NULLIFY (work_matrices%S_eigenvectors)
      NULLIFY (work_matrices%slambda)

      nspins = SIZE(gs_mos)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)

      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
      END DO

      IF (do_admm) THEN
         CPASSERT(do_hfx)
         CPASSERT(ASSOCIATED(sub_env%admm_A))
         CALL get_qs_env(qs_env, matrix_s_aux_fit=matrix_s_aux_fit)
         CALL dbcsr_get_info(matrix_s_aux_fit(1)%matrix, nfullrows_total=nao_aux)
      END IF

      NULLIFY (fm_struct)
      ALLOCATE (work_matrices%fm_pool_ao_mo_occ(nspins))
      DO ispin = 1, nspins
         NULLIFY (work_matrices%fm_pool_ao_mo_occ(ispin)%pool)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo_occ(ispin), context=blacs_env)
         CALL fm_pool_create(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
      END DO

      ALLOCATE (work_matrices%S_C0_C0T(nspins))
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)
      DO ispin = 1, nspins
         NULLIFY (work_matrices%S_C0_C0T(ispin)%matrix)
         CALL cp_fm_create(work_matrices%S_C0_C0T(ispin)%matrix, fm_struct)
      END DO
      CALL cp_fm_struct_release(fm_struct)

      ALLOCATE (work_matrices%S_C0(nspins))
      DO ispin = 1, nspins
         NULLIFY (work_matrices%S_C0(ispin)%matrix)
         CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, work_matrices%S_C0(ispin)%matrix)

         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, gs_mos(ispin)%mos_occ, work_matrices%S_C0(ispin)%matrix, &
                                      ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
         CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 1.0_dp, work_matrices%S_C0(ispin)%matrix, &
                      gs_mos(ispin)%mos_occ, 0.0_dp, work_matrices%S_C0_C0T(ispin)%matrix)
      END DO

      IF (sub_env%is_split) THEN
         DO ispin = 1, nspins
            NULLIFY (fm_struct_evects(ispin)%struct)
            CALL cp_fm_struct_create(fm_struct_evects(ispin)%struct, nrow_global=nao, &
                                     ncol_global=nmo_occ(ispin), context=sub_env%blacs_env)
         END DO

         ALLOCATE (work_matrices%evects_sub(nspins, nstates), work_matrices%Aop_evects_sub(nspins, nstates))
         DO istate = 1, nstates
            DO ispin = 1, nspins
               NULLIFY (work_matrices%evects_sub(ispin, istate)%matrix)
               NULLIFY (work_matrices%Aop_evects_sub(ispin, istate)%matrix)
            END DO
         END DO

         CALL get_blacs_info(blacs_env, para_env=para_env)
         igroup = sub_env%group_distribution(para_env%mepos)
         ngroups = sub_env%ngroups

         DO istate = ngroups - igroup, nstates, ngroups
            DO ispin = 1, nspins
               CALL cp_fm_create(work_matrices%evects_sub(ispin, istate)%matrix, fm_struct_evects(ispin)%struct)
               CALL cp_fm_create(work_matrices%Aop_evects_sub(ispin, istate)%matrix, fm_struct_evects(ispin)%struct)
            END DO
         END DO

         DO ispin = nspins, 1, -1
            CALL cp_fm_struct_release(fm_struct_evects(ispin)%struct)
         END DO
      END IF

      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=sub_env%blacs_env)
      NULLIFY (work_matrices%rho_ao_orb_fm_sub)
      CALL cp_fm_create(work_matrices%rho_ao_orb_fm_sub, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      NULLIFY (work_matrices%rho_ao_aux_fit_fm_sub, work_matrices%wfm_aux_orb_sub)
      IF (do_admm) THEN
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao_aux, ncol_global=nao_aux, context=sub_env%blacs_env)
         CALL cp_fm_create(work_matrices%rho_ao_aux_fit_fm_sub, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         CALL cp_fm_struct_create(fm_struct, nrow_global=nao_aux, ncol_global=nao, context=sub_env%blacs_env)
         CALL cp_fm_create(work_matrices%wfm_aux_orb_sub, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
      END IF

      ! group-specific dbcsr matrices
      NULLIFY (work_matrices%A_ia_munu_sub)
      CALL dbcsr_allocate_matrix_set(work_matrices%A_ia_munu_sub, nspins)
      DO ispin = 1, nspins
         CALL dbcsr_init_p(work_matrices%A_ia_munu_sub(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(work_matrices%A_ia_munu_sub(ispin)%matrix, template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, sab=sub_env%sab_orb)
      END DO

      ! group-specific response density
      NULLIFY (rho_ia_ao)
      CALL dbcsr_allocate_matrix_set(rho_ia_ao, nspins)
      DO ispin = 1, nspins
         CALL dbcsr_init_p(rho_ia_ao(ispin)%matrix)
         CALL tddfpt_dbcsr_create_by_dist(rho_ia_ao(ispin)%matrix, template=matrix_s(1)%matrix, &
                                          dbcsr_dist=sub_env%dbcsr_dist, sab=sub_env%sab_orb)
      END DO

      NULLIFY (work_matrices%rho_orb_struct_sub)
      CALL qs_rho_create(work_matrices%rho_orb_struct_sub)
      CALL qs_rho_set(work_matrices%rho_orb_struct_sub, rho_ao=rho_ia_ao)
      CALL qs_rho_rebuild(work_matrices%rho_orb_struct_sub, qs_env, rebuild_ao=.FALSE., &
                          rebuild_grids=.TRUE., pw_env_external=sub_env%pw_env)

      NULLIFY (work_matrices%rho_aux_fit_struct_sub)
      IF (do_admm) THEN
         NULLIFY (rho_ia_ao)
         CALL dbcsr_allocate_matrix_set(rho_ia_ao, nspins)
         DO ispin = 1, nspins
            CALL dbcsr_init_p(rho_ia_ao(ispin)%matrix)
            CALL tddfpt_dbcsr_create_by_dist(rho_ia_ao(ispin)%matrix, template=matrix_s_aux_fit(1)%matrix, &
                                             dbcsr_dist=sub_env%dbcsr_dist, sab=sub_env%sab_aux_fit)
         END DO

         CALL qs_rho_create(work_matrices%rho_aux_fit_struct_sub)
         CALL qs_rho_set(work_matrices%rho_aux_fit_struct_sub, rho_ao=rho_ia_ao)
         CALL qs_rho_rebuild(work_matrices%rho_aux_fit_struct_sub, qs_env, rebuild_ao=.FALSE., &
                             rebuild_grids=.TRUE., pw_env_external=sub_env%pw_env)
      END IF

      ! work plain-wave grids
      CALL pw_env_get(sub_env%pw_env, auxbas_pw_pool=auxbas_pw_pool)
      ALLOCATE (work_matrices%A_ia_rspace_sub(nspins))
      ALLOCATE (work_matrices%wpw_gspace_sub(nspins), work_matrices%wpw_rspace_sub(nspins))
      DO ispin = 1, nspins
         NULLIFY (work_matrices%A_ia_rspace_sub(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, work_matrices%A_ia_rspace_sub(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)

         NULLIFY (work_matrices%wpw_gspace_sub(ispin)%pw, work_matrices%wpw_rspace_sub(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, work_matrices%wpw_gspace_sub(ispin)%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, work_matrices%wpw_rspace_sub(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
      END DO

      ! HFX-related globally distributed matrices
      NULLIFY (work_matrices%hfx_fm_ao_ao, work_matrices%hfx_rho_ao, work_matrices%hfx_hmat)
      IF (do_hfx) THEN
         IF (do_admm) THEN
            CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist, sab_aux_fit=sab_hfx)
            dbcsr_template_hfx => matrix_s_aux_fit(1)%matrix
         ELSE
            CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist, sab_orb=sab_hfx)
            dbcsr_template_hfx => matrix_s(1)%matrix
         END IF

         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)
         CALL cp_fm_create(work_matrices%hfx_fm_ao_ao, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         CALL dbcsr_allocate_matrix_set(work_matrices%hfx_rho_ao, nspins)
         DO ispin = 1, nspins
            CALL dbcsr_init_p(work_matrices%hfx_rho_ao(ispin)%matrix)
            CALL tddfpt_dbcsr_create_by_dist(work_matrices%hfx_rho_ao(ispin)%matrix, &
                                             template=dbcsr_template_hfx, dbcsr_dist=dbcsr_dist, sab=sab_hfx)
         END DO

         CALL dbcsr_allocate_matrix_set(work_matrices%hfx_hmat, nspins)
         DO ispin = 1, nspins
            CALL dbcsr_init_p(work_matrices%hfx_hmat(ispin)%matrix)
            CALL tddfpt_dbcsr_create_by_dist(work_matrices%hfx_hmat(ispin)%matrix, &
                                             template=dbcsr_template_hfx, dbcsr_dist=dbcsr_dist, sab=sab_hfx)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_create_work_matrices

! **************************************************************************************************
!> \brief Allocate work matrices for sTDA kernel
!> \param work_matrices  work matrices (allocated on exit)
!> \param gs_mos         occupied and virtual molecular orbitals optimised for the ground state
!> \param nstates        number of excited states to converge
!> \param qs_env         Quickstep environment
!> \param sub_env        parallel group environment
!> \par History
!>    * 04.2019 created from full kernel version [JHU]
! **************************************************************************************************
   SUBROUTINE stda_create_work_matrices(work_matrices, gs_mos, nstates, qs_env, sub_env)
      TYPE(tddfpt_work_matrices), INTENT(out)            :: work_matrices
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      INTEGER, INTENT(in)                                :: nstates
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'stda_create_work_matrices'

      INTEGER                                            :: handle, igroup, ispin, istate, nao, &
                                                            ngroups, nspins
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ, nmo_virt
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_p_type), DIMENSION(maxspins)     :: fm_struct_evects
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s

      CALL timeset(routineN, handle)

      NULLIFY (work_matrices%gamma_exchange, work_matrices%ctransformed)

      nspins = SIZE(gs_mos)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)

      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
      END DO

      NULLIFY (fm_struct)
      ALLOCATE (work_matrices%fm_pool_ao_mo_occ(nspins))
      DO ispin = 1, nspins
         NULLIFY (work_matrices%fm_pool_ao_mo_occ(ispin)%pool)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo_occ(ispin), context=blacs_env)
         CALL fm_pool_create(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
      END DO

      ALLOCATE (work_matrices%S_C0_C0T(nspins))
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)
      DO ispin = 1, nspins
         NULLIFY (work_matrices%S_C0_C0T(ispin)%matrix)
         CALL cp_fm_create(work_matrices%S_C0_C0T(ispin)%matrix, fm_struct)
      END DO
      CALL cp_fm_struct_release(fm_struct)

      ALLOCATE (work_matrices%S_C0(nspins))
      DO ispin = 1, nspins
         NULLIFY (work_matrices%S_C0(ispin)%matrix)
         CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, work_matrices%S_C0(ispin)%matrix)

         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, gs_mos(ispin)%mos_occ, work_matrices%S_C0(ispin)%matrix, &
                                      ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
         CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 1.0_dp, work_matrices%S_C0(ispin)%matrix, &
                      gs_mos(ispin)%mos_occ, 0.0_dp, work_matrices%S_C0_C0T(ispin)%matrix)
      END DO

      DO ispin = 1, nspins
         NULLIFY (fm_struct_evects(ispin)%struct)
         CALL cp_fm_struct_create(fm_struct_evects(ispin)%struct, nrow_global=nao, &
                                  ncol_global=nmo_occ(ispin), context=sub_env%blacs_env)
      END DO

      IF (sub_env%is_split) THEN
         ALLOCATE (work_matrices%evects_sub(nspins, nstates), work_matrices%Aop_evects_sub(nspins, nstates))
         DO istate = 1, nstates
            DO ispin = 1, nspins
               NULLIFY (work_matrices%evects_sub(ispin, istate)%matrix)
               NULLIFY (work_matrices%Aop_evects_sub(ispin, istate)%matrix)
            END DO
         END DO

         CALL get_blacs_info(blacs_env, para_env=para_env)
         igroup = sub_env%group_distribution(para_env%mepos)
         ngroups = sub_env%ngroups

         DO istate = ngroups - igroup, nstates, ngroups
            DO ispin = 1, nspins
               CALL cp_fm_create(work_matrices%evects_sub(ispin, istate)%matrix, fm_struct_evects(ispin)%struct)
               CALL cp_fm_create(work_matrices%Aop_evects_sub(ispin, istate)%matrix, fm_struct_evects(ispin)%struct)
            END DO
         END DO
      END IF

      ! sTDA specific work arrays
      ALLOCATE (work_matrices%ctransformed(nspins))
      DO ispin = 1, nspins
         NULLIFY (work_matrices%ctransformed(ispin)%matrix)
         CALL cp_fm_create(work_matrices%ctransformed(ispin)%matrix, fm_struct_evects(ispin)%struct)
      END DO
      NULLIFY (work_matrices%shalf)
      CALL dbcsr_init_p(work_matrices%shalf)
      CALL dbcsr_create(work_matrices%shalf, template=matrix_s(1)%matrix)
      ! forces
      ALLOCATE (work_matrices%S_eigenvalues(nao))
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)
      CALL cp_fm_create(work_matrices%S_eigenvectors, fm_struct)
      CALL cp_fm_create(work_matrices%slambda, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      DO ispin = nspins, 1, -1
         CALL cp_fm_struct_release(fm_struct_evects(ispin)%struct)
      END DO

      NULLIFY (work_matrices%rho_ao_orb_fm_sub)
      NULLIFY (work_matrices%rho_ao_aux_fit_fm_sub, work_matrices%wfm_aux_orb_sub)
      NULLIFY (work_matrices%rho_aux_fit_struct_sub)
      NULLIFY (work_matrices%rho_orb_struct_sub)
      NULLIFY (work_matrices%hfx_fm_ao_ao, work_matrices%hfx_rho_ao, work_matrices%hfx_hmat)
      NULLIFY (work_matrices%A_ia_rspace_sub, work_matrices%wpw_gspace_sub, &
               work_matrices%wpw_rspace_sub)
      NULLIFY (work_matrices%A_ia_munu_sub)

      NULLIFY (work_matrices%ewald_env)
      NULLIFY (work_matrices%ewald_pw)

      CALL timestop(handle)

   END SUBROUTINE stda_create_work_matrices

! **************************************************************************************************
!> \brief Release work matrices.
!> \param work_matrices  work matrices (destroyed on exit)
!> \param sub_env        parallel group environment
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_release_work_matrices(work_matrices, sub_env)
      TYPE(tddfpt_work_matrices), INTENT(inout)          :: work_matrices
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_release_work_matrices'

      INTEGER                                            :: handle, ispin, istate
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      ! HFX-related matrices
      IF (ASSOCIATED(work_matrices%hfx_hmat)) THEN
         DO ispin = SIZE(work_matrices%hfx_hmat), 1, -1
            CALL dbcsr_deallocate_matrix(work_matrices%hfx_hmat(ispin)%matrix)
         END DO
         DEALLOCATE (work_matrices%hfx_hmat)
      END IF

      IF (ASSOCIATED(work_matrices%hfx_rho_ao)) THEN
         DO ispin = SIZE(work_matrices%hfx_rho_ao), 1, -1
            CALL dbcsr_deallocate_matrix(work_matrices%hfx_rho_ao(ispin)%matrix)
         END DO
         DEALLOCATE (work_matrices%hfx_rho_ao)
      END IF

      IF (ASSOCIATED(work_matrices%hfx_fm_ao_ao)) &
         CALL cp_fm_release(work_matrices%hfx_fm_ao_ao)

      ! real-space and reciprocal-space grids
      IF (ASSOCIATED(sub_env%pw_env)) THEN
         CALL pw_env_get(sub_env%pw_env, auxbas_pw_pool=auxbas_pw_pool)
         DO ispin = SIZE(work_matrices%wpw_rspace_sub), 1, -1
            CALL pw_pool_give_back_pw(auxbas_pw_pool, work_matrices%wpw_rspace_sub(ispin)%pw)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, work_matrices%wpw_gspace_sub(ispin)%pw)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, work_matrices%A_ia_rspace_sub(ispin)%pw)
         END DO
         DEALLOCATE (work_matrices%A_ia_rspace_sub, work_matrices%wpw_gspace_sub, work_matrices%wpw_rspace_sub)
      END IF

      IF (ASSOCIATED(work_matrices%rho_aux_fit_struct_sub)) &
         CALL qs_rho_release(work_matrices%rho_aux_fit_struct_sub)
      CALL qs_rho_release(work_matrices%rho_orb_struct_sub)

      IF (ASSOCIATED(work_matrices%A_ia_munu_sub)) THEN
         DO ispin = SIZE(work_matrices%A_ia_munu_sub), 1, -1
            CALL dbcsr_deallocate_matrix(work_matrices%A_ia_munu_sub(ispin)%matrix)
         END DO
         DEALLOCATE (work_matrices%A_ia_munu_sub)
      END IF

      IF (ASSOCIATED(work_matrices%wfm_aux_orb_sub)) &
         CALL cp_fm_release(work_matrices%wfm_aux_orb_sub)
      IF (ASSOCIATED(work_matrices%rho_ao_aux_fit_fm_sub)) &
         CALL cp_fm_release(work_matrices%rho_ao_aux_fit_fm_sub)
      CALL cp_fm_release(work_matrices%rho_ao_orb_fm_sub)

      IF (ALLOCATED(work_matrices%evects_sub)) THEN
         DO istate = SIZE(work_matrices%evects_sub, 2), 1, -1
            DO ispin = SIZE(work_matrices%evects_sub, 1), 1, -1
               CALL cp_fm_release(work_matrices%Aop_evects_sub(ispin, istate)%matrix)
               CALL cp_fm_release(work_matrices%evects_sub(ispin, istate)%matrix)
            END DO
         END DO
         DEALLOCATE (work_matrices%Aop_evects_sub, work_matrices%evects_sub)
      END IF

      DO ispin = SIZE(work_matrices%fm_pool_ao_mo_occ), 1, -1
         CALL cp_fm_release(work_matrices%S_C0(ispin)%matrix)
         CALL cp_fm_release(work_matrices%S_C0_C0T(ispin)%matrix)
      END DO
      DEALLOCATE (work_matrices%S_C0, work_matrices%S_C0_C0T)

      DO ispin = SIZE(work_matrices%fm_pool_ao_mo_occ), 1, -1
         CALL fm_pool_release(work_matrices%fm_pool_ao_mo_occ(ispin)%pool)
      END DO
      DEALLOCATE (work_matrices%fm_pool_ao_mo_occ)

      ! sTDA
      IF (ASSOCIATED(work_matrices%gamma_exchange)) THEN
         CALL dbcsr_deallocate_matrix_set(work_matrices%gamma_exchange)
         NULLIFY (work_matrices%gamma_exchange)
      END IF
      IF (ASSOCIATED(work_matrices%ctransformed)) THEN
         CALL cp_fm_vect_dealloc(work_matrices%ctransformed)
         NULLIFY (work_matrices%ctransformed)
      END IF
      CALL dbcsr_release_p(work_matrices%shalf)
      !
      IF (ASSOCIATED(work_matrices%S_eigenvectors)) &
         CALL cp_fm_release(work_matrices%S_eigenvectors)
      IF (ASSOCIATED(work_matrices%slambda)) &
         CALL cp_fm_release(work_matrices%slambda)
      IF (ASSOCIATED(work_matrices%S_eigenvalues)) &
         DEALLOCATE (work_matrices%S_eigenvalues)
      ! Ewald
      CALL ewald_env_release(work_matrices%ewald_env)
      CALL ewald_pw_release(work_matrices%ewald_pw)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_release_work_matrices

END MODULE qs_tddfpt2_types
