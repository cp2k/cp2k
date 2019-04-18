!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
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
                                              cp_fm_get_info,&
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
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
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
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_tddfpt2_stda_types,           ONLY: stda_env_type
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_dbcsr_create_by_dist,&
                                              tddfpt_subgroup_env_type
   USE xc,                              ONLY: xc_prep_2nd_deriv
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_release
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_types'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   ! number of first derivative components (3: d/dx, d/dy, d/dz)
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC :: tddfpt_ground_state_mos, tddfpt_work_matrices, tddfpt_kernel_env_type, &
             full_kernel_env_type
   PUBLIC :: tddfpt_create_work_matrices, stda_create_work_matrices, tddfpt_release_work_matrices, &
             tddfpt_create_kernel_env, tddfpt_release_kernel_env

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
      !> globally distributed dense matrices with shape [nao x nmo_occ(spin)]
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: wfm_ao_mo_occ
      !> globally distributed dense matrices with shape [nmo_occ(spin) x nmo_occ(spin)]
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: wfm_mo_occ_mo_occ
      !> globally distributed dense matrices with shape [nmo_virt(spin) x nmo_occ(spin)]
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: wfm_mo_virt_mo_occ
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
      !Ewald environments
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
   END TYPE tddfpt_work_matrices

! **************************************************************************************************
!> \brief Collection of variables required to evaluate adiabatic TDDFPT kernel.
!> \par History
!>   * 12.2016 created [Sergey Chulkov]
! **************************************************************************************************
   TYPE full_kernel_env_type
      ! ground state electron density
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho_set
      ! response density
      TYPE(xc_rho_set_type), POINTER                     :: xc_rho1_set
      !> first and second derivatives of exchange-correlation functional
      TYPE(xc_derivative_set_type), POINTER              :: xc_deriv_set
      !> XC input section
      TYPE(section_vals_type), POINTER                   :: xc_section
      !> flags which indicate required components of the exchange-correlation functional
      !> (density, gradient, etc)
      TYPE(xc_rho_cflags_type)                           :: xc_rho1_cflags
      !> the method used to compute derivatives
      INTEGER                                            :: deriv_method_id
      !> the density smoothing method
      INTEGER                                            :: rho_smooth_id
      !> scaling coefficients in the linear combination:
      !> K = alpha * K_{\alpha,\alpha} + beta * K_{\alpha,\beta}
      REAL(kind=dp)                                      :: alpha, beta
   END TYPE full_kernel_env_type
! **************************************************************************************************
!> \brief Type to hold environments for the different kernels
!> \par History
!>   * 04.2019 created [JHU]
! **************************************************************************************************
   TYPE tddfpt_kernel_env_type
      TYPE(full_kernel_env_type), POINTER                :: full_kernel => Null()
      TYPE(full_kernel_env_type), POINTER                :: admm_kernel => Null()
      TYPE(stda_env_type), POINTER                       :: stda_kernel => Null()
   END TYPE tddfpt_kernel_env_type

CONTAINS

! **************************************************************************************************
!> \brief Allocate work matrices for full kernel
!> \param work_matrices  work matrices (allocated on exit)
!> \param gs_mos         occupied and virtual molecular orbitals optimised for the ground state
!> \param nstates        number of excited states to converge
!> \param do_hfx         flag that requested to allocate work matrices required for computation
!>                       of exact-exchange terms
!> \param qs_env         Quickstep environment
!> \param sub_env        parallel group environment
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_create_work_matrices(work_matrices, gs_mos, nstates, do_hfx, qs_env, sub_env)
      TYPE(tddfpt_work_matrices), INTENT(out)            :: work_matrices
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      INTEGER, INTENT(in)                                :: nstates
      LOGICAL, INTENT(in)                                :: do_hfx
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_create_work_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, igroup, ispin, istate, nao, &
                                                            nao_aux, ngroups, nspins
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ, nmo_virt
      LOGICAL                                            :: do_admm
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

      NULLIFY (work_matrices%gamma_exchange, work_matrices%ctransformed)

      nspins = SIZE(gs_mos)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)

      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
      END DO

      do_admm = do_hfx .AND. ASSOCIATED(sub_env%admm_A)
      IF (do_admm) THEN
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

      ALLOCATE (work_matrices%S_C0(nspins), work_matrices%wfm_ao_mo_occ(nspins))
      ALLOCATE (work_matrices%wfm_mo_occ_mo_occ(nspins), work_matrices%wfm_mo_virt_mo_occ(nspins))
      DO ispin = 1, nspins
         NULLIFY (work_matrices%S_C0(ispin)%matrix, work_matrices%wfm_ao_mo_occ(ispin)%matrix)
         CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, work_matrices%S_C0(ispin)%matrix)
         CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, work_matrices%wfm_ao_mo_occ(ispin)%matrix)

         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, gs_mos(ispin)%mos_occ, work_matrices%S_C0(ispin)%matrix, &
                                      ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
         CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 1.0_dp, work_matrices%S_C0(ispin)%matrix, &
                      gs_mos(ispin)%mos_occ, 0.0_dp, work_matrices%S_C0_C0T(ispin)%matrix)

         NULLIFY (work_matrices%wfm_mo_occ_mo_occ(ispin)%matrix)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_occ(ispin), ncol_global=nmo_occ(ispin), context=blacs_env)
         CALL cp_fm_create(work_matrices%wfm_mo_occ_mo_occ(ispin)%matrix, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         NULLIFY (work_matrices%wfm_mo_virt_mo_occ(ispin)%matrix)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_virt(ispin), ncol_global=nmo_occ(ispin), context=blacs_env)
         CALL cp_fm_create(work_matrices%wfm_mo_virt_mo_occ(ispin)%matrix, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
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

         DO istate = ngroups-igroup, nstates, ngroups
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'stda_create_work_matrices', &
         routineP = moduleN//':'//routineN

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

      ALLOCATE (work_matrices%S_C0(nspins), work_matrices%wfm_ao_mo_occ(nspins))
      ALLOCATE (work_matrices%wfm_mo_occ_mo_occ(nspins), work_matrices%wfm_mo_virt_mo_occ(nspins))
      DO ispin = 1, nspins
         NULLIFY (work_matrices%S_C0(ispin)%matrix, work_matrices%wfm_ao_mo_occ(ispin)%matrix)
         CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, work_matrices%S_C0(ispin)%matrix)
         CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, work_matrices%wfm_ao_mo_occ(ispin)%matrix)

         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, gs_mos(ispin)%mos_occ, work_matrices%S_C0(ispin)%matrix, &
                                      ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
         CALL cp_gemm('N', 'T', nao, nao, nmo_occ(ispin), 1.0_dp, work_matrices%S_C0(ispin)%matrix, &
                      gs_mos(ispin)%mos_occ, 0.0_dp, work_matrices%S_C0_C0T(ispin)%matrix)

         NULLIFY (work_matrices%wfm_mo_occ_mo_occ(ispin)%matrix)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_occ(ispin), ncol_global=nmo_occ(ispin), &
                                  context=blacs_env)
         CALL cp_fm_create(work_matrices%wfm_mo_occ_mo_occ(ispin)%matrix, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         NULLIFY (work_matrices%wfm_mo_virt_mo_occ(ispin)%matrix)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_virt(ispin), ncol_global=nmo_occ(ispin), &
                                  context=blacs_env)
         CALL cp_fm_create(work_matrices%wfm_mo_virt_mo_occ(ispin)%matrix, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
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

         DO istate = ngroups-igroup, nstates, ngroups
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_release_work_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, istate
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CALL timeset(routineN, handle)

      ! HFX-ralated matrices
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
         CALL cp_fm_release(work_matrices%wfm_mo_virt_mo_occ(ispin)%matrix)
         CALL cp_fm_release(work_matrices%wfm_mo_occ_mo_occ(ispin)%matrix)
         CALL cp_fm_release(work_matrices%wfm_ao_mo_occ(ispin)%matrix)
         CALL cp_fm_release(work_matrices%S_C0(ispin)%matrix)
         CALL cp_fm_release(work_matrices%S_C0_C0T(ispin)%matrix)
      END DO
      DEALLOCATE (work_matrices%S_C0, work_matrices%S_C0_C0T, work_matrices%wfm_ao_mo_occ)
      DEALLOCATE (work_matrices%wfm_mo_occ_mo_occ, work_matrices%wfm_mo_virt_mo_occ)

      DO ispin = SIZE(work_matrices%fm_pool_ao_mo_occ), 1, -1
         CALL fm_pool_release(work_matrices%fm_pool_ao_mo_occ(ispin)%pool)
      END DO
      DEALLOCATE (work_matrices%fm_pool_ao_mo_occ)

      ! sTDA
      IF (ASSOCIATED(work_matrices%gamma_exchange)) THEN
         CALL dbcsr_deallocate_matrix_set(work_matrices%gamma_exchange)
         NULLIFY (work_matrices%gamma_exchange)
      ENDIF
      IF (ASSOCIATED(work_matrices%ctransformed)) THEN
         CALL cp_fm_vect_dealloc(work_matrices%ctransformed)
         NULLIFY (work_matrices%ctransformed)
      ENDIF
      CALL dbcsr_release_p(work_matrices%shalf)
      ! Ewald
      CALL ewald_env_release(work_matrices%ewald_env)
      CALL ewald_pw_release(work_matrices%ewald_pw)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_release_work_matrices

! **************************************************************************************************
!> \brief Create kernel environment.
!> \param kernel_env       kernel environment (allocated and initialised on exit)
!> \param rho_struct_sub   ground state charge density
!> \param xc_section       input section which defines an exchange-correlation functional
!> \param is_rks_triplets  indicates that the triplet excited states calculation using
!>                         spin-unpolarised molecular orbitals has been requested
!> \param sub_env          parallel group environment
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
!>    * 06.2018 the charge density needs to be provided via a dummy argument [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_create_kernel_env(kernel_env, rho_struct_sub, xc_section, is_rks_triplets, sub_env)
      TYPE(full_kernel_env_type), INTENT(out)            :: kernel_env
      TYPE(qs_rho_type), POINTER                         :: rho_struct_sub
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(in)                                :: is_rks_triplets
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_create_kernel_env', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nao, nspins
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ
      LOGICAL                                            :: lsd
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_ij_r, rho_ij_r2, tau_ij_r, tau_ij_r2
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(section_vals_type), POINTER                   :: xc_fun_section

      CALL timeset(routineN, handle)

      nspins = SIZE(sub_env%mos_occ)
      lsd = (nspins > 1) .OR. is_rks_triplets

      DO ispin = 1, nspins
         CALL cp_fm_get_info(sub_env%mos_occ(ispin)%matrix, nrow_global=nao, ncol_global=nmo_occ(ispin))
      END DO

      CALL pw_env_get(sub_env%pw_env, auxbas_pw_pool=auxbas_pw_pool)

      CALL qs_rho_get(rho_struct_sub, rho_r=rho_ij_r, tau_r=tau_ij_r)

      NULLIFY (kernel_env%xc_rho_set, kernel_env%xc_rho1_set, kernel_env%xc_deriv_set)

      IF (is_rks_triplets) THEN
         ! we are about to compute triplet states using spin-restricted reference MOs;
         ! we still need the beta-spin density component in order to compute the TDDFT kernel
         ALLOCATE (rho_ij_r2(2))
         rho_ij_r2(1)%pw => rho_ij_r(1)%pw
         rho_ij_r2(2)%pw => rho_ij_r(1)%pw

         IF (ASSOCIATED(tau_ij_r)) THEN
            ALLOCATE (tau_ij_r2(2))
            tau_ij_r2(1)%pw => tau_ij_r(1)%pw
            tau_ij_r2(2)%pw => tau_ij_r(1)%pw
         END IF

         CALL xc_prep_2nd_deriv(kernel_env%xc_deriv_set, kernel_env%xc_rho_set, rho_ij_r2, &
                                auxbas_pw_pool, xc_section=xc_section, tau_r=tau_ij_r2)

         IF (ASSOCIATED(tau_ij_r)) &
            DEALLOCATE (tau_ij_r2)

         DEALLOCATE (rho_ij_r2)
      ELSE
         CALL xc_prep_2nd_deriv(kernel_env%xc_deriv_set, kernel_env%xc_rho_set, rho_ij_r, &
                                auxbas_pw_pool, xc_section=xc_section, tau_r=tau_ij_r)
      END IF

      ! ++ allocate structure for response density
      kernel_env%xc_section => xc_section
      kernel_env%deriv_method_id = section_get_ival(xc_section, "XC_GRID%XC_DERIV")
      kernel_env%rho_smooth_id = section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO")

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      kernel_env%xc_rho1_cflags = xc_functionals_get_needs(functionals=xc_fun_section, lsd=lsd, add_basic_components=.TRUE.)

      CALL xc_rho_set_create(kernel_env%xc_rho1_set, auxbas_pw_pool%pw_grid%bounds_local, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      kernel_env%alpha = 1.0_dp
      kernel_env%beta = 0.0_dp

      ! kernel_env%beta is taken into account in spin-restricted case only
      IF (nspins == 1) THEN
         IF (is_rks_triplets) THEN
            ! K_{triplets} = K_{alpha,alpha} - K_{alpha,beta}
            kernel_env%beta = -1.0_dp
         ELSE
            !                                                 alpha                 beta
            ! K_{singlets} = K_{alpha,alpha} + K_{alpha,beta} = 2 * K_{alpha,alpha} + 0 * K_{alpha,beta},
            ! due to the following relation : K_{alpha,alpha,singlets} == K_{alpha,beta,singlets}
            kernel_env%alpha = 2.0_dp
         END IF
      END IF

      CALL timestop(handle)
   END SUBROUTINE tddfpt_create_kernel_env

! **************************************************************************************************
!> \brief Release kernel environment.
!> \param kernel_env  kernel environment (destroyed on exit)
!> \par History
!>    * 02.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_release_kernel_env(kernel_env)
      TYPE(full_kernel_env_type), INTENT(inout)          :: kernel_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_release_kernel_env', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL xc_rho_set_release(kernel_env%xc_rho1_set)
      CALL xc_dset_release(kernel_env%xc_deriv_set)
      CALL xc_rho_set_release(kernel_env%xc_rho_set)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_release_kernel_env

END MODULE qs_tddfpt2_types
