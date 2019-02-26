!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief All the kernel specific subroutines for XAS TDP calculations
!> \author A. Bussy (03.2019)
! **************************************************************************************************

MODULE xas_tdp_kernel
   USE basis_set_types,                 ONLY: gto_basis_set_type, get_gto_basis_set
   USE cp_array_utils,                  ONLY: cp_2d_r_p_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply, dbcsr_deallocate_matrix_set,&
                                              dbcsr_allocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create, cp_fm_struct_release, &
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create, cp_fm_release, cp_fm_type, &
                                              cp_fm_get_submatrix
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_type, dbcsr_p_type, dbcsr_distribution_type, &
                                              dbcsr_get_info, dbcsr_create, dbcsr_get_block_p, &
                                              dbcsr_print, dbcsr_put_block, dbcsr_release, &
                                              dbcsr_finalize, dbcsr_copy, dbcsr_set, &
                                              dbcsr_multiply, dbcsr_desymmetrize, &
                                              dbcsr_complete_redistribute, dbcsr_iterator_type, &
                                              dbcsr_iterator_blocks_left, dbcsr_iterator_start, &
                                              dbcsr_iterator_stop, dbcsr_iterator_next_block, &
                                              dbcsr_add, dbcsr_reserve_block2d, &
                                              dbcsr_transposed, dbcsr_get_stored_coordinates, &
                                              dbcsr_type_symmetric
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_bcast, mp_isend, mp_irecv, mp_waitall
   USE qs_environment_types,            ONLY: qs_environment_type, get_qs_env 
   USE qs_kind_types,                   ONLY: get_qs_kind, qs_kind_type
   USE qs_o3c_methods,                  ONLY: contract3_o3c
   USE qs_o3c_types,                    ONLY: get_o3c_container, o3c_vec_type, get_o3c_vec , &
                                              o3c_iterate, o3c_iterator_create, &
                                              o3c_iterator_release, o3c_iterator_type, &
                                              get_o3c_iterator_info, o3c_container_type, &
                                              o3c_vec_create, o3c_vec_release
   USE util,                            ONLY: locate, get_limit
   USE xas_tdp_types,                   ONLY: donor_state_type, xas_tdp_env_type, &
                                              xas_tdp_control_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "xas_tdp_kernel"

   PUBLIC :: kernel_coulomb_xc, kernel_exchange

CONTAINS

! **************************************************************************************************
!> \brief Computes, if asked for it, the Coulomb and XC kernel matrices, in the usuall matrix format
!> \param coul_ker pointer the the Coulomb kernel matrix (can be void pointer)
!> \param xc_ker array of pointer to the different xc kernels (5 of them):
!>               1) the restricted closed-shell singlet kernel
!>               2) the restricted closed-shell triplet kernel
!>               3) the spin-conserving open-shell xc kernel
!>               4) the on-diagonal spin-flip open-shell xc kernel
!>               5) the off-diagonal spin-flip open-shell xc kernel
!> \param dist the dbcsr_dist of the kernel matrices
!> \param blk_size the block size of the kernel matrices
!> \param donor_state ...
!> \param xas_tdp_control ...
!> \param xas_tdp_env ...
!> \param qs_env ...
!> \note Coulomb and xc kernel are put together in the same routine because they use the same RI
!>       Coulomb: (aI|Jb) = (aI|P) (P|Q)^-1 (Q|Jb)
!>       XC : (aI|fxc|Jb) = (aI|P) (P|Q)^-1 (Q|fxc|R) (R|S)^-1 (S|Jb)
!>       In the above formula, a,b label the sgfs
!>       The routine analyses the xas_tdp_control to know which kernel must be computed and how
!>       (open-shell, singlet, triplet, ROKS, LSD, etc...)
!>       On entry, the pointers should be allocated
! **************************************************************************************************
   SUBROUTINE kernel_coulomb_xc(coul_ker, xc_ker, dist, blk_size, donor_state, xas_tdp_env, &
                                xas_tdp_control, qs_env)

      TYPE(dbcsr_type), INTENT(INOUT)                 :: coul_ker
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: xc_ker
      TYPE(dbcsr_distribution_type), POINTER          :: dist
      INTEGER, DIMENSION(:), POINTER                  :: blk_size
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "kernel_coulomb_xc", routineP = moduleN//":"//routineN

      INTEGER                                         :: handle, ndo_mo, ri_atom, iex, &
                                                         natom, nsgfp, source, ip, nex_atom, bo(2),&
                                                         ndo_so, i, ub, lb
      INTEGER, DIMENSION(:), POINTER                  :: ri_blk_size
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: lhs_int, rhs_int, matrices
      REAL(dp), DIMENSION(:,:), POINTER               :: PQ
      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist
      TYPE(cp_para_env_type), POINTER                 :: para_env
      LOGICAL                                         :: found, do_coulomb, do_xc, do_sg, do_tp, &
                                                         do_sc, do_sf

      NULLIFY(lhs_int, rhs_int, PQ, ri_blk_size, dbcsr_dist, matrices, para_env)

      CALL timeset(routineN, handle)

!  Initialization
      ndo_mo = donor_state%ndo_mo
      do_xc = xas_tdp_control%do_xc
      do_sg = xas_tdp_control%do_singlet
      do_tp = xas_tdp_control%do_triplet
      do_sc = xas_tdp_control%do_spin_cons
      do_sf = xas_tdp_control%do_spin_flip
      ndo_so = ndo_mo; IF (xas_tdp_control%do_uks) ndo_so = 2*ndo_mo
      ri_atom = donor_state%at_index
      CALL get_qs_env(qs_env, natom=natom, dbcsr_dist=dbcsr_dist, para_env=para_env)
      !we always need the coulomb kernel, except for triplets only or spin-flip only
      do_coulomb = .TRUE.
      IF ((do_tp .AND. .NOT. do_sg) .OR. (do_sf .AND. .NOT. do_sc)) do_coulomb = .FALSE.

! Pre-kernel calculation RI preparation

   !  Get the RI inverse coulomb for the excited kind (P|Q)^-1
      PQ => xas_tdp_env%ri_inv_coul
      nsgfp = SIZE(PQ, 1)

   !  Get the contracted 3-center integrals in lhs_int and copy them into rhs_int
      CALL contract_o3c_int(lhs_int, "COULOMB", donor_state, xas_tdp_env, xas_tdp_control, qs_env)
      ALLOCATE(rhs_int(ndo_so))
      CALL copy_ri_contr_int(rhs_int, lhs_int)

   !  In both Coulomb and XC need to multiply lhs_int by (P|Q)^-1
   !  Coming out of contract_o3c_int, the intrgrals are (aI_c|P), so still need to sum over c
   !  => create a symmetric matrix with the same block dist as contr_int and (P|Q)^-1 everywhere
      CALL dbcsr_allocate_matrix_set(matrices, 1)
      ALLOCATE(ri_blk_size(natom), matrices(1)%matrix)
      ri_blk_size = nsgfp
      CALL dbcsr_create(matrices(1)%matrix, name="WORK", matrix_type="S", dist=dbcsr_dist, &
                        row_blk_size=ri_blk_size, col_blk_size=ri_blk_size)
      CALL reserve_contraction_blocks(matrices, xas_tdp_env%lb_o3c_coul, xas_tdp_env%ub_o3c_coul, &
                                      qs_env, block_content=PQ)

      ! After this each row of LHS with index a contains (aI|P)*(P|Q)^-1 (identitcal block columns)
      DO i = 1, ndo_so
         CALL dbcsr_multiply('N', 'N', 1.0_dp, lhs_int(i)%matrix, matrices(1)%matrix, 0.0_dp, &
                             lhs_int(i)%matrix)
      END DO

      CALL dbcsr_deallocate_matrix_set(matrices)
      DEALLOCATE(ri_blk_size)

!  Deal with the Coulomb case  
      IF (do_coulomb) CALL coulomb(coul_ker, lhs_int, rhs_int, dist, blk_size, &
                                   xas_tdp_control, qs_env)

!  Deal with the XC case
      IF (do_xc) THEN

         ! In the end, we compute: (aI|fxc|Jb) = (aI|P) (P|Q)^-1 (Q|fxc|R) (R|S)^-1 (S|Jb)
         ! where fxc can take different spin contributions => Keep lhs_int as (aI|P) (P|Q)^-1

         ! Precompute the product (R|S)^-1 (S|Jb) and store it in rhs_int
         ! Careful: rhs_int is still in the (aI_c|P) format, the contraction over I is made in the
         !          final MM product
         CALL ri_all_blocks_mm(rhs_int, PQ)

      ! Find on which processor the integrals (Q|fxc|R) for this atom are stored
         CALL get_qs_env(qs_env, para_env=para_env)
         found = .FALSE.
         nex_atom = SIZE(xas_tdp_env%ex_atom_indices)
         DO ip = 0, para_env%num_pe-1

            bo = get_limit(nex_atom, para_env%num_pe, ip)
            DO iex = bo(1), bo(2)

               IF (xas_tdp_env%ex_atom_indices(iex) == ri_atom) THEN
                  source = ip
                  found = .TRUE.
                  EXIT
               END IF
            END DO !iex
            IF (found) EXIT
         END DO !ip

      ! Broadcast the integrals to all procs (delete it later)
         lb = 1; IF (do_sf .AND. .NOT. do_sc) lb = 4
         ub = 2; IF (do_sc) ub = 3; IF (do_sf) ub = 4
         DO i = lb,ub
            IF (.NOT. ASSOCIATED(xas_tdp_env%ri_fxc(ri_atom,i)%array)) THEN
               ALLOCATE(xas_tdp_env%ri_fxc(ri_atom,i)%array(nsgfp, nsgfp))
            END IF
            CALL mp_bcast(xas_tdp_env%ri_fxc(ri_atom,i)%array, source, para_env%group)
         END DO

      ! Case study on the calculation type
         IF (do_sg .OR. do_tp) THEN
            CALL rcs_xc(xc_ker(1)%matrix, xc_ker(2)%matrix, lhs_int, rhs_int, dist, blk_size, &
                        donor_state, xas_tdp_env, xas_tdp_control, qs_env)
         END IF

         IF (do_sc) THEN
            CALL sc_os_xc(xc_ker(3)%matrix, lhs_int, rhs_int, dist, blk_size, donor_state, &
                          xas_tdp_env, xas_tdp_control, qs_env)
         END IF

         IF (do_sf) THEN
            CALL ondiag_sf_os_xc(xc_ker(4)%matrix, lhs_int, rhs_int, dist, blk_size, donor_state, &
                                 xas_tdp_env, xas_tdp_control, qs_env)
         END IF

         IF (do_sf .AND. .NOT. xas_tdp_control%tamm_dancoff) THEN
            CPABORT("NYI")
         END IF
 
      END IF ! do_xc

!  Clean-up
      CALL dbcsr_deallocate_matrix_set(lhs_int)
      CALL dbcsr_deallocate_matrix_set(rhs_int)

      CALL timestop(handle)
   
   END SUBROUTINE kernel_coulomb_xc
   
! **************************************************************************************************
!> \brief Create the matrix containing the Coulomb kernel, which is:
!>        (aI_sigma|J_tau b) ~= (aI_sigma|P) * (P|Q) * (Q|J_tau b)
!> \param coul_ker the Coulomb kernel
!> \param lhs_int the LHS part of the RI: (aI_sigma|P)*(P|Q)^-1, each element of the array is the
!>                contraction for a specific donor spin-orbital
!> \param rhs_int the RHS part of the RI: (bJ_tau|Q), not fully contracted yet, to be transposed too
!> \param dist the inherited dbcsr ditribution
!> \param blk_size the inherited block sizes
!> \param xas_tdp_control ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE  coulomb(coul_ker, lhs_int, rhs_int, dist, blk_size, xas_tdp_control, qs_env)
      
      TYPE(dbcsr_type), INTENT(INOUT)                 :: coul_ker
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: lhs_int, rhs_int
      TYPE(dbcsr_distribution_type), POINTER          :: dist
      INTEGER, DIMENSION(:), POINTER                  :: blk_size
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "coulomb", routineP = moduleN//":"//routineN

      TYPE(dbcsr_type)                                :: work_mat
      LOGICAL                                         :: quadrants(3)

      ! Create a normal type work matrix
      CALL dbcsr_create(work_mat, name="WORK", matrix_type="N", dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)

      !Compute the product (rhs_int is still in (aI_b|P) format, that is not fully contracted over b
      !the sum over b is done during the MM
      !In the special case of ROKS, same MOs for each spin => put same (aI|Jb) product on the
      !alpha-alpha, alpha-beta and beta-beta quadrants of the kernel matrix.
      quadrants = .FALSE.
      IF (xas_tdp_control%do_roks) THEN
         quadrants(1:3) = .TRUE.
      ELSE
         quadrants(1) = .TRUE.
      END IF
      CALL ri_int_product(work_mat, lhs_int, rhs_int, quadrants, qs_env)
      CALL dbcsr_finalize(work_mat)

      !Create the symmetric kernel matrix and redistribute work_mat into it
      CALL dbcsr_create(coul_ker, name="COULOMB KERNEL", matrix_type="S", dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)
      CALL dbcsr_complete_redistribute(work_mat, coul_ker)

      !clean-up
      CALL dbcsr_release(work_mat)

   END SUBROUTINE coulomb

! **************************************************************************************************
!> \brief Create the matrix containing the XC kenrel in the spin-conserving open-shell case:
!>        (aI_sigma|fxc|J_tau b) ~= (aI_sigma|P) (P|Q)^-1 (Q|fxc|R) (R|S)^-1 (S|J_tau b)
!> \param xc_ker the kernel matrix
!> \param lhs_int the LHS part of the RI: (aI_sigma|P) (P|Q)^-1
!> \param rhs_int the RHS part of the RI: (bJ_tau|S) (S|R)^-1, not fully contracted yet
!> \param dist inherited dbcsr dist
!> \param blk_size inherited block sizes
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> note Prior to calling this function, the (Q|fxc|R) integral must be brodcasted to all procs
! **************************************************************************************************
   SUBROUTINE sc_os_xc(xc_ker, lhs_int, rhs_int, dist, blk_size, donor_state, xas_tdp_env, &
                       xas_tdp_control, qs_env)

      TYPE(dbcsr_type), INTENT(INOUT)                 :: xc_ker
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: lhs_int, rhs_int
      TYPE(dbcsr_distribution_type), POINTER          :: dist
      INTEGER, DIMENSION(:), POINTER                  :: blk_size
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "sc_os_xc", routineP = moduleN//":"//routineN

      INTEGER                                         :: ndo_mo, ri_atom
      TYPE(dbcsr_type)                                :: work_mat
      LOGICAL                                         :: quadrants(3)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: lhs_work

      NULLIFY(lhs_work)

   ! Initialization
      ndo_mo = donor_state%ndo_mo
      ri_atom = donor_state%at_index
      !normal type work matrix such that distribution of all spin quadrants match
      CALL dbcsr_create(work_mat, name="WORK", matrix_type="N", dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)
      !Create a lhs_work, in which the whole (aI_sigma|P) (P|Q)^-1 (Q|fxc|R) will be put
      ALLOCATE(lhs_work(SIZE(lhs_int)))

   ! Case study: UKS or ROKS ?
      IF (xas_tdp_control%do_uks) THEN

         ! In the case of UKS, donor MOs might be different for different spins. Moreover, the 
         ! fxc itself might change since fxc = fxc_sigma,tau
         ! => Carfully treat each spin-quadrant separately

      ! alpha-alpha spin quadrant (upper-lefet)
         quadrants = .FALSE.; quadrants(1) = .TRUE.

         ! Copy the alpha part of lhs_int, multiply by the alpha-alpha (Q|fxc|R) and then
         ! by the alpha part of rhs_int
         CALL copy_ri_contr_int(lhs_work(1:ndo_mo), lhs_int(1:ndo_mo))
         CALL ri_all_blocks_mm(lhs_work(1:ndo_mo), xas_tdp_env%ri_fxc(ri_atom,1)%array)
         CALL ri_int_product(work_mat, lhs_work(1:ndo_mo), rhs_int(1:ndo_mo), quadrants, qs_env)

      ! alpha-beta spin quadrant (upper-right)
         quadrants = .FALSE.; quadrants(2) = .TRUE.

         !Copy the alpha part of LHS, multiply by the alpha-beta kernel and the beta part of RHS
         CALL copy_ri_contr_int(lhs_work(1:ndo_mo), lhs_int(1:ndo_mo))
         CALL ri_all_blocks_mm(lhs_work(1:ndo_mo), xas_tdp_env%ri_fxc(ri_atom,2)%array)
         CALL ri_int_product(work_mat, lhs_work(1:ndo_mo), rhs_int(ndo_mo+1:2*ndo_mo), &
                             quadrants, qs_env)

      ! beta-beta spin quadrant (lower-right)
         quadrants = .FALSE.; quadrants(3) = .TRUE.

         !Copy the beta part of LHS, multiply by the beta-beta kernel and the beta part of RHS
         CALL copy_ri_contr_int(lhs_work(ndo_mo+1:2*ndo_mo), lhs_int(ndo_mo+1:2*ndo_mo))
         CALL ri_all_blocks_mm(lhs_work(ndo_mo+1:2*ndo_mo), xas_tdp_env%ri_fxc(ri_atom,3)%array)
         CALL ri_int_product(work_mat, lhs_work(ndo_mo+1:2*ndo_mo), rhs_int(ndo_mo+1:2*ndo_mo), &
                             quadrants, qs_env)

      ELSE IF (xas_tdp_control%do_roks) THEN

         ! In the case of ROKS, fxc = fxc_sigma,tau is different for each spin quadrant, but the
         ! donor MOs remain the same

      ! alpha-alpha kernel in the upper left quadrant
         quadrants = .FALSE.; quadrants(1) = .TRUE.

         !Copy the LHS and multiply by alpha-alpha kernel
         CALL copy_ri_contr_int(lhs_work, lhs_int)
         CALL ri_all_blocks_mm(lhs_work, xas_tdp_env%ri_fxc(ri_atom,1)%array)
         CALL ri_int_product(work_mat, lhs_work, rhs_int, quadrants, qs_env)

      ! alpha-beta kernel in the upper-right quadrant
         quadrants = .FALSE.; quadrants(2) = .TRUE.

         !Copy LHS and multiply by the alpha-beta kernel
         CALL copy_ri_contr_int(lhs_work, lhs_int)
         CALL ri_all_blocks_mm(lhs_work, xas_tdp_env%ri_fxc(ri_atom,2)%array)
         CALL ri_int_product(work_mat, lhs_work, rhs_int, quadrants, qs_env)

      ! beta-beta kernel in the lower-right quadrant
         quadrants = .FALSE.; quadrants(3) = .TRUE.

         !Copy the LHS and multiply by the beta-beta kernel
         CALL copy_ri_contr_int(lhs_work, lhs_int)
         CALL ri_all_blocks_mm(lhs_work, xas_tdp_env%ri_fxc(ri_atom,3)%array)
         CALL ri_int_product(work_mat, lhs_work, rhs_int, quadrants, qs_env)

      END IF
      CALL dbcsr_finalize(work_mat)

   ! Create a symmetric kernel matrix and redistribute the normal work matrix into it
      CALL dbcsr_create(xc_ker, name="SC OS XC KERNEL", matrix_type="S", dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)
      CALL dbcsr_complete_redistribute(work_mat, xc_ker)

   !clean-up
      CALL dbcsr_deallocate_matrix_set(lhs_work)
      CALL dbcsr_release(work_mat)

   END SUBROUTINE sc_os_xc

! **************************************************************************************************
!> \brief Create the matrix containing the on-diagonal spin-flip XC kernel (open-shell), which is:
!>        (a I_sigma|fxc|J_tau b) * delta_sigma,tau, fxc = 1/(rhoa-rhob) * (dE/drhoa - dE/drhob)
!>        with RI: (a I_sigma|fxc|J_tau b) ~= (aI_sigma|P) (P|Q)^-1 (Q|fxc|R) (R|S)^-1 (S|J_tau b)
!> \param xc_ker the kernel matrix
!> \param lhs_int the LHS part of the RI: (aI_sigma|P) (P|Q)^-1
!> \param rhs_int the RHS part of the RI: (bJ_tau|S) (S|R)^-1, not fully contracted yet
!> \param dist inherited dbcsr dist
!> \param blk_size inherited block sizes
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note  It must be later on multiplied by the spin-swapped Q projector
!>        Prior to calling this function, the (Q|fxc|R) integral must be brodcasted to all procs
! **************************************************************************************************
   SUBROUTINE ondiag_sf_os_xc(xc_ker, lhs_int, rhs_int, dist, blk_size, donor_state, xas_tdp_env, &
                              xas_tdp_control, qs_env)

      TYPE(dbcsr_type), INTENT(INOUT)                 :: xc_ker
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: lhs_int, rhs_int
      TYPE(dbcsr_distribution_type), POINTER          :: dist
      INTEGER, DIMENSION(:), POINTER                  :: blk_size
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "ondiag_sf_os_xc", routineP = moduleN//":"//routineN

      INTEGER                                         :: ndo_mo, ri_atom
      TYPE(dbcsr_type)                                :: work_mat
      LOGICAL                                         :: quadrants(3)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: lhs_work

      NULLIFY(lhs_work)

   ! Initialization
      ndo_mo = donor_state%ndo_mo
      ri_atom = donor_state%at_index
      !normal type work matrix such that distribution of all spin quadrants match
      CALL dbcsr_create(work_mat, name="WORK", matrix_type="N", dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)
      !Create a lhs_work, in which the whole (aI_sigma|P) (P|Q)^-1 (Q|fxc|R) will be put
      !because in spin-flip fxc is spin-independant, can take the product once and for all
      ALLOCATE(lhs_work(SIZE(lhs_int)))
      CALL copy_ri_contr_int(lhs_work(:), lhs_int(:))
      CALL ri_all_blocks_mm(lhs_work(:), xas_tdp_env%ri_fxc(ri_atom,4)%array)

   ! Case study: UKS or ROKS ?
      IF (xas_tdp_control%do_uks) THEN

         ! In the case of UKS, donor MOs might be different for different spins
         ! => Carfully treat each spin-quadrant separately
         ! NO alpha-beta because of the delta_sigma,tau

      ! alpha-alpha spin quadrant (upper-lefet)
         quadrants = .FALSE.; quadrants(1) = .TRUE.
         CALL ri_int_product(work_mat, lhs_work(1:ndo_mo), rhs_int(1:ndo_mo), quadrants, qs_env)

      ! beta-beta spin quadrant (lower-right)
         quadrants = .FALSE.; quadrants(3) = .TRUE.
         CALL ri_int_product(work_mat, lhs_work(ndo_mo+1:2*ndo_mo), rhs_int(ndo_mo+1:2*ndo_mo), &
                             quadrants, qs_env)

      ELSE IF (xas_tdp_control%do_roks) THEN

         ! In the case of ROKS, same donor MOs for both spins => can do it all at once
         ! But NOT the alpha-beta quadrant because of delta_sigma,tau

         quadrants(1) = .TRUE.; quadrants(2) = .FALSE.; quadrants(3) = .TRUE.
         CALL ri_int_product(work_mat, lhs_work, rhs_int, quadrants, qs_env)

      END IF
      CALL dbcsr_finalize(work_mat)

   ! Create a symmetric kernel matrix and redistribute the normal work matrix into it
      CALL dbcsr_create(xc_ker, name="ON-DIAG SF OS XC KERNEL", matrix_type="S", dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)
      CALL dbcsr_complete_redistribute(work_mat, xc_ker)

   !clean-up
      CALL dbcsr_deallocate_matrix_set(lhs_work)
      CALL dbcsr_release(work_mat)

   END SUBROUTINE ondiag_sf_os_xc

! **************************************************************************************************
!> \brief Create the matrix containing the XC kernel in the restricted closed-shell case, for
!>        singlets: (aI|fxc|Jb) = (aI|P) (P|Q)^-1 (Q|fxc_alp,alp+fxc_alp,bet|R) (R|S)^-1 (S|Jb)
!>        triplets: (aI|fxc|Jb) = (aI|P) (P|Q)^-1 (Q|fxc_alp,alp-fxc_alp,bet|R) (R|S)^-1 (S|Jb)
!> \param sg_xc_ker the singlet kernel matrix
!> \param tp_xc_ker the triplet kernel matrix
!> \param lhs_int the LHS part of the RI: (aI|P) (P|Q)^-1
!> \param rhs_int the RHS part of the RI: (bJ|S) (S|R)^-1, not fully contracted yet
!> \param dist inherited dbcsr dist
!> \param blk_size inherited block sizes
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE rcs_xc(sg_xc_ker, tp_xc_ker, lhs_int, rhs_int, dist, blk_size, donor_state, &
                     xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(dbcsr_type), INTENT(INOUT)                 :: sg_xc_ker, tp_xc_ker
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: lhs_int, rhs_int
      TYPE(dbcsr_distribution_type), POINTER          :: dist
      INTEGER, DIMENSION(:), POINTER                  :: blk_size
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "rcx_xc", routineP = moduleN//":"//routineN

      INTEGER                                         :: nsgfp, ri_atom
      LOGICAL                                         :: quadrants(3)
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: fxc
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: lhs_work
      TYPE(dbcsr_type)                                :: work_mat

      NULLIFY(lhs_work)

   ! Initialization
      ri_atom = donor_state%at_index
      nsgfp = SIZE(xas_tdp_env%ri_fxc(ri_atom,1)%array, 1)
      !Create a lhs_work in which the wholw (aI|P) (P|Q)^-1 (Q|fxc|R) will be placed
      ALLOCATE(lhs_work(SIZE(lhs_int)))

   ! Case study: singlet and/or triplet ?
      IF (xas_tdp_control%do_singlet) THEN

         ! singlet work structures
         ALLOCATE(fxc(nsgfp, nsgfp))
         CALL dbcsr_create(work_mat, name="WORK", matrix_type="N", dist=dist, &
                           row_blk_size=blk_size, col_blk_size=blk_size)

         ! Take the sum of fxc for alpha-alpha and alpha-beta
         CALL dcopy(nsgfp*nsgfp, xas_tdp_env%ri_fxc(ri_atom,1)%array, 1, fxc, 1)
         CALL daxpy(nsgfp*nsgfp, 1.0_dp, xas_tdp_env%ri_fxc(ri_atom,2)%array, 1, fxc, 1)

         ! Copy the fresh lhs_int = (aI|P) (P|Q)^-1 and multiply by (Q|fxc|R)
         CALL copy_ri_contr_int(lhs_work, lhs_int)
         CALL ri_all_blocks_mm(lhs_work, fxc)

         ! Compute the final LHS RHS product => spin-restricted, only upper-left quadrant
         quadrants = .FALSE.; quadrants(1) = .TRUE.
         CALL ri_int_product(work_mat, lhs_work, rhs_int, quadrants, qs_env)
         CALL dbcsr_finalize(work_mat)

         !Create the symmetric kernel matrix and redistribute work_mat into it
         CALL dbcsr_create(sg_xc_ker, name="XC SINGLET KERNEL", matrix_type="S", dist=dist, &
                           row_blk_size=blk_size, col_blk_size=blk_size)
         CALL dbcsr_complete_redistribute(work_mat, sg_xc_ker)

         !singlet clean-up
         CALL dbcsr_release(work_mat)
         DEALLOCATE(fxc)

      END IF

      IF (xas_tdp_control%do_triplet) THEN

         ! triplet work structures
         ALLOCATE(fxc(nsgfp, nsgfp))
         CALL dbcsr_create(work_mat, name="WORK", matrix_type="N", dist=dist, &
                           row_blk_size=blk_size, col_blk_size=blk_size)

         ! Take the difference of fxc for alpha-alpha and alpha-beta
         CALL dcopy(nsgfp*nsgfp, xas_tdp_env%ri_fxc(ri_atom,1)%array, 1, fxc, 1)
         CALL daxpy(nsgfp*nsgfp, -1.0_dp, xas_tdp_env%ri_fxc(ri_atom,2)%array, 1, fxc, 1)

         ! Copy the fresh lhs_int = (aI|P) (P|Q)^-1 and multiply by (Q|fxc|R)
         CALL copy_ri_contr_int(lhs_work, lhs_int)
         CALL ri_all_blocks_mm(lhs_work, fxc)

         ! Compute the final LHS RHS product => spin-restricted, only upper-left quadrant
         quadrants = .FALSE.; quadrants(1) = .TRUE.
         CALL ri_int_product(work_mat, lhs_work, rhs_int, quadrants, qs_env)
         CALL dbcsr_finalize(work_mat)

         !Create the symmetric kernel matrix and redistribute work_mat into it
         CALL dbcsr_create(tp_xc_ker, name="XC TRIPLET KERNEL", matrix_type="S", dist=dist, &
                           row_blk_size=blk_size, col_blk_size=blk_size)
         CALL dbcsr_complete_redistribute(work_mat, tp_xc_ker)

         ! triplet clean-up
         CALL dbcsr_release(work_mat)
         DEALLOCATE(fxc)

      END IF

   ! clean-up
      CALL dbcsr_deallocate_matrix_set(lhs_work)

   END SUBROUTINE rcs_xc

! **************************************************************************************************
!> \brief Computes the exact exchange kernel matrix using RI. Returns an array of 3 matrices, 
!>        which are:
!>        1) the on-diagonal kernel: (ab|I_sigma J_tau) * delta_sigma,tau 
!>        2) the off-diagonal spin-conserving kernel: (aJ_sigma|I_tau b) * delta_sigma,tau
!>        3) the off-diagonal spin-flip kernel: (aJ_tau|I_sigma b) * (1 - delta_sigma,tau)
!>        An internal analysis determines which of the above are computed (can range from 0 to 3),
!>        Only to be called if some HFX is required in the first place
!> \param dist the usual dbcsr distribution
!> \param blk_size the block size
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note In the case of spin-conserving excitation, the kernel must later be multiplied by the 
!>       usual Q projector. In the case of spin-flip, one needs to project the excitations coming
!>       from alpha donor MOs on the unoccupied beta MOs. This is done by multiplying by a Q 
!>       projector where the alpha-alpha and beta-beta quadrants are swapped
!>       The ex_ker array should be allocated on entry (not the internals)
! **************************************************************************************************
   SUBROUTINE kernel_exchange(ex_ker, dist, blk_size, donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: ex_ker
      TYPE(dbcsr_distribution_type), POINTER          :: dist
      INTEGER, DIMENSION(:), POINTER                  :: blk_size
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env
      
      CHARACTER(len=*), PARAMETER :: routineN = "kernel_exchange", routineP = moduleN//":"//routineN

      INTEGER                                         :: handle, nsgfp, natom, i, ndo_so
      LOGICAL                                         :: do_off_sf, do_off_sc
      INTEGER, DIMENSION(:), POINTER                  :: ri_blk_size
      REAL(dp), DIMENSION(:,:), POINTER               :: PQ
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: lhs_int, rhs_int, matrices
      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist

      NULLIFY(PQ, lhs_int, rhs_int, dbcsr_dist, matrices)

      !Don't do anything if no hfx
      IF (.NOT. xas_tdp_control%do_hfx) RETURN

      CALL timeset(routineN, handle)
      
!  Initialization
      CALL get_qs_env(qs_env, natom=natom, dbcsr_dist=dbcsr_dist)
      ndo_so = donor_state%ndo_mo; IF (xas_tdp_control%do_uks) ndo_so = 2*ndo_so

      !compute the off-dig spin-conserving only if not TDA and anything that is spin-conserving
      do_off_sc = (.NOT. xas_tdp_control%tamm_dancoff) .AND. (xas_tdp_control%do_spin_cons &
                   .OR. xas_tdp_control%do_singlet .OR. xas_tdp_control%do_triplet)
      !compute the off-diag spin-flip if not TDA and spin-flip
      do_off_sf = (.NOT. xas_tdp_control%tamm_dancoff) .AND. xas_tdp_control%do_spin_flip

!  Preparing for on- and off-diagonal exchange kernels
      ! Get the RI inverse coulomb for the excited kind (P|Q)^-1
      PQ => xas_tdp_env%ri_inv_coul
      nsgfp = SIZE(PQ, 1)

      ! Need the once contracted integrals (aI_c|P)  (not fully contracted over c yet)
      CALL contract_o3c_int(lhs_int, "EXCHANGE", donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      ! Also need a copy if off-diag
      IF (do_off_sc .OR. do_off_sf) THEN
         ALLOCATE(rhs_int(ndo_so))
         CALL copy_ri_contr_int(rhs_int, lhs_int)
      END IF

      ! Multiply lhs by (P|Q)^-1, ensuring the full contraction over b too
      CALL dbcsr_allocate_matrix_set(matrices, 1)
      ALLOCATE(ri_blk_size(natom), matrices(1)%matrix)
      ri_blk_size = nsgfp
      CALL dbcsr_create(matrices(1)%matrix, name="WORK", matrix_type="S", dist=dbcsr_dist, &
                        row_blk_size=ri_blk_size, col_blk_size=ri_blk_size)
      CALL reserve_contraction_blocks(matrices, xas_tdp_env%lb_o3c_coul, xas_tdp_env%ub_o3c_coul, &
                                      qs_env, block_content=PQ)

      ! After this each row of LHS with index a contains (aI|P)*(P|Q)^-1 (identitcal block columns)
      DO i = 1, ndo_so
         CALL dbcsr_multiply('N', 'N', 1.0_dp, lhs_int(i)%matrix, matrices(1)%matrix, 0.0_dp, &
                             lhs_int(i)%matrix)
      END DO

      CALL dbcsr_deallocate_matrix_set(matrices)
      DEALLOCATE(ri_blk_size)

!  The on-diagonal exchange : (ab|P) * (P|Q)^-1 * (Q|I_sigma J_tau) * delta_sigma,tau
      CALL ondiag_ex(ex_ker(1)%matrix, lhs_int, dist, blk_size, donor_state, xas_tdp_env, &
                     xas_tdp_control, qs_env)

!  The off-diag spin-conserving case: (aJ_sigma|P) * (P|Q)^-1 * (Q|I_tau b) * delta_sigma,tau
      IF (do_off_sc) THEN
         CALL offdiag_ex_sc(ex_ker(2)%matrix, lhs_int, rhs_int, dist, blk_size, &
                            donor_state, xas_tdp_control, qs_env)
      END IF
         
!  The off-diag spin-flip case: (aJ_tau|P) * (P|Q)^-1 * (Q|I_sigma b) * (1 - delta_sigma,tau)
      IF (do_off_sf) THEN
         CPABORT("NYI")
      END IF

      !clean-up 
      CALL dbcsr_deallocate_matrix_set(lhs_int)
      IF (do_off_sc .OR. do_off_sf) CALL dbcsr_deallocate_matrix_set(rhs_int)

      CALL timestop(handle)

   END SUBROUTINE kernel_exchange

! **************************************************************************************************
!> \brief Create the matrix containing the on-diagonal exact exchange kernel, which is:
!>        (ab|I_sigma J_tau) * delta_sigma,tau, where a,b are AOs, I_sigma and J_tau are the donor
!>        spin-orbitals. A RI is done: (ab|I_sigma J_tau) = (ab|P) * (P|Q)^-1 * (Q|I_sigma J_tau)
!> \param ondiag_ex_ker the on-diagonal exchange kernel in dbcsr format
!> \param contr1_int the already once-contracted RI 3-center integrals (aI_sigma|P)*(P|Q)^-1, 
!>        where each matrix of the array contains the contraction for the donor spin-orbital I_sigma,
!>        which is already multiply by the inverse RI 2-electron Coulomb integral. For a given row
!>        with index a, all columns are identical and contain (aI_sigma|P)*(P|Q)^-1
!> \param dist the inherited dbcsr distribution
!> \param blk_size the inherited dbcsr block sizes
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE ondiag_ex(ondiag_ex_ker, contr1_int, dist, blk_size, donor_state, xas_tdp_env, &
                            xas_tdp_control, qs_env)

      TYPE(dbcsr_type), INTENT(INOUT)                 :: ondiag_ex_ker
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: contr1_int
      TYPE(dbcsr_distribution_type), POINTER          :: dist
      INTEGER, DIMENSION(:), POINTER                  :: blk_size
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "ondiag_ex", routineP = moduleN//":"//routineN

      INTEGER                                         :: s, ndo_so, ndo_mo, ri_atom, nblk, iso, &
                                                         jso, sr, er, iblk, jblk, blk
      LOGICAL                                         :: do_uks, do_roks, found
      TYPE(cp_para_env_type), POINTER                 :: para_env
      TYPE(cp_blacs_env_type), POINTER                :: blacs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      INTEGER, DIMENSION(:), POINTER                  :: ri_blk_size, vec_blk_size
      TYPE(cp_fm_struct_type), POINTER                :: fm_struct
      TYPE(cp_fm_type), POINTER                       :: contr2_int
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: coeffs
      TYPE(o3c_vec_type), DIMENSION(:), POINTER       :: o3c_vec
      TYPE(dbcsr_type), POINTER                       :: abIJ
      TYPE(dbcsr_type)                                :: contr1_int_t, abIJ_desymm, work_mat
      REAL(dp), DIMENSION(:), POINTER                 :: v
      REAL(dp), DIMENSION(:,:), POINTER               :: pblock
      TYPE(dbcsr_iterator_type)                       :: iter

      NULLIFY(para_env, blacs_env, matrix_s, ri_blk_size, fm_struct, contr2_int, o3c_vec, v)
      NULLIFY(vec_blk_size, pblock, abIJ)

      ! Initialization
      ndo_mo = donor_state%ndo_mo
      ri_atom = donor_state%at_index
      do_roks = xas_tdp_control%do_roks
      do_uks = xas_tdp_control%do_uks
      ndo_so = ndo_mo; IF (do_uks) ndo_so = 2*ndo_mo !if not UKS, same donor MOs for both spins

      ! We want to compute (ab|I_sigma J_tau) = (ab|P) * (P|Q)^-1 * (Q|I_sigma J_tau) 
      ! Already have (cJ_tau|P) * (P|Q)^-1 stored in contr1_int. Need to further contract the
      ! AOs c with the coeff of the I_alpha spin-orbital. Do that via matrix multiplication
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(contr1_int(1)%matrix, nfullcols_total=s, col_blk_size=ri_blk_size)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, para_env=para_env, &
                               nrow_global=s, ncol_global=ndo_so)
      CALL cp_fm_create(contr2_int, fm_struct)
      ALLOCATE(coeffs(ri_blk_size(ri_atom), ndo_so))

      nblk = SIZE(ri_blk_size)
      ALLOCATE(o3c_vec(nblk), vec_blk_size(nblk))
      vec_blk_size = 1; vec_blk_size(ri_atom) = ri_blk_size(ri_atom)
      CALL o3c_vec_create(o3c_vec, vec_blk_size)

      ! a and b need to overlap for non-zero (ab|IJ) => same block structure as overlap S
      ALLOCATE(abIJ)
      CALL dbcsr_copy(abIJ, matrix_s(1)%matrix, name="(ab|IJ)")
      CALL dbcsr_create(work_mat, name="WORK", matrix_type="N", dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)

      ! Loop over donor spin-orbitals. End matrix is symmetric => span only upper half
      DO iso = 1, ndo_so

         ! Contract the second time with coeffs of spin-orbital iso
         CALL dbcsr_transposed(contr1_int_t, contr1_int(iso)%matrix)
         CALL cp_dbcsr_sm_fm_multiply(contr1_int_t, donor_state%gs_coeffs, contr2_int, ncol=ndo_so)
         sr = 1
         IF (ri_atom > 1) sr = SUM(ri_blk_size(1:ri_atom-1))+1
         er = SUM(ri_blk_size(1:ri_atom))

         ! take (P|Q)^-1 * (Q|IJ) and put that in coeffs 
         CALL cp_fm_get_submatrix(contr2_int, coeffs, start_row=sr, start_col=1, n_rows=er-sr+1, &
                                  n_cols=ndo_so)

         DO jso = iso, ndo_so

            ! There is no alpha-beta exchange. In case of UKS, iso,jso span all spin-orbitals
            ! => CYCLE if iso and jso are indexing MOs with different spin (and we have UKS)
            IF (do_uks .AND. (iso <= ndo_mo .AND. jso > ndo_mo)) CYCLE

            ! compute (ab|IJ) = sum_P (ab|P) * (P|Q)^-1 * (Q|IJ)
            CALL get_o3c_vec(o3c_vec, ri_atom, v)
            v(:) = coeffs(:,jso)
            CALL dbcsr_set(abIJ, 0.0_dp)
            CALL contract3_o3c(xas_tdp_env%ri_o3c_ex, o3c_vec, abIJ)

            ! (ab|IJ) is symmetric, but need it as normal for distributions to match
            CALL dbcsr_desymmetrize(abIJ, abIJ_desymm)

            CALL dbcsr_iterator_start(iter, abIJ_desymm)
            DO WHILE (dbcsr_iterator_blocks_left(iter))

               CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
               IF (iso == jso .AND. jblk<iblk) CYCLE

               CALL dbcsr_get_block_p(abIJ_desymm, iblk, jblk, pblock, found)

               IF (found) THEN
                  CALL dbcsr_put_block(work_mat, (iso-1)*nblk+iblk, (jso-1)*nblk+jblk, pblock)

                  !In case of ROKS, we have (ab|IJ) for alpha-alpha spin, but it is the same for
                  !beta-beta => replicate the blocks (alpha-beta is zero)
                  IF (do_roks) THEN
                     !the beta-beta block
                     CALL dbcsr_put_block(work_mat, (ndo_so+iso-1)*nblk+iblk, (ndo_so+jso-1)*nblk+jblk, pblock)
                  END IF
               END IF

            END DO !iterator
            CALL dbcsr_iterator_stop(iter)
            CALL dbcsr_release(abIJ_desymm)

         END DO !jso
         CALL dbcsr_release(contr1_int_t)
      END DO !iso

      CALL dbcsr_finalize(work_mat)
      CALL dbcsr_create(ondiag_ex_ker, name="ONDIAG EX KERNEL", matrix_type="S", dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)
      CALL dbcsr_complete_redistribute(work_mat, ondiag_ex_ker)

      !Clean-up
      CALL cp_fm_release(contr2_int)
      CALL cp_fm_struct_release(fm_struct)
      CALL dbcsr_release(work_mat)
      CALL dbcsr_release(abIJ)
      CALL o3c_vec_release(o3c_vec)
      DEALLOCATE(o3c_vec, vec_blk_size, abIJ)

   END SUBROUTINE ondiag_ex

! **************************************************************************************************
!> \brief Create the matrix containing the off-diagonal exact exchange kernel in the spin-conserving
!>        case (which also includes excitations from the closed=shell ref state ) This matrix reads:
!>        (aJ_sigma|I_tau b) * delta_sigma,tau , where a, b are AOs and J_sigma, I_tau are the donor
!>        spin-orbital. A RI is done: (aJ_sigma|I_tau b) = (aJ_sigma|P) * (P|Q)^-1 * (Q|I_tau b)
!> \param offdiag_ex_ker the off-diagonal, spin-conserving exchange kernel in dbcsr format
!> \param lhs_int the LHS part of the RI: (aJ_sigma|P)*(P|Q)^-1, each element of the array is the
!>                contraction for a specific donor spin-orbital
!> \param rhs_int the RHS part of the RI: (b I_tau|Q), still needs transposing
!> \param dist the inherited dbcsr ditribution
!> \param blk_size the inherited block sizes 
!> \param donor_state ...
!> \param xas_tdp_control ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE offdiag_ex_sc(offdiag_ex_ker, lhs_int, rhs_int, dist, blk_size, donor_state, &
                            xas_tdp_control, qs_env)

      TYPE(dbcsr_type), INTENT(INOUT)                 :: offdiag_ex_ker
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: lhs_int
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: rhs_int
      TYPE(dbcsr_distribution_type), POINTER          :: dist
      INTEGER, DIMENSION(:), POINTER                  :: blk_size
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "offdiag_ex_sc", routineP = moduleN//":"//routineN

      INTEGER                                         :: ndo_mo
      TYPE(dbcsr_type)                                :: work_mat
      LOGICAL                                         :: quadrants(3), do_roks, do_uks

      !Initialization
      ndo_mo = donor_state%ndo_mo
      do_roks = xas_tdp_control%do_roks
      do_uks = xas_tdp_control%do_uks

      !Given the lhs_int and rhs_int, all we need to do is multiply elements from the former by
      !the transpose of the later, and put the result in the correct spin quadrants

      !Create a normal type work matrix
      CALL dbcsr_create(work_mat, name="WORK", matrix_type="N", dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)

      !Case study on closed-shell, ROKS or UKS
      quadrants = .FALSE.
      IF (do_roks) THEN
         !In ROKS, the donor MOs for each spin are the same => copy the product in both the
         !alpha-alpha and the beta-beta quadrants
         quadrants(1) = .TRUE.
         quadrants(3) = .TRUE.
         CALL ri_int_product(work_mat, lhs_int, rhs_int, quadrants, qs_env, mo_transpose=.TRUE.)
      ELSE IF (do_uks) THEN
         !In UKS, the donor MOs are possibly different for each spin => start with the 
         !alpha-alpha product and the perform the beta-beta product separately
         quadrants(1) = .TRUE.
         CALL ri_int_product(work_mat, lhs_int(1:ndo_mo), rhs_int(1:ndo_mo), quadrants, &
                             qs_env, mo_transpose=.TRUE.)
         quadrants(1) = .FALSE.; quadrants(3) = .TRUE.
         CALL ri_int_product(work_mat, lhs_int(ndo_mo+1:2*ndo_mo), rhs_int(ndo_mo+1:2*ndo_mo), &
                             quadrants, qs_env, mo_transpose=.TRUE.)
      ELSE
         !In the restricted closed-shell case, only have one spin and a single qudarant
         quadrants(1) = .TRUE.
         CALL ri_int_product(work_mat, lhs_int, rhs_int, quadrants, qs_env, mo_transpose=.TRUE.)
      END IF
      CALL dbcsr_finalize(work_mat)

      !Create the symmetric kernel matrix and redistribute work_mat into it
      CALL dbcsr_create(offdiag_ex_ker, name="OFFDIAG EX KERNEL", matrix_type="S", dist=dist, &
                        row_blk_size=blk_size, col_blk_size=blk_size)
      CALL dbcsr_complete_redistribute(work_mat, offdiag_ex_ker)

      !clean-up
      CALL dbcsr_release(work_mat)

   END SUBROUTINE offdiag_ex_sc

! **************************************************************************************************
!> \brief Contract the ri 3-center integrals contained in o3c with repect to the donor MOs coeffs, 
!>        for a given excited atom k => (aI|k) = sum_b c_Ib (ab|k)
!> \param contr_int the contracted integrals as array of dbcsr matrices
!> \param op_type for which operator type we contract (COULOMB or EXCHANGE)
!> \param donor_state ...
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note  In the output matrices, (aI_b|k) is stored at block a,b where I_b is the partial 
!>        contraction that only includes coeffs from atom b. Note that the contracted matrix is
!>        not symmetric, but the blocks (aI_b|P) and (bI_a|P) are strictly identical
! **************************************************************************************************
   SUBROUTINE contract_o3c_int(contr_int, op_type, donor_state, xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: contr_int
      CHARACTER(len=*), INTENT(IN)                    :: op_type
      TYPE(donor_state_type), POINTER                 :: donor_state
      TYPE(xas_tdp_env_type), POINTER                 :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER             :: xas_tdp_control
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "contract_o3c_int", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: handle, natom, ndo_mo, imo, nsgfk, kkind, &
                                                         ub, lb, sr, er, iat, katom, nspins, ndo_so,&
                                                         ispin, i, sr2, er2
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s, matrices
      TYPE(dbcsr_type), POINTER                       :: aI_bP, PI_ab, work1, work2 
      CHARACTER                                       :: my_op
      INTEGER, DIMENSION(:), POINTER                  :: std_blk_size, ri_blk_size
      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(gto_basis_set_type), POINTER               :: ri_basis
      TYPE(o3c_vec_type), DIMENSION(:), POINTER       :: o3c_vec
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: coeffs
      REAL(dp), DIMENSION(:), POINTER                 :: v
      TYPE(o3c_container_type), POINTER               :: o3c_coul, o3c_ex
      INTEGER, DIMENSION(:), ALLOCATABLE              :: nsgf_orb
      LOGICAL                                         :: do_uks
      TYPE(cp_para_env_type), POINTER                 :: para_env

      NULLIFY(matrix_s, std_blk_size, ri_blk_size, dbcsr_dist, qs_kind_set, ri_basis, o3c_vec)
      NULLIFY(aI_bP, PI_ab, work1, work2, matrices)

      CALL timeset(routineN, handle)

!  Initialization
      CALL get_qs_env(qs_env, natom=natom, matrix_s=matrix_s, qs_kind_set=qs_kind_set, &
                      dbcsr_dist=dbcsr_dist, para_env=para_env)
      ndo_mo = donor_state%ndo_mo
      kkind = donor_state%kind_index
      katom = donor_state%at_index
      !by default contract for Coulomb
      my_op = "C"
      o3c_coul => xas_tdp_env%ri_o3c_coul
      IF (op_type == "EXCHANGE") THEN
         my_op = "E"
         CPASSERT(ASSOCIATED(xas_tdp_env%ri_o3c_ex))
         o3c_ex => xas_tdp_env%ri_o3c_ex
      END IF
      do_uks = xas_tdp_control%do_uks
      nspins = 1; IF (do_uks) nspins = 2
      ndo_so = nspins*ndo_mo

!  contracted integrals block sizes
      CALL dbcsr_get_info(matrix_s(1)%matrix, col_blk_size=std_blk_size)
      ! getting the block dimensions for the RI basis
      CALL get_qs_kind(qs_kind_set(kkind), basis_set=ri_basis, basis_type="RI_XAS")
      CALL get_gto_basis_set(ri_basis, nsgf=nsgfk)
      ALLOCATE(ri_blk_size(natom))
      ri_blk_size = nsgfk

!  Get some info about the o3c_coul. Whether we contract for coulomb or exchange, MOs are assumed
!  to be localized as in the Coulomb integrals => get the lowest/highest a,b such that (ab|P) != 0
      lb = xas_tdp_env%lb_o3c_coul
      ub = xas_tdp_env%ub_o3c_coul

!  Create  work matrices. They must have the block distribution of a symmetric matric to enter the 
!  o3c routines. HACK => create a symmetric matrix with non-symmetric block size (aI_bP, PI_ab)
!  Also create their normal matrix conterparts (work1, work2)
      ALLOCATE(aI_bP, PI_ab, work1, work2, matrices(2))
      CALL dbcsr_create(aI_bP, dist=dbcsr_dist, matrix_type="S", name="(aI_b|P)", &
                        row_blk_size=std_blk_size, col_blk_size=ri_blk_size)
      CALL dbcsr_create(work1, dist=dbcsr_dist, matrix_type="N", name="WORK 1", &
                        row_blk_size=std_blk_size, col_blk_size=ri_blk_size)

      CALL dbcsr_create(PI_ab, dist=dbcsr_dist, matrix_type="S", name="(P|I_ab)", &
                        row_blk_size=ri_blk_size, col_blk_size=std_blk_size)
      CALL dbcsr_create(work2, dist=dbcsr_dist, matrix_type="N", name="WORK 2", &
                        row_blk_size=ri_blk_size, col_blk_size=std_blk_size)

      !reserve the blocks accordingly
      matrices(1)%matrix => aI_bP; matrices(2)%matrix => PI_ab
      CALL reserve_contraction_blocks(matrices, lb, ub, qs_env)

      matrices(1)%matrix => work1; matrices(2)%matrix => work2
      CALL reserve_contraction_blocks(matrices, lb, ub, qs_env , matrix_type="N")
      DEALLOCATE(matrices)

   !  Create the contracted integral matrices
      ALLOCATE(contr_int(ndo_so))
      DO i = 1, ndo_so
         ALLOCATE(contr_int(i)%matrix)
      END DO

!  Prepare the o3c_vec with the MO coeffs
      ALLOCATE(nsgf_orb(natom), o3c_vec(natom))
      nsgf_orb = 1
      nsgf_orb(lb:ub) = std_blk_size(lb:ub)
      CALL o3c_vec_create(o3c_vec, nsgf_orb)

   ! Only take the coeffs for atoms that are present in o3c_coul
      sr = 1; IF (lb > 1) sr = SUM(std_blk_size(1:lb-1))+1
      er = SUM(std_blk_size(1:ub))
      ALLOCATE(coeffs(er-sr+1,ndo_mo))

      DO ispin = 1, nspins
         CALL cp_fm_get_submatrix(fm=donor_state%gs_coeffs, target_m=coeffs, start_row=sr, &
                                  n_rows=er-sr+1, start_col=(ispin-1)*ndo_mo+1 , n_cols=ndo_mo)

!     Loop over the donor MOs, fill the o3c_vec and contract
         DO imo = 1, ndo_mo

            DO iat = lb,ub
               sr2 = 1; IF(iat > lb) sr2 = SUM(std_blk_size(lb:iat-1))+1
               er2 = SUM(std_blk_size(lb:iat))

               CALL get_o3c_vec(o3c_vec, iat, v)
               v(:) = coeffs(sr2:er2,imo)
            END DO !iat

            ! do the contraction
            CALL dbcsr_set(aI_bP, 0.0_dp); CALL dbcsr_set(PI_ab, 0.0_dp)
            IF (my_op == "C") THEN
               CALL contract_o3c_once(o3c_coul, o3c_vec, aI_bP, PI_ab, katom)
            ELSE
               CALL contract_o3c_once(o3c_ex, o3c_vec, aI_bP, PI_ab, katom)
            END IF
         
            ! Get the full "normal" (aI_b|P) contracted integrals
            CALL change_dist(PI_ab, work2, para_env)
            CALL dbcsr_transposed(contr_int((ispin-1)*ndo_mo+imo)%matrix, work2)
            CALL change_dist(aI_bP, work1, para_env)
            CALL dbcsr_add(contr_int((ispin-1)*ndo_mo+imo)%matrix, work1, 1.0_dp, 1.0_dp)
            
         END DO !imo
      END DO !ispin
      
!  Clean-up
      CALL o3c_vec_release(o3c_vec)
      CALL dbcsr_release(aI_bP)
      CALL dbcsr_release(PI_ab)
      CALL dbcsr_release(work1)
      CALL dbcsr_release(work2)
      DEALLOCATE(aI_bP, PI_ab, o3c_vec, ri_blk_size, work1, work2)

      CALL timestop(handle)

   END SUBROUTINE contract_o3c_int
   
! **************************************************************************************************
!> \brief Home made utility to transfer a matrix from a given proc distribution to another
!> \param mat_in ...
!> \param mat_out ...
!> \param para_env ...
!> \note It is assumed both matrices have the same block sizes and structure, and that all blocks
!>       are at least reserved
! **************************************************************************************************
   SUBROUTINE change_dist(mat_in, mat_out, para_env)

      TYPE(dbcsr_type), INTENT(INOUT)                 :: mat_in, mat_out
      TYPE(cp_para_env_type), POINTER                 :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = "change_dist", routineP = moduleN//":"//routineN

      INTEGER                                         :: handle, iblk, jblk, blk, s1, s2, dest, &
                                                         tag, nblk, group, source, num_pe, &
                                                         is, ir
      TYPE(dbcsr_iterator_type)                       :: iter
      REAL(dp), DIMENSION(:,:), POINTER               :: pbin, pbout
      LOGICAL                                         :: found_in, found_out
      TYPE(cp_2d_r_p_type), DIMENSION(:), ALLOCATABLE :: send_buff, recv_buff
      INTEGER, DIMENSION(:), ALLOCATABLE              :: send_req, recv_req

      NULLIFY(pbin, pbout)

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat_in, nblkrows_total=nblk)
      group = para_env%group
      num_pe = para_env%num_pe

      !Allocate a pointer array which will point on the block => we won't need to send more than
      !nblk**2/(2*num_pe) message per processor (only look at upper diagonal, blocks well distributed)
      ALLOCATE(send_buff(nblk**2/(2*num_pe)+1), recv_buff(nblk**2/(2*num_pe)+1))
      ALLOCATE(send_req(nblk**2/(2*num_pe)+1), recv_req(nblk**2/(2*num_pe)+1))
      is = 0; ir = 0

      !Iterate over input matrix
      CALL dbcsr_iterator_start(iter, mat_in)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         
         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
         CALL dbcsr_get_block_p(mat_in, iblk, jblk, pbin, found_in)

         IF (found_in) THEN

            !Check if the block is on the same proc on mat_out
            CALL dbcsr_get_block_p(mat_out, iblk, jblk, pbout, found_out)
            IF (found_out) THEN
               !Simply copy the block from mat_in to mat_out
               s1 = SIZE(pbout,1); s2 = SIZE(pbout,2)
               CALL dcopy(s1*s2, pbin, 1, pbout, 1)
            ELSE
               !If block not on the same processor, need to send it
               CALL dbcsr_get_stored_coordinates(mat_out, iblk, jblk, dest) 
               !unique tag
               tag = nblk*iblk + jblk
               !point on the block such that the buffer is not changed after isend
               is = is + 1
               send_buff(is)%array => pbin
               CALL mp_isend(msgin=send_buff(is)%array, dest=dest, comm=group, request=send_req(is), &
                             tag=tag)
            END IF

         END IF !found_in

      END DO !iterator
      CALL dbcsr_iterator_stop(iter)

      !Iterate over the output matrix
      CALL dbcsr_iterator_start(iter, mat_out)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
         CALL dbcsr_get_block_p(mat_out, iblk, jblk, pbout, found_out)

         IF (found_out) THEN
            !Check if the block is on the same proc on mat_in
            CALL dbcsr_get_block_p(mat_in, iblk, jblk, pbin, found_in)
            !If not, need to receiv it
            IF (.NOT. found_in) THEN
               CALL dbcsr_get_stored_coordinates(mat_in, iblk, jblk, source)
               tag = nblk*iblk + jblk
               ir = ir + 1
               recv_buff(ir)%array => pbout
               CALL mp_irecv(msgout=recv_buff(ir)%array, source=source, request=recv_req(ir), &
                             tag=tag, comm=group)
            END IF
         END IF

      END DO !iterator
      CALL dbcsr_iterator_stop(iter)

      !Make sure all communication is over. is, ir are 0 is no communication
      CALL mp_waitall(send_req(1:is))
      CALL mp_waitall(recv_req(1:ir))

      CALL timestop(handle)

   END SUBROUTINE change_dist
   
! **************************************************************************************************
!> \brief Reserves the blocks in of a dbcsr matrix as needed for RI 3-center contraction (aI_b|P)
!> \param matrices the matrices for which blocks are reserved
!> \param lb the lower bond of atoms a,b interaction with basis P
!> \param ub the upper bond of atoms a,b interaction with basis P
!> \param qs_env ...
!> \param matrix_type whether the matrices are symmetric or not
!> \param block_content in the case where all the blocks have the same dimension, can fill them
!>        with the values taken from this array
!> \note Even if the matrix_type is set to "normal", only reserve blocks on the upper-diagonal
! **************************************************************************************************
   SUBROUTINE reserve_contraction_blocks(matrices, lb, ub, qs_env, matrix_type, block_content)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrices
      INTEGER, INTENT(IN)                             :: lb, ub
      TYPE(qs_environment_type), POINTER              :: qs_env
      CHARACTER, INTENT(IN), OPTIONAL                 :: matrix_type
      REAL(dp), DIMENSION(:,:), INTENT(IN), OPTIONAL  :: block_content

      CHARACTER(len=*), PARAMETER :: routineN = "reserve_contraction_blocks", &
                                     routineP = moduleN//":"//routineN


      INTEGER                                         :: iblk, jblk, blk, s1, s2, i, n
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(dbcsr_type), POINTER                       :: template
      TYPE(dbcsr_iterator_type)                       :: iter
      LOGICAL                                         :: found
      REAL(dp), DIMENSION(:,:), POINTER               :: pblock_s, pblock_m
      CHARACTER                                       :: my_type

      NULLIFY(matrix_s, pblock_s, pblock_m, template)

!  Initialization
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      IF (PRESENT(block_content)) THEN
         s1 = SIZE(block_content,1)
         s2 = SIZE(block_content,2)
      END IF
      n = SIZE(matrices)
      !assume symmetric type
      my_type = dbcsr_type_symmetric
      IF (PRESENT(matrix_type)) my_type = matrix_type

!  If non-symmetric, need to desymmetrize matrix as as a template
      IF (my_type .NE. dbcsr_type_symmetric) THEN
         ALLOCATE(template)
         CALL dbcsr_desymmetrize(matrix_s(1)%matrix, template)
      ELSE
         template => matrix_s(1)%matrix
      END IF

!  Loop over matrix_s as need a,b to overlap
      CALL dbcsr_iterator_start(iter, template)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
         IF (iblk < lb .OR. jblk < lb) CYCLE
         IF (iblk > ub .OR. jblk > ub) CYCLE
         IF (jblk < iblk) CYCLE 
         
         CALL dbcsr_get_block_p(template, iblk, jblk, pblock_s, found) 

         IF (found) THEN
            DO i = 1,n
               NULLIFY(pblock_m)
               CALL dbcsr_reserve_block2d(matrices(i)%matrix, iblk, jblk, pblock_m)
               IF (PRESENT(block_content)) THEN
                  CALL dcopy(s1*s2, block_content, 1, pblock_m, 1)
               ELSE
                  pblock_m = 0.0_dp
               END IF
            END DO
         END IF

      END DO !dbcsr iter
      CALL dbcsr_iterator_stop(iter)
      DO i = 1,n
         CALL dbcsr_finalize(matrices(i)%matrix)
      END DO

!  Clean-up
      IF (my_type .NE. dbcsr_type_symmetric) THEN
         CALL dbcsr_release(template)
         DEALLOCATE(template)
      END IF

   END SUBROUTINE reserve_contraction_blocks
   
! **************************************************************************************************
!> \brief Contraction of the 3-center integrals over index 1 and 2, for a given atom_k. The results
!>        are stored in two matrices, such that (a,b are block indices):
!>        mat_aIb(ab) = mat_aIb(ab) + sum j_b (i_aj_b|k)*v(j_b) and
!>        mat_bIa(ba) = mat_bIa(ba) + sum i_a (i_aj_b|k)*v(i_a)
!>        The block size of the columns of mat_aIb and the rows of mat_bIa are the size of k
!> \param o3c the container for the integrals
!> \param vec the contraction coefficients
!> \param mat_aIb ...
!> \param mat_bIa ...
!> \param atom_k the atom for which we contract
!> \note To get the full (aI|P) matrix, one has to sum mat_aIb + mat_bIa^T (the latter has empty diag)
! **************************************************************************************************
   SUBROUTINE contract_o3c_once(o3c, vec, mat_aIb, mat_bIa, atom_k)
      TYPE(o3c_container_type), POINTER                  :: o3c
      TYPE(o3c_vec_type), DIMENSION(:), POINTER          :: vec
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_aIb, mat_bIa
      INTEGER, INTENT(IN)                                :: atom_k

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_o3c_once', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, icol, i, irow, jatom, &
                                                            katom, mepos, n, nthread, s1, s2
      LOGICAL                                            :: found, ijsymmetric, trans
      REAL(dp), DIMENSION(:), POINTER                    :: v
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock
      REAL(dp), DIMENSION(:, :), ALLOCATABLE             :: work
      REAL(dp), DIMENSION(:, :, :), POINTER              :: iabc
      TYPE(o3c_iterator_type)                            :: o3c_iterator
      INTEGER, DIMENSION(:), POINTER                     :: atom_blk_size

      CALL timeset(routineN, handle)

      CALL get_o3c_container(o3c, ijsymmetric=ijsymmetric)
      CPASSERT(ijsymmetric)

      CALL dbcsr_get_info(mat_aIb, row_blk_size=atom_blk_size)

      nthread = 1
!$    nthread = omp_get_max_threads()
      CALL o3c_iterator_create(o3c, o3c_iterator, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (nthread,o3c_iterator,vec,mat_aIb,mat_bIa,atom_k,atom_blk_size)&
!$OMP PRIVATE (mepos,iabc,iatom,jatom,katom,irow,icol,trans,pblock,v,found,i,n,s1,s2,work)

      mepos = 0
!$    mepos = omp_get_thread_num()

      DO WHILE (o3c_iterate(o3c_iterator, mepos=mepos) == 0)
         CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, iatom=iatom, jatom=jatom, katom=katom, &
                                    integral=iabc)

         IF (katom .NE. atom_k) CYCLE

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
            trans = .FALSE.
         ELSE
            irow = jatom
            icol = iatom
            trans = .TRUE.
         END IF

      ! Deal with mat_aIb
         CALL get_o3c_vec(vec, icol, v)
         n = SIZE(v)
         s1 = atom_blk_size(irow)
         s2 = SIZE(iabc,3)
         ALLOCATE(work(s1,s2))
         work = 0.0_dp

         IF (trans) THEN
            DO i = 1, n
               CALL daxpy(s1*s2, v(i), iabc(i, :, :), 1, work(:, :), 1)
            END DO
         ELSE
            DO i = 1, n
               CALL daxpy(s1*s2, v(i), iabc(:, i, :), 1, work(:, :), 1)
            END DO
         END IF

!$OMP CRITICAL(mat_aIb)
         CALL dbcsr_get_block_p(matrix=mat_aIb, row=irow, col=icol, BLOCK=pblock, found=found)
         IF (found) THEN
            CALL daxpy(s1*s2, 1.0_dp, work(:,:), 1, pblock(:, :), 1)
         END IF
!$OMP END CRITICAL(mat_aIb) 
         DEALLOCATE(work)

         ! Deal with mat_bIa, keep block diagonal empty
         IF (irow == icol) CYCLE

         CALL get_o3c_vec(vec, irow, v)
         n = SIZE(v)
         s1 = SIZE(iabc,3)
         s2 = atom_blk_size(icol)
         ALLOCATE(work(s1,s2))
         work = 0.0_dp

         IF (trans) THEN
            DO i = 1, n
               CALL daxpy(s1*s2, v(i), TRANSPOSE(iabc(:, i, :)), 1, work(:, :), 1)
            END DO
         ELSE 
            DO i = 1, n
               CALL daxpy(s1*s2, v(i), TRANSPOSE(iabc(i, :, :)), 1, work(:, :), 1)
            END DO
         END IF

!$OMP CRITICAL(mat_bIa)
         CALL dbcsr_get_block_p(matrix=mat_bIa, row=irow, col=icol, BLOCK=pblock, found=found)
         IF (found) THEN
            CALL daxpy(s1*s2, 1.0_dp, work(:,:), 1, pblock(:, :), 1)
         END IF
!$OMP END CRITICAL(mat_bIa)
         DEALLOCATE(work)

      END DO !o3c iterator
!$OMP END PARALLEL
      CALL o3c_iterator_release(o3c_iterator)

      CALL timestop(handle)

   END SUBROUTINE contract_o3c_once
   
! **************************************************************************************************
!> \brief Multiply all the blocks of a contracted RI integral (aI|P) by a matrix of type (P|...|Q)
!> \param contr_int the integral array
!> \param PQ the smaller matrix to multiply all blocks
!> \note  It is assumed that all non-zero blocks have the same number of columns. Can pass partial 
!>        arrays, e.g. contr_int(1:3)
! **************************************************************************************************
   SUBROUTINE ri_all_blocks_mm(contr_int, PQ)
      
      TYPE(dbcsr_p_type), DIMENSION(:)                :: contr_int
      REAL(dp), DIMENSION(:,:), INTENT(IN)            :: PQ

      CHARACTER(len=*), PARAMETER :: routineN = "ri_all_blocks_mm", routineP = moduleN//":"//routineN

      INTEGER                                         :: imo, ndo_mo, iblk, jblk, blk, s1, s2
      LOGICAL                                         :: found
      REAL(dp), DIMENSION(:,:), POINTER               :: pblock
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: work
      TYPE(dbcsr_iterator_type)                       :: iter

      NULLIFY(pblock)

      ndo_mo = SIZE(contr_int)

      DO imo = 1, ndo_mo
         CALL dbcsr_iterator_start(iter, contr_int(imo)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
            CALL dbcsr_get_block_p(contr_int(imo)%matrix, iblk, jblk, pblock, found)

            IF (found) THEN
               s1 = SIZE(pblock,1)
               s2 = SIZE(pblock,2)
               ALLOCATE(work(s1, s2))
               CALL dgemm('N', 'N', s1, s2, s2, 1.0_dp, pblock, s1, PQ, s2, 0.0_dp, work, s1)
               CALL dcopy(s1*s2, work, 1, pblock, 1)
               DEALLOCATE(work)
            END IF

         END DO ! dbcsr iterator
         CALL dbcsr_iterator_stop(iter)
      END DO !imo

   END SUBROUTINE ri_all_blocks_mm
   
! **************************************************************************************************
!> \brief Copies an (partial) array of contracted RI integrals into anoter one
!> \param new_int where the copy is stored
!> \param ref_int what is copied
!> \note Allocate the matrices of new_int if not done already
! **************************************************************************************************
   SUBROUTINE copy_ri_contr_int(new_int, ref_int)

      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT) :: new_int
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)    :: ref_int

      CHARACTER(len=*), PARAMETER :: routineN = "copy_ri_contr_int", routineP = moduleN//":"//routineN

      INTEGER                                         :: ndo_so, iso

      CPASSERT(SIZE(new_int) == SIZE(ref_int))
      ndo_so = SIZE(ref_int)

      DO iso = 1, ndo_so
         IF (.NOT. ASSOCIATED(new_int(iso)%matrix)) ALLOCATE(new_int(iso)%matrix)
         CALL dbcsr_copy(new_int(iso)%matrix, ref_int(iso)%matrix)
      END DO

   END SUBROUTINE copy_ri_contr_int
   
! **************************************************************************************************
!> \brief Takes the product of contracted integrals and put them in a kernel matrix
!> \param kernel the matrix where the products are stored
!> \param lhs_int the left-hand side contracted integrals
!> \param rhs_int the right-hand side contracted integrals
!> \param quadrants on which quadrant(s) on the kernel matrix the product is stored
!> \param qs_env ...
!> \param mo_transpose whether the MO blocks should be transpose, i.e. (aI|Jb) => (aJ|Ib)
!> \Note It is assumed that the kerenl matrix is NOT symmetric
!>       There are three quadrants, corresponding to 1: the upper-left (diagonal), 2: the
!>       upper-right (off-diagonal) and 3: the lower-right (diagonal). 
!>       Need to finalize the kernel matrix after calling this routine (possibly multiple times)
! **************************************************************************************************
   SUBROUTINE ri_int_product(kernel, lhs_int, rhs_int, quadrants, qs_env, mo_transpose)

      TYPE(dbcsr_type), INTENT(INOUT)                 :: kernel
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)    :: lhs_int
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)    :: rhs_int
      TYPE(qs_environment_type), POINTER              :: qs_env
      LOGICAL, DIMENSION(3), INTENT(IN)               :: quadrants
      LOGICAL, INTENT(IN), OPTIONAL                   :: mo_transpose

      CHARACTER(len=*), PARAMETER :: routineN ="ri_int_product", routineP = moduleN//":"//routineN

      INTEGER                                         :: ndo_so, iso, jso, iblk, jblk, blk, nblk, &
                                                         i,j
      TYPE(dbcsr_type)                                :: prod
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(dbcsr_iterator_type)                       :: iter
      LOGICAL                                         :: found, my_mt
      REAL(dp), DIMENSION(:,:), POINTER               :: pblock

      NULLIFY(matrix_s, pblock)

!  Initialization
      CPASSERT(SIZE(lhs_int) == SIZE(rhs_int))
      CPASSERT(ANY(quadrants))
      ndo_so = SIZE(lhs_int)
      CALL get_qs_env(qs_env, matrix_s=matrix_s, natom=nblk)
      CALL dbcsr_create(prod, template=matrix_s(1)%matrix, matrix_type="N") 
      my_mt = .FALSE.
      IF (PRESENT(mo_transpose)) my_mt = mo_transpose

      ! The kernel matrix is symmetric (even if normal type) => only fill upper half on diagonal
      ! quadrants, but the whole thing on upper-right quadrant
      DO iso = 1, ndo_so
         DO jso = 1, ndo_so

            ! If on-diagonal quadrants only, can skip jso < iso
            IF (.NOT. quadrants(2) .AND. jso<iso) CYCLE

            i = iso; j = jso;
            IF (my_mt) THEN
               i = jso; j = iso;
            END IF

            ! Take the product lhs*rhs^T
            CALL dbcsr_multiply('N', 'T', 1.0_dp, lhs_int(i)%matrix, rhs_int(j)%matrix, &
                                0.0_dp, prod)

            ! Loop over blocks of prod and fill kernel matrix => ok cuz same (but replicated) dist
            CALL dbcsr_iterator_start(iter, prod)
            DO WHILE (dbcsr_iterator_blocks_left(iter))

               CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
               IF ((iso == jso .AND. jblk<iblk) .AND. .NOT. quadrants(2)) CYCLE

               CALL dbcsr_get_block_p(prod, iblk, jblk, pblock, found)

               IF (found) THEN

                  ! Case study on quadrant
                  !upper-left
                  IF (quadrants(1)) THEN
                     CALL dbcsr_put_block(kernel, (iso-1)*nblk+iblk, (jso-1)*nblk+jblk, pblock) 
                  END IF

                  !upper-right 
                  IF (quadrants(2)) THEN
                     CALL dbcsr_put_block(kernel, (iso-1)*nblk+iblk, (ndo_so+jso-1)*nblk+jblk, pblock)
                  END IF

                  !lower-right
                  IF (quadrants(3)) THEN
                     CALL dbcsr_put_block(kernel, (ndo_so+iso-1)*nblk+iblk, (ndo_so+jso-1)*nblk+jblk, pblock)
                  END IF

               END IF

            END DO ! dbcsr iterator
            CALL dbcsr_iterator_stop(iter)

         END DO !jso
      END DO !iso

!  Clean-up
      CALL dbcsr_release(prod)

   END SUBROUTINE ri_int_product
   
END MODULE xas_tdp_kernel
