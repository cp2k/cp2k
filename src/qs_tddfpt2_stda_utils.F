!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Simplified Tamm Dancoff approach (sTDA).
! **************************************************************************************************
MODULE qs_tddfpt2_stda_utils

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: stda_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_row_scale,&
                                              cp_fm_schur_product
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type,&
                                              cp_fm_vectorssum
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add_on_diag, dbcsr_create, dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, &
        dbcsr_get_block_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_release, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE ewald_environment_types,         ONLY: ewald_env_create,&
                                              ewald_env_get,&
                                              ewald_env_set,&
                                              ewald_environment_type,&
                                              read_ewald_section_tb
   USE ewald_methods_tb,                ONLY: tb_ewald_overlap,&
                                              tb_spme_evaluate
   USE ewald_pw_types,                  ONLY: ewald_pw_create,&
                                              ewald_pw_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE iterate_matrix,                  ONLY: matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: oorootpi
   USE message_passing,                 ONLY: mp_sum
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_tddfpt2_stda_types,           ONLY: stda_env_type
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_work_matrices
   USE scf_control_types,               ONLY: scf_control_type
   USE util,                            ONLY: get_limit
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_stda_utils'

   PUBLIC:: stda_init_matrices, stda_calculate_kernel

CONTAINS

! **************************************************************************************************
!> \brief Calculate sTDA matrices
!> \param qs_env ...
!> \param stda_kernel ...
!> \param sub_env ...
!> \param work ...
!> \param tddfpt_control ...
! **************************************************************************************************
   SUBROUTINE stda_init_matrices(qs_env, stda_kernel, sub_env, work, tddfpt_control)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_env_type)                                :: stda_kernel
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control

      CHARACTER(len=*), PARAMETER :: routineN = 'stda_init_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: do_coulomb
      TYPE(cell_type), POINTER                           :: cell, cell_ref
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(section_vals_type), POINTER                   :: ewald_section, poisson_section, &
                                                            print_section

      CALL timeset(routineN, handle)

      do_coulomb = .NOT. tddfpt_control%rks_triplets
      IF (do_coulomb) THEN
         ! calculate exchange gamma matrix
         CALL setup_gamma(qs_env, stda_kernel, sub_env, work%gamma_exchange)
      END IF

      ! calculate S_half and Lowdin MO coefficients
      CALL get_lowdin_mo_coefficients(qs_env, sub_env, work)

      ! initialize Ewald for sTDA
      IF (tddfpt_control%stda_control%do_ewald) THEN
         NULLIFY (ewald_env, ewald_pw)
         CALL ewald_env_create(ewald_env, sub_env%para_env)
         poisson_section => section_vals_get_subs_vals(qs_env%input, "DFT%POISSON")
         CALL ewald_env_set(ewald_env, poisson_section=poisson_section)
         ewald_section => section_vals_get_subs_vals(poisson_section, "EWALD")
         print_section => section_vals_get_subs_vals(qs_env%input, "PRINT%GRID_INFORMATION")
         CALL get_qs_env(qs_env, cell=cell, cell_ref=cell_ref)
         CALL read_ewald_section_tb(ewald_env, ewald_section, cell_ref%hmat)
         CALL ewald_pw_create(ewald_pw, ewald_env, cell, cell_ref, print_section=print_section)
         work%ewald_env => ewald_env
         work%ewald_pw => ewald_pw
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE stda_init_matrices
! **************************************************************************************************
!> \brief Calculate sTDA exchange-type contributions
!> \param qs_env ...
!> \param stda_env ...
!> \param sub_env ...
!> \param gamma_matrix sTDA exchange-type contributions
!> \note  Note the specific sTDA notation exchange-type integrals (ia|jb) refer to Coulomb interaction
! **************************************************************************************************
   SUBROUTINE setup_gamma(qs_env, stda_env, sub_env, gamma_matrix)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_env_type)                                :: stda_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gamma_matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'setup_gamma', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: rsmooth = 1.0_dp

      CHARACTER(len=20)                                  :: cstring
      INTEGER                                            :: handle, iatom, icol, ikind, irow, jatom, &
                                                            jkind, natom, nmat
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: dr, eta, fcut, r, ralpha, rcut, rcuta, &
                                                            rcutb, x
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gblock
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list

      CALL timeset(routineN, handle)

      cstring = "GAMMA EXCHANGE MATRIX"
      dbcsr_dist => sub_env%dbcsr_dist
      n_list => sub_env%sab_orb
      CALL get_qs_env(qs_env=qs_env, natom=natom)

      nmat = 1
      NULLIFY (gamma_matrix)
      CALL dbcsr_allocate_matrix_set(gamma_matrix, nmat)

      ALLOCATE (row_blk_sizes(natom))
      row_blk_sizes(1:natom) = 1
      ALLOCATE (gamma_matrix(1)%matrix)

      CALL dbcsr_create(gamma_matrix(1)%matrix, name=TRIM(ADJUSTL(cstring)), dist=dbcsr_dist, &
                        matrix_type=dbcsr_type_symmetric, row_blk_size=row_blk_sizes, &
                        col_blk_size=row_blk_sizes, nze=0, mutable_work=.TRUE.)

      DEALLOCATE (row_blk_sizes)

      ! setup the matrices using the neighbor list
      CALL cp_dbcsr_alloc_block_from_nbl(gamma_matrix(1)%matrix, n_list)
      CALL dbcsr_set(gamma_matrix(1)%matrix, 0.0_dp)

      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, n_list)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij)

         dr = SQRT(SUM(rij(:)**2)) ! interatomic distance

         eta = (stda_env%kind_param_set(ikind)%kind_param%hardness_param + &
                stda_env%kind_param_set(jkind)%kind_param%hardness_param)/2.0_dp

         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)

         NULLIFY (gblock)
         CALL dbcsr_get_block_p(matrix=gamma_matrix(1)%matrix, &
                                row=irow, col=icol, BLOCK=gblock, found=found)
         CPASSERT(found)

         ! get rcuta and rcutb
         rcuta = stda_env%kind_param_set(ikind)%kind_param%rcut
         rcutb = stda_env%kind_param_set(jkind)%kind_param%rcut
         rcut = rcuta + rcutb

         !>   Computes the short-range gamma parameter from
         !>   Nataga-Mishimoto-Ohno-Klopman formula equivalently as it is done for xTB
         IF (dr < 1.e-6) THEN
            ! on site terms
            gblock(:, :) = gblock(:, :) + eta
         ELSEIF (dr > rcut) THEN
            ! do nothing
         ELSE
            IF (dr < rcut - rsmooth) THEN
               fcut = 1.0_dp
            ELSE
               r = dr - (rcut - rsmooth)
               x = r/rsmooth
               fcut = -6._dp*x**5 + 15._dp*x**4 - 10._dp*x**3 + 1._dp
            END IF
            ralpha = dr**stda_env%alpha_param
            gblock(:, :) = gblock(:, :) + &
                           fcut*(1._dp/(dr**(stda_env%alpha_param) + eta**(-stda_env%alpha_param))) &
                           **(1._dp/stda_env%alpha_param) - fcut/dr
         END IF

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)
      CALL dbcsr_finalize(gamma_matrix(1)%matrix)

      CALL timestop(handle)

   END SUBROUTINE setup_gamma

! **************************************************************************************************
!> \brief Calculate Lowdin MO coefficients
!> \param qs_env ...
!> \param sub_env ...
!> \param work ...
! **************************************************************************************************
   SUBROUTINE get_lowdin_mo_coefficients(qs_env, sub_env, work)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work

      CHARACTER(len=*), PARAMETER :: routineN = 'get_lowdin_mo_coefficients', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iounit, ispin, max_iter_lanczos, &
                                                            nactive, ndep, nsgf, nspins, &
                                                            order_lanczos
      LOGICAL                                            :: converged
      REAL(KIND=dp)                                      :: eps_lanczos, threshold
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type), POINTER                          :: fm_s_half, fm_work1
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrixkp_s
      TYPE(dbcsr_type)                                   :: sm_hinv
      TYPE(dbcsr_type), POINTER                          :: sm_h, sm_s
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (logger) !get output_unit
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      ! Calculate S^1/2 matrix
      IF (iounit > 0) THEN
         WRITE (iounit, "(1X,A)") "", &
            "-------------------------------------------------------------------------------", &
            "-                      sTDA Kernel: Create Matrix SQRT(S)                     -", &
            "-------------------------------------------------------------------------------"
      END IF

      IF (sub_env%is_split) THEN
         CPABORT('SPLIT')
      ELSE
         CALL get_qs_env(qs_env=qs_env, matrix_s_kp=matrixkp_s)
         CPASSERT(ASSOCIATED(matrixkp_s))
         IF (SIZE(matrixkp_s, 2) > 1) CPWARN("not implemented for k-points.")
         sm_s => matrixkp_s(1, 1)%matrix
      END IF
      sm_h => work%shalf

      CALL dbcsr_create(sm_hinv, template=sm_s)
      CALL dbcsr_add_on_diag(sm_h, 1.0_dp)
      threshold = 1.0e-8_dp
      order_lanczos = 3
      eps_lanczos = 1.0e-4_dp
      max_iter_lanczos = 40
      CALL matrix_sqrt_Newton_Schulz(sm_h, sm_hinv, sm_s, &
                                     threshold, order_lanczos, eps_lanczos, max_iter_lanczos, &
                                     converged=converged)
      CALL dbcsr_release(sm_hinv)
      !
      NULLIFY (qs_kind_set)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      ! Get the total number of contracted spherical Gaussian basis functions
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)
      !
      IF (.NOT. converged) THEN
         IF (iounit > 0) THEN
            WRITE (iounit, "(T3,A)") "STDA| Newton-Schulz iteration did not converge"
            WRITE (iounit, "(T3,A)") "STDA| Calculate SQRT(S) from diagonalization"
         END IF
         CALL get_qs_env(qs_env=qs_env, scf_control=scf_control)
         ! Provide full size work matrices
         CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                  para_env=sub_env%para_env, &
                                  context=sub_env%blacs_env, &
                                  nrow_global=nsgf, &
                                  ncol_global=nsgf)
         CALL cp_fm_create(matrix=fm_s_half, matrix_struct=fmstruct, name="S^(1/2) MATRIX")
         CALL cp_fm_create(matrix=fm_work1, matrix_struct=fmstruct, name="TMP MATRIX")
         CALL cp_fm_struct_release(fmstruct=fmstruct)
         CALL copy_dbcsr_to_fm(sm_s, fm_s_half)
         CALL cp_fm_power(fm_s_half, fm_work1, 0.5_dp, scf_control%eps_eigval, ndep)
         IF (ndep /= 0) &
            CALL cp_warn(__LOCATION__, &
                         "Overlap matrix exhibits linear dependencies. At least some "// &
                         "eigenvalues have been quenched.")
         CALL copy_fm_to_dbcsr(fm_s_half, sm_h)
         CALL cp_fm_release(matrix=fm_s_half)
         CALL cp_fm_release(matrix=fm_work1)
         IF (iounit > 0) WRITE (iounit, *)
      END IF

      nspins = SIZE(sub_env%mos_occ)

      DO ispin = 1, nspins
         CALL cp_fm_get_info(work%ctransformed(ispin)%matrix, ncol_global=nactive)
         CALL cp_dbcsr_sm_fm_multiply(work%shalf, sub_env%mos_occ(ispin)%matrix, &
                                      work%ctransformed(ispin)%matrix, nactive, alpha=1.0_dp, beta=0.0_dp)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE get_lowdin_mo_coefficients

! **************************************************************************************************
!> \brief Calculate Lowdin transformed Davidson trial vector X
!>        shalf (dbcsr), xvec, xt (fm) are defined in the same sub_env
!> \param shalf ...
!> \param xvec ...
!> \param xt ...
! **************************************************************************************************
   SUBROUTINE get_lowdin_x(shalf, xvec, xt)

      TYPE(dbcsr_type)                                   :: shalf
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: xvec
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(INOUT)    :: xt

      CHARACTER(len=*), PARAMETER :: routineN = 'get_lowdin_x', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nactive, nspins

      CALL timeset(routineN, handle)

      nspins = SIZE(xvec)

      ! Build Lowdin transformed tilde(X)= S^1/2 X for each spin
      DO ispin = 1, nspins
         CALL cp_fm_get_info(xt(ispin)%matrix, ncol_global=nactive)
         CALL cp_dbcsr_sm_fm_multiply(shalf, xvec(ispin)%matrix, &
                                      xt(ispin)%matrix, nactive, alpha=1.0_dp, beta=0.0_dp)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE get_lowdin_x

! **************************************************************************************************
!> \brief ...Calculate the sTDA kernel contribution by contracting the Lowdin MO coefficients --
!>           transition charges with the Coulomb-type or exchange-type integrals
!> \param qs_env ...
!> \param stda_control ...
!> \param stda_env ...
!> \param sub_env ...
!> \param work ...
!> \param is_rks_triplets ...
!> \param X ...
!> \param res ... vector AX with A being the sTDA matrix and X the Davidson trial vector of the
!>                eigenvalue problem A*X = omega*X
! **************************************************************************************************
   SUBROUTINE stda_calculate_kernel(qs_env, stda_control, stda_env, sub_env, &
                                    work, is_rks_triplets, X, res)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_control_type)                            :: stda_control
      TYPE(stda_env_type)                                :: stda_env
      TYPE(tddfpt_subgroup_env_type)                     :: sub_env
      TYPE(tddfpt_work_matrices)                         :: work
      LOGICAL, INTENT(IN)                                :: is_rks_triplets
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: X
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(INOUT)    :: res

      CHARACTER(len=*), PARAMETER :: routineN = 'stda_calculate_kernel', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, ewald_type, handle, ia, iatom, &
                                                            ikind, is, ispin, jatom, jkind, jspin, &
                                                            natom, nsgf, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, kind_of, last_sgf
      INTEGER, DIMENSION(2)                              :: nactive, nlim
      LOGICAL                                            :: calculate_forces, do_coulomb, do_ewald, &
                                                            use_virial
      REAL(KIND=dp)                                      :: alpha, bp, dr, eta, gabr, hfx, rbeta, &
                                                            spinfac
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tcharge, tv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gtcharge
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gab, pblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: xtransformed
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct, fmstructjspin
      TYPE(cp_fm_type), POINTER                          :: ct, ctjspin, cvec, cvecjspin, xt
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: pdens
      TYPE(dbcsr_type), POINTER                          :: tempmat
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      nactive(:) = stda_env%nactive(:)
      nspins = SIZE(X)
      spinfac = 2.0_dp
      IF (nspins == 2) spinfac = 1.0_dp

      IF (nspins == 1 .AND. is_rks_triplets) THEN
         do_coulomb = .FALSE.
      ELSE
         do_coulomb = .TRUE.
      END IF
      do_ewald = stda_control%do_ewald

      para_env => sub_env%para_env

      CALL get_qs_env(qs_env, natom=natom, cell=cell, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set)
      ALLOCATE (first_sgf(natom))
      ALLOCATE (last_sgf(natom))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, last_sgf=last_sgf)

      ! calculate Loewdin transformed Davidson trial vector tilde(X)=S^1/2*X
      ! and tilde(tilde(X))=S^1/2_A*tilde(X)_A
      ALLOCATE (xtransformed(nspins))
      DO ispin = 1, nspins
         NULLIFY (fmstruct)
         ct => work%ctransformed(ispin)%matrix
         CALL cp_fm_get_info(ct, matrix_struct=fmstruct)
         CALL cp_fm_create(matrix=xtransformed(ispin)%matrix, matrix_struct=fmstruct, name="XTRANSFORMED")
      ENDDO
      CALL get_lowdin_x(work%shalf, X, xtransformed)

      ALLOCATE (tcharge(natom), gtcharge(natom, 1))

      DO ispin = 1, nspins
         CALL cp_fm_set_all(res(ispin)%matrix, 0.0_dp)
      END DO

      DO ispin = 1, nspins
         ct => work%ctransformed(ispin)%matrix
         CALL cp_fm_get_info(ct, matrix_struct=fmstruct, nrow_global=nsgf)
         ALLOCATE (tv(nsgf))
         CALL cp_fm_create(cvec, fmstruct)
         !
         ! *** Coulomb contribution
         !
         IF (do_coulomb) THEN
            tcharge(:) = 0.0_dp
            DO jspin = 1, nspins
               ctjspin => work%ctransformed(jspin)%matrix
               xt => xtransformed(jspin)%matrix
               CALL cp_fm_get_info(ctjspin, matrix_struct=fmstructjspin)
               CALL cp_fm_get_info(ctjspin, matrix_struct=fmstructjspin, nrow_global=nsgf)
               CALL cp_fm_create(cvecjspin, fmstructjspin)
               ! CV(mu,j) = CT(mu,j)*XT(mu,j)
               CALL cp_fm_schur_product(ctjspin, xt, cvecjspin)
               ! TV(mu) = SUM_j CV(mu,j)
               CALL cp_fm_vectorssum(cvecjspin, tv, "R")
               ! contract charges
               ! TC(a) = SUM_(mu of a) TV(mu)
               DO ia = 1, natom
                  DO is = first_sgf(ia), last_sgf(ia)
                     tcharge(ia) = tcharge(ia) + tv(is)
                  END DO
               END DO
               CALL cp_fm_release(cvecjspin)
            ENDDO !jspin
            ! Apply tcharge*gab -> gtcharge
            ! gT(b) = SUM_a g(a,b)*TC(a)
            ! gab = work%gamma_exchange(1)%matrix
            gtcharge = 0.0_dp
            ! short range contribution
            tempmat => work%gamma_exchange(1)%matrix
            CALL dbcsr_iterator_start(iter, tempmat)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, iatom, jatom, gab, blk)
               gtcharge(iatom, 1) = gtcharge(iatom, 1) + gab(1, 1)*tcharge(jatom)
               IF (iatom /= jatom) THEN
                  gtcharge(jatom, 1) = gtcharge(jatom, 1) + gab(1, 1)*tcharge(iatom)
               END IF
            END DO
            CALL dbcsr_iterator_stop(iter)
            ! Ewald long range contribution
            IF (do_ewald) THEN
               ewald_env => work%ewald_env
               ewald_pw => work%ewald_pw
               CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type)
               CALL get_qs_env(qs_env=qs_env, virial=virial, atprop=atprop)
               use_virial = .FALSE.
               calculate_forces = .FALSE.
               n_list => sub_env%sab_orb
               CALL tb_ewald_overlap(gtcharge, tcharge, alpha, n_list, virial, use_virial, atprop)
               CALL tb_spme_evaluate(ewald_env, ewald_pw, particle_set, cell, &
                                     gtcharge, tcharge, calculate_forces, virial, use_virial, atprop)
               ! add self charge interaction contribution
               IF (para_env%ionode) THEN
                  gtcharge(:, 1) = gtcharge(:, 1) - 2._dp*alpha*oorootpi*tcharge(:)
               END IF
            ELSE
               nlim = get_limit(natom, para_env%num_pe, para_env%mepos)
               DO iatom = nlim(1), nlim(2)
                  DO jatom = 1, iatom - 1
                     rij = particle_set(iatom)%r - particle_set(jatom)%r
                     rij = pbc(rij, cell)
                     dr = SQRT(SUM(rij(:)**2))
                     IF (dr > 1.e-6_dp) THEN
                        gtcharge(iatom, 1) = gtcharge(iatom, 1) + tcharge(jatom)/dr
                        gtcharge(jatom, 1) = gtcharge(jatom, 1) + tcharge(iatom)/dr
                     END IF
                  END DO
               END DO
            END IF
            CALL mp_sum(gtcharge, para_env%group)
            ! expand charges
            ! TV(mu) = TC(a of mu)
            tv(1:nsgf) = 0.0_dp
            DO ia = 1, natom
               DO is = first_sgf(ia), last_sgf(ia)
                  tv(is) = gtcharge(ia, 1)
               END DO
            END DO
            ! CV(mu,i) = TV(mu)*CV(mu,i)
            ct => work%ctransformed(ispin)%matrix
            CALL cp_fm_to_fm(ct, cvec)
            CALL cp_fm_row_scale(cvec, tv)
            ! rho(nu,i) = rho(nu,i) + Shalf(nu,mu)*CV(mu,i)
            CALL cp_dbcsr_sm_fm_multiply(work%shalf, cvec, res(ispin)%matrix, nactive(ispin), spinfac, 1.0_dp)
         END IF
         !
         ! *** Exchange contribution
         !
         IF (stda_env%hfx_fraction > 0.0_dp) THEN
            ALLOCATE (kind_of(natom))
            CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
            CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
            !
            tempmat => work%shalf
            CALL dbcsr_create(pdens, template=tempmat, matrix_type=dbcsr_type_no_symmetry)
            ! P(nu,mu) = SUM_j XT(nu,j)*CT(mu,j)
            ct => work%ctransformed(ispin)%matrix
            xt => xtransformed(ispin)%matrix
            CALL cp_dbcsr_plus_fm_fm_t(pdens, xt, ct, nactive(ispin), &
                                       1.0_dp, keep_sparsity=.FALSE.)
            CALL dbcsr_filter(pdens, stda_env%eps_td_filter)
            ! Apply PP*gab -> PP; gab = gamma_coulomb
            ! P(nu,mu) = P(nu,mu)*g(a of nu,b of mu)
            bp = stda_env%beta_param
            hfx = stda_env%hfx_fraction
            CALL dbcsr_iterator_start(iter, pdens)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, iatom, jatom, pblock, blk)
               rij = particle_set(iatom)%r - particle_set(jatom)%r
               rij = pbc(rij, cell)
               dr = SQRT(SUM(rij(:)**2))
               ikind = kind_of(iatom)
               jkind = kind_of(jatom)
               eta = (stda_env%kind_param_set(ikind)%kind_param%hardness_param + &
                      stda_env%kind_param_set(jkind)%kind_param%hardness_param)/2.0_dp
               rbeta = dr**bp
               gabr = (1._dp/(rbeta + (hfx*eta)**(-bp)))**(1._dp/bp)
               pblock = gabr*pblock
            END DO
            CALL dbcsr_iterator_stop(iter)
            ! CV(mu,i) = P(nu,mu)*CT(mu,i)
            CALL cp_dbcsr_sm_fm_multiply(pdens, ct, cvec, nactive(ispin), 1.0_dp, 0.0_dp)
            ! rho(nu,i) = rho(nu,i) + ShalfP(nu,mu)*CV(mu,i)
            CALL cp_dbcsr_sm_fm_multiply(work%shalf, cvec, res(ispin)%matrix, nactive(ispin), -1.0_dp, 1.0_dp)
            !
            CALL dbcsr_release(pdens)
            DEALLOCATE (kind_of)
         END IF
         !
         CALL cp_fm_release(cvec)
         DEALLOCATE (tv)
      ENDDO

      DO ispin = 1, nspins
         CALL cp_fm_release(matrix=xtransformed(ispin)%matrix)
      END DO
      DEALLOCATE (xtransformed)
      DEALLOCATE (tcharge, gtcharge)
      DEALLOCATE (first_sgf, last_sgf)

      CALL timestop(handle)

   END SUBROUTINE stda_calculate_kernel

! **************************************************************************************************

END MODULE qs_tddfpt2_stda_utils
