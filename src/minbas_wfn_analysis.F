!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate localized minimal basis and analyze wavefunctions
!> \par History
!>      12.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE minbas_wfn_analysis
   USE atomic_charges,                  ONLY: print_atomic_charges,&
                                              print_bond_orders
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Lu2004,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_allocate_matrix_set, cp_dbcsr_copy, cp_dbcsr_create, &
        cp_dbcsr_deallocate_matrix_set, cp_dbcsr_filter, cp_dbcsr_get_block_p, &
        cp_dbcsr_get_occupation, cp_dbcsr_init, cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
        cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, &
        cp_dbcsr_multiply, cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_reserve_diag_blocks, &
        cp_dbcsr_set, cp_dbcsr_trace, cp_dbcsr_type, dbcsr_distribution_obj, &
        dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
   USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                              cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE iterate_matrix,                  ONLY: invert_Hotelling
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE lapack,                          ONLY: lapack_ssyev
   USE mao_basis,                       ONLY: mao_generate_basis
   USE message_passing,                 ONLY: mp_sum
   USE molden_utils,                    ONLY: write_mos_molden
   USE mulliken,                        ONLY: compute_bond_order,&
                                              mulliken_charges
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_mo_methods,                   ONLY: make_basis_lowdin
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type,&
                                              set_mo_set
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE block_type
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE  :: mat
   END TYPE block_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'minbas_wfn_analysis'

   PUBLIC ::  minbas_analysis

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param input_section ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE minbas_analysis(qs_env, input_section, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input_section
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=*), PARAMETER :: routineN = 'minbas_analysis', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, homo, i, iab, ispin, n2, nao, &
                                                            natom, ndep, nimages, nmao, nmo, np, &
                                                            np1, nspin, nvirt
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: ecount
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, row_blk_sizes
      LOGICAL                                            :: do_bondorder, explicit, full_ortho, occeq
      REAL(KIND=dp)                                      :: alpha, amax, eps_filter, filter_eps, &
                                                            trace
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dval, dvalo, dvalv, eigval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: border, fnorm, mcharge, prmao
      REAL(KIND=dp), DIMENSION(:), POINTER               :: occupation_numbers
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mao_coef, pqmat, quambo, sqmat
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: matrix_s
      TYPE(cp_dbcsr_type)                                :: psmat, sinv, smao, smaox, sortho, spmat
      TYPE(cp_dbcsr_type), POINTER                       :: smat
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_a, fm_struct_b, fm_struct_c, &
                                                            fm_struct_d, fm_struct_e
      TYPE(cp_fm_type), POINTER                          :: fm1, fm2, fm3, fm4, fm5, fm6, fm_mos, &
                                                            fma, fmb, fmwork
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_obj), POINTER              :: dbcsr_dist
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mbas, mos
      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

!     TYPE(mo_set_type), POINTER                         :: mo_set

      ! only do MINBAS analysis if explicitely requested
      CALL section_vals_get(input_section, explicit=explicit)
      IF (.NOT. explicit) RETURN

      ! k-points?
      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspin = dft_control%nspins
      nimages = dft_control%nimages
      IF (nimages > 1) THEN
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T2,A)") &
               "K-Points: Localized Minimal Basis Analysis not available."
         END IF
      END IF
      IF (nimages > 1) RETURN

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(/,T2,A)') '!-----------------------------------------------------------------------------!'
         WRITE (UNIT=unit_nr, FMT="(T26,A)") "LOCALIZED MINIMAL BASIS ANALYSIS"
         WRITE (UNIT=unit_nr, FMT="(T18,A)") "W.C. Lu et al, J. Chem. Phys. 120, 2629 (2004)"
         WRITE (unit_nr, '(T2,A)') '!-----------------------------------------------------------------------------!'
      END IF
      CALL cite_reference(Lu2004)

      ! input options
      CALL section_vals_val_get(input_section, "EPS_FILTER", r_val=eps_filter)
      CALL section_vals_val_get(input_section, "FULL_ORTHOGONALIZATION", l_val=full_ortho)
      CALL section_vals_val_get(input_section, "BOND_ORDER", l_val=do_bondorder)

      ! MAO matrices
      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, natom=natom)
      CALL get_ks_env(ks_env=ks_env, particle_set=particle_set, dbcsr_dist=dbcsr_dist)
      NULLIFY (mao_coef, quambo, pqmat, sqmat)
      CALL cp_dbcsr_allocate_matrix_set(quambo, nspin)
      CALL cp_dbcsr_allocate_matrix_set(sqmat, nspin)
      CALL cp_dbcsr_allocate_matrix_set(pqmat, nspin)
      ALLOCATE (row_blk_sizes(natom), col_blk_sizes(natom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes)
      CALL get_particle_set(particle_set, qs_kind_set, nmao=col_blk_sizes)
      ! check if MAOs have been specified
      DO iab = 1, natom
         IF (col_blk_sizes(iab) < 0) &
            CPABORT("Number of MAOs has to be specified in KIND section for all elements")
      END DO

      ! basis size, orbitals, maos
      CALL get_qs_env(qs_env, mos=mos)
      CALL get_mo_set(mo_set=mos(1)%mo_set, nao=nao, nmo=nmo)
      nmao = SUM(col_blk_sizes)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A,T71,I10)') 'Total Number of Atomic Basis Set Functions   :', nao
         WRITE (unit_nr, '(T2,A,T71,I10)') 'Total Number of Minimal Basis Set Functions  :', nmao
         IF (nspin == 1) THEN
            WRITE (unit_nr, '(T2,A,T71,I10)') 'Total Number of Molecular Orbitals available :', nmo
         ELSE
            DO ispin = 1, nspin
               CALL get_mo_set(mo_set=mos(ispin)%mo_set, nmo=n2)
               WRITE (unit_nr, '(T2,A,i2,T71,I10)') 'Total Number of Molecular Orbitals available for Spin ', ispin, n2
            END DO
         END IF
      END IF
      CPASSERT(nmao <= nao)

      DO ispin = 1, nspin
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, nmo=n2)
         IF (n2 /= nmo) EXIT
      END DO
      IF (nmao > nmo) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A)') 'Localized Minimal Basis Analysis not possible'
         END IF
      ELSEIF (nmo /= n2) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A)') 'Different Number of Alpha and Beta MOs'
            WRITE (unit_nr, '(T2,A)') 'Localized Minimal Basis Analysis not possible'
         END IF
      ELSE
         IF (nao > nmo) THEN
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A)') 'WARNING: Only a subset of MOs is available: Analysis depends on MOs'
            END IF
         END IF

         CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)
         NULLIFY (fm_struct_a, fm_struct_b)
         CALL cp_fm_struct_create(fm_struct_a, nrow_global=nao, ncol_global=nmao, &
                                  para_env=para_env, context=blacs_env)
         CALL cp_fm_struct_create(fm_struct_b, nrow_global=nmo, ncol_global=nmao, &
                                  para_env=para_env, context=blacs_env)
         CALL cp_fm_create(fm1, fm_struct_a)
         CALL cp_fm_create(fm2, fm_struct_b)
         CALL cp_fm_create(fma, fm_struct_b)
         CALL cp_fm_create(fmb, fm_struct_b)

         DO ispin = 1, nspin
            ! coeficients
            ALLOCATE (quambo(ispin)%matrix)
            CALL cp_dbcsr_init(quambo(ispin)%matrix)
            CALL cp_dbcsr_create(matrix=quambo(ispin)%matrix, &
                                 name="QUAMBO", dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                                 row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
            ALLOCATE (sqmat(ispin)%matrix)
            CALL cp_dbcsr_init(sqmat(ispin)%matrix)
            CALL cp_dbcsr_create(matrix=sqmat(ispin)%matrix, &
                                 name="SQMAT", dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                                 row_blk_size=col_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
            ALLOCATE (pqmat(ispin)%matrix)
            CALL cp_dbcsr_init(pqmat(ispin)%matrix)
            CALL cp_dbcsr_create(matrix=pqmat(ispin)%matrix, &
                                 name="PQMAT", dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                                 row_blk_size=col_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
         END DO
         CALL cp_dbcsr_init(sortho)
         CALL cp_dbcsr_create(sortho, name="SORTHO", dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                              row_blk_size=col_blk_sizes, col_blk_size=col_blk_sizes, nze=0)
         CALL cp_dbcsr_reserve_diag_blocks(matrix=sortho)
         DEALLOCATE (row_blk_sizes, col_blk_sizes)

         ! optimze MAOs (mao_coef is allocated in the routine)
         CALL mao_generate_basis(qs_env, mao_coef)

         CALL get_qs_env(qs_env, matrix_s_kp=matrix_s)
         smat => matrix_s(1, 1)%matrix
         DO ispin = 1, nspin
            CALL cp_dbcsr_init(smao)
            CALL cp_dbcsr_create(smao, name="S*MAO", template=mao_coef(1)%matrix)

            ! SMAO = Overlap*MAOs
            CALL cp_dbcsr_multiply("N", "N", 1.0_dp, smat, mao_coef(ispin)%matrix, 0.0_dp, smao)
            ! a(nj)* = C(vn)(T) * SMAO(vj)
            CALL copy_dbcsr_to_fm(smao, fm1)
            CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=fm_mos)
            CALL cp_gemm("T", "N", nmo, nmao, nao, 1.0_dp, fm_mos, fm1, 0.0_dp, fm2)
            CALL cp_dbcsr_release(smao)
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, homo=homo)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,T51,A,i2,T71,I10)') 'MOs in Occupied Valence Set', 'Spin ', ispin, homo
            END IF
            nvirt = nmo-homo
            NULLIFY (fm_struct_c)
            CALL cp_fm_struct_create(fm_struct_c, nrow_global=nvirt, ncol_global=nvirt, &
                                     para_env=para_env, context=blacs_env)
            CALL cp_fm_create(fm3, fm_struct_c)
            CALL cp_fm_create(fm4, fm_struct_c)
            ! B(vw) = a(vj)* * a(wj)*
            CALL cp_gemm("N", "T", nvirt, nvirt, nmao, 1.0_dp, fm2, fm2, 0.0_dp, fm3, &
                         a_first_row=homo+1, b_first_row=homo+1)
            ALLOCATE (eigval(nvirt))
            CALL choose_eigv_solver(fm3, fm4, eigval)
            ! SVD(B) -> select p largest eigenvalues and vectors
            np = nmao-homo
            np1 = nvirt-np+1
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,T51,A,i2,T71,I10)') 'MOs in Virtual Valence Set', 'Spin ', ispin, np
            END IF
            ! R(vw) = SUM_p T(vp)*T(wp)
            CALL cp_gemm("N", "T", nvirt, nvirt, np, 1.0_dp, fm4, fm4, 0.0_dp, fm3, &
                         a_first_col=np1, b_first_col=np1)
            !
            ALLOCATE (dval(nmao), dvalo(nmao), dvalv(nmao))
            NULLIFY (fm_struct_d)
            CALL cp_fm_struct_create(fm_struct_d, nrow_global=nvirt, ncol_global=nmao, &
                                     para_env=para_env, context=blacs_env)
            CALL cp_fm_create(fm5, fm_struct_d)
            NULLIFY (fm_struct_e)
            CALL cp_fm_struct_create(fm_struct_e, nrow_global=nmao, ncol_global=nmao, &
                                     para_env=para_env, context=blacs_env)
            CALL cp_fm_create(fm6, fm_struct_e)
            ! D(j) = SUM_n (a(nj)*)^2 + SUM_vw R(vw) * a(vj)* * a(wj)*
            CALL cp_gemm("N", "N", nvirt, nmao, nvirt, 1.0_dp, fm3, fm2, 0.0_dp, fm5, &
                         b_first_row=homo+1)
            CALL cp_gemm("T", "N", nmao, nmao, nvirt, 1.0_dp, fm2, fm5, 0.0_dp, fm6, &
                         a_first_row=homo+1)
            CALL cp_fm_get_diag(fm6, dvalv(1:nmao))
            CALL cp_gemm("T", "N", nmao, nmao, homo, 1.0_dp, fm2, fm2, 0.0_dp, fm6)
            CALL cp_fm_get_diag(fm6, dvalo(1:nmao))
            DO i = 1, nmao
               dval(i) = 1.0_dp/SQRT(dvalo(i)+dvalv(i))
            END DO
            ! scale intermediate expansion
            CALL cp_fm_to_fm_submat(fm2, fma, homo, nmao, 1, 1, 1, 1)
            CALL cp_fm_to_fm_submat(fm5, fma, nvirt, nmao, 1, 1, homo+1, 1)
            CALL cp_fm_column_scale(fma, dval)
            ! Orthogonalization
            CALL cp_fm_create(fmwork, fm_struct_e)
            CALL cp_gemm("T", "N", nmao, nmao, nmo, 1.0_dp, fma, fma, 0.0_dp, fm6)
            IF (full_ortho) THEN
               ! full orthogonalization
               CALL cp_fm_power(fm6, fmwork, -0.5_dp, 1.0e-12_dp, ndep)
               IF (ndep > 0 .AND. unit_nr > 0) THEN
                  WRITE (unit_nr, '(T2,A,T71,I10)') 'Warning: linear dependent basis   ', ndep
               END IF
               CALL cp_gemm("N", "N", nmo, nmao, nmao, 1.0_dp, fma, fm6, 0.0_dp, fmb)
            ELSE
               ! orthogonalize on-atom blocks
               CALL copy_fm_to_dbcsr(fm6, sortho, keep_sparsity=.TRUE.)
               CALL diag_sqrt_invert(sortho)
               CALL copy_dbcsr_to_fm(sortho, fm6)
               CALL cp_gemm("N", "N", nmo, nmao, nmao, 1.0_dp, fma, fm6, 0.0_dp, fmb)
            END IF
            ! store as QUAMBO
            CALL cp_gemm("N", "N", nao, nmao, nmo, 1.0_dp, fm_mos, fmb, 0.0_dp, fm1)
            CALL copy_fm_to_dbcsr(fm1, quambo(ispin)%matrix)
            CALL cp_dbcsr_filter(quambo(ispin)%matrix, eps_filter)
            !
            DEALLOCATE (eigval, dval, dvalo, dvalv)
            CALL cp_fm_release(fm3)
            CALL cp_fm_release(fm4)
            CALL cp_fm_release(fm5)
            CALL cp_fm_release(fm6)
            CALL cp_fm_release(fmwork)
            CALL cp_fm_struct_release(fm_struct_c)
            CALL cp_fm_struct_release(fm_struct_d)
            CALL cp_fm_struct_release(fm_struct_e)
         END DO

         ! clean up
         CALL cp_fm_release(fm1)
         CALL cp_fm_release(fm2)
         CALL cp_fm_release(fma)
         CALL cp_fm_release(fmb)
         CALL cp_fm_struct_release(fm_struct_a)
         CALL cp_fm_struct_release(fm_struct_b)

         ! Start wfn analysis
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(/,T2,A)') 'Localized Minimal Basis Wavefunction Analysis'
         END IF

         ! localization of basis
         DO ispin = 1, nspin
            amax = cp_dbcsr_get_occupation(quambo(ispin)%matrix)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(/,T2,A,I2,T69,F10.3,A2,/)') &
                  'Occupation of Basis Function Representation (Spin) ', ispin, amax*100._dp, ' %'
            END IF
         END DO

         CALL cp_fm_struct_create(fm_struct_a, nrow_global=nao, ncol_global=nmao, &
                                  para_env=para_env, context=blacs_env)
         CALL cp_fm_create(fm1, fm_struct_a)
         CALL cp_fm_struct_create(fm_struct_b, nrow_global=nmao, ncol_global=nmo, &
                                  para_env=para_env, context=blacs_env)
         CALL cp_fm_create(fm2, fm_struct_b)
         CALL cp_fm_create(fm3, fm_struct_b)
         CALL cp_fm_struct_create(fm_struct_c, nrow_global=nmo, ncol_global=nmo, &
                                  para_env=para_env, context=blacs_env)
         CALL cp_fm_create(fm4, fm_struct_c)
         ALLOCATE (fnorm(nmo, nspin), ecount(natom, 3, nspin), prmao(natom, nspin))
         ecount = 0
         prmao = 0.0_dp
         DO ispin = 1, nspin
            CALL cp_dbcsr_init(smao)
            CALL cp_dbcsr_create(smao, name="S*QM", template=mao_coef(1)%matrix)
            smat => matrix_s(1, 1)%matrix
            CALL cp_dbcsr_multiply("N", "N", 1.0_dp, smat, quambo(ispin)%matrix, 0.0_dp, smao)
            ! calculate atomic extend of basis
            CALL pm_extend(quambo(ispin)%matrix, smao, ecount(:, :, ispin))
            CALL cp_dbcsr_init(sinv)
            CALL cp_dbcsr_create(sinv, name="QM*S*QM", template=sqmat(ispin)%matrix)
            CALL cp_dbcsr_multiply("T", "N", 1.0_dp, quambo(ispin)%matrix, smao, 0.0_dp, sqmat(ispin)%matrix)
            ! atomic MAO projection
            CALL project_mao(mao_coef(ispin)%matrix, smao, sqmat(ispin)%matrix, prmao(:, ispin))
            ! invert overlap
            CALL invert_Hotelling(sinv, sqmat(ispin)%matrix, 1.e-6_dp, silent=.TRUE.)
            CALL cp_dbcsr_init(smaox)
            CALL cp_dbcsr_create(smaox, name="S*QM*SINV", template=smao)
            CALL cp_dbcsr_multiply("N", "N", 1.0_dp, smao, sinv, 0.0_dp, smaox)
            CALL copy_dbcsr_to_fm(smaox, fm1)
            CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=fm_mos, homo=homo)
            CALL cp_gemm("T", "N", nmao, nmo, nao, 1.0_dp, fm1, fm_mos, 0.0_dp, fm2)
            CALL cp_dbcsr_sm_fm_multiply(sqmat(ispin)%matrix, fm2, fm3, nmo)
            CALL cp_gemm("T", "N", nmo, nmo, nmao, 1.0_dp, fm2, fm3, 0.0_dp, fm4)
            CALL cp_fm_get_diag(fm4, fnorm(1:nmo, ispin))
            ! fm2 are the projected MOs (in MAO basis); orthogonalize the occupied subspace
            CALL make_basis_lowdin(vmatrix=fm2, ncol=homo, matrix_s=sqmat(ispin)%matrix)
            ! pmat
            CALL get_mo_set(mos(ispin)%mo_set, occupation_numbers=occupation_numbers, maxocc=alpha)
            occeq = ALL(occupation_numbers(1:homo) == alpha)
            CALL cp_dbcsr_copy(pqmat(ispin)%matrix, sqmat(ispin)%matrix)
            CALL cp_dbcsr_set(pqmat(ispin)%matrix, 0.0_dp)
            IF (occeq) THEN
               CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=pqmat(ispin)%matrix, matrix_v=fm2, &
                                          ncol=homo, alpha=alpha, keep_sparsity=.FALSE.)
            ELSE
               CALL cp_fm_to_fm(fm2, fm3)
               CALL cp_fm_column_scale(fm3, occupation_numbers(1:homo))
               alpha = 1.0_dp
               CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=pqmat(ispin)%matrix, matrix_v=fm2, &
                                          matrix_g=fm3, ncol=homo, alpha=alpha, keep_sparsity=.TRUE.)
            END IF

            CALL cp_dbcsr_release(smao)
            CALL cp_dbcsr_release(smaox)
            CALL cp_dbcsr_release(sinv)
         END DO
         ! Basis extension
         CALL mp_sum(ecount, para_env%group)
         IF (unit_nr > 0) THEN
            IF (nspin == 1) THEN
               WRITE (unit_nr, '(T2,A,T20,A,T40,A,T60,A)') 'Ref. Atom', ' # > 0.100 ', ' # > 0.010 ', ' # > 0.001 '
               DO i = 1, natom
                  WRITE (unit_nr, '(T2,I8,T20,I10,T40,I10,T60,I10)') i, ecount(i, 1:3, 1)
               END DO
            ELSE
               WRITE (unit_nr, '(T2,A,T20,A,T40,A,T60,A)') 'Ref. Atom', ' # > 0.100 ', ' # > 0.010 ', ' # > 0.001 '
               DO i = 1, natom
                  WRITE (unit_nr, '(T2,I8,T20,2I6,T40,2I6,T60,2I6)') &
                     i, ecount(i, 1, 1:2), ecount(i, 2, 1:2), ecount(i, 3, 1:2)
               END DO
            END IF
         END IF
         ! MAO projection
         CALL mp_sum(prmao, para_env%group)
         IF (unit_nr > 0) THEN
            DO ispin = 1, nspin
               WRITE (unit_nr, '(/,T2,A,I2)') 'Projection on same atom MAO orbitals: Spin ', ispin
               DO i = 1, natom, 2
                  IF (i < natom) THEN
                     WRITE (unit_nr, '(T2,A,I8,T20,A,F10.6,T42,A,I8,T60,A,F10.6)') &
                        " Atom:", i, "Projection:", prmao(i, ispin), " Atom:", i+1, "Projection:", prmao(i+1, ispin)
                  ELSE
                     WRITE (unit_nr, '(T2,A,I8,T20,A,F10.6)') " Atom:", i, "Projection:", prmao(i, ispin)
                  END IF
               END DO
            END DO
         END IF
         ! MO expansion completness
         DO ispin = 1, nspin
            CALL get_mo_set(mos(ispin)%mo_set, homo=homo, nmo=nmo)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(/,T2,A,I2)') 'MO expansion in Localized Minimal Basis: Spin ', ispin
               WRITE (unit_nr, '(T64,A)') 'Occupied Orbitals'
               WRITE (unit_nr, '(8F10.6)') fnorm(1:homo, ispin)
               WRITE (unit_nr, '(T65,A)') 'Virtual Orbitals'
               WRITE (unit_nr, '(8F10.6)') fnorm(homo+1:nmo, ispin)
            END IF
         END DO
         ! Mulliken population
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(/,T2,A)') 'Mulliken Population Analysis '
         END IF
         ALLOCATE (mcharge(natom, nspin))
         DO ispin = 1, nspin
            CALL cp_dbcsr_trace(pqmat(ispin)%matrix, sqmat(ispin)%matrix, trace)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, '(T2,A,I2,T66,F15.4)') 'Number of Electrons: Trace(PS) Spin ', ispin, trace
            END IF
            CALL mulliken_charges(pqmat(ispin)%matrix, sqmat(ispin)%matrix, para_env, mcharge(:, ispin))
         END DO
         CALL print_atomic_charges(particle_set, qs_kind_set, unit_nr, "Minimal Basis Mulliken Charges", &
                                   electronic_charges=mcharge)
         ! Mayer bond orders
         IF (do_bondorder) THEN
            ALLOCATE (border(natom, natom))
            border = 0.0_dp
            CALL cp_dbcsr_init(psmat)
            CALL cp_dbcsr_create(psmat, name="PS", template=sqmat(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
            CALL cp_dbcsr_init(spmat)
            CALL cp_dbcsr_create(spmat, name="SP", template=sqmat(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
            filter_eps = 1.e-6_dp
            DO ispin = 1, nspin
               CALL cp_dbcsr_multiply("N", "N", 1.0_dp, pqmat(ispin)%matrix, sqmat(ispin)%matrix, 0.0_dp, psmat, &
                                      filter_eps=filter_eps)
               CALL cp_dbcsr_multiply("N", "N", 1.0_dp, sqmat(ispin)%matrix, pqmat(ispin)%matrix, 0.0_dp, spmat, &
                                      filter_eps=filter_eps)
               CALL compute_bond_order(psmat, spmat, border)
            END DO
            CALL mp_sum(border, para_env%group)
            border = border*REAL(nspin, KIND=dp)
            CALL cp_dbcsr_release(psmat)
            CALL cp_dbcsr_release(spmat)
            CALL print_bond_orders(particle_set, unit_nr, border)
            DEALLOCATE (border)
         END IF

         ! for printing purposes we now copy the QUAMBOs into MO format
         ALLOCATE (mbas(nspin))
         DO ispin = 1, nspin
            CALL allocate_mo_set(mbas(ispin)%mo_set, nao, nmao, nmao, 0.0_dp, 1.0_dp, 0.0_dp)
            mo_set => mbas(ispin)%mo_set
            CALL set_mo_set(mo_set, homo=nmao)
            ALLOCATE (mo_set%eigenvalues(nmao))
            mo_set%eigenvalues = 0.0_dp
            ALLOCATE (mo_set%occupation_numbers(nmao))
            mo_set%occupation_numbers = 1.0_dp
            CALL cp_fm_create(mo_set%mo_coeff, fm_struct_a)
            CALL copy_dbcsr_to_fm(quambo(ispin)%matrix, mo_set%mo_coeff)
         END DO

         ! Print basis functions: cube files
         DO ispin = 1, nspin
            mo_set => mbas(ispin)%mo_set
            CALL get_mo_set(mo_set, mo_coeff=fm_mos)
            CALL post_minbas_cubes(qs_env, input_section, fm_mos, ispin)
         END DO
         ! Print basis functions: molden format
         CALL write_mos_molden(mbas, qs_kind_set, particle_set, input_section)
         DO ispin = 1, nspin
            CALL deallocate_mo_set(mbas(ispin)%mo_set)
         END DO
         DEALLOCATE (mbas)

         DEALLOCATE (fnorm, ecount, prmao, mcharge)
         CALL cp_fm_release(fm1)
         CALL cp_fm_release(fm2)
         CALL cp_fm_release(fm3)
         CALL cp_fm_release(fm4)
         CALL cp_fm_struct_release(fm_struct_a)
         CALL cp_fm_struct_release(fm_struct_b)
         CALL cp_fm_struct_release(fm_struct_c)

         ! clean up
         CALL cp_dbcsr_release(sortho)
         CALL cp_dbcsr_deallocate_matrix_set(mao_coef)
         CALL cp_dbcsr_deallocate_matrix_set(quambo)
         CALL cp_dbcsr_deallocate_matrix_set(sqmat)
         CALL cp_dbcsr_deallocate_matrix_set(pqmat)

      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(/,T2,A)') &
            '!--------------------------END OF MINBAS ANALYSIS-----------------------------!'
      END IF

      CALL timestop(handle)

   END SUBROUTINE minbas_analysis

! **************************************************************************************************
!> \brief ...
!> \param sortho ...
! **************************************************************************************************
   SUBROUTINE diag_sqrt_invert(sortho)
      TYPE(cp_dbcsr_type)                                :: sortho

      CHARACTER(len=*), PARAMETER :: routineN = 'diag_sqrt_invert', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iatom, info, jatom, lwork, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: w, work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: amat, bmat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sblock
      TYPE(cp_dbcsr_iterator)                            :: dbcsr_iter

      CALL cp_dbcsr_iterator_start(dbcsr_iter, sortho)
      DO WHILE (cp_dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL cp_dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, sblock)
         CPASSERT(iatom == jatom)
         n = SIZE(sblock, 1)
         lwork = MAX(n*n, 100)
         ALLOCATE (amat(n, n), bmat(n, n), w(n), work(lwork))
         amat(1:n, 1:n) = sblock(1:n, 1:n)
         info = 0
         CALL lapack_ssyev("V", "U", n, amat, n, w, work, lwork, info)
         CPASSERT(info == 0)
         w(1:n) = 1._dp/SQRT(w(1:n))
         DO i = 1, n
            bmat(1:n, i) = amat(1:n, i)*w(i)
         END DO
         sblock(1:n, 1:n) = MATMUL(amat, TRANSPOSE(bmat))
         DEALLOCATE (amat, bmat, w, work)
      END DO
      CALL cp_dbcsr_iterator_stop(dbcsr_iter)

   END SUBROUTINE diag_sqrt_invert

! **************************************************************************************************
!> \brief ...
!> \param quambo ...
!> \param smao ...
!> \param ecount ...
! **************************************************************************************************
   SUBROUTINE pm_extend(quambo, smao, ecount)
      TYPE(cp_dbcsr_type)                                :: quambo, smao
      INTEGER, DIMENSION(:, :), INTENT(INOUT)            :: ecount

      CHARACTER(len=*), PARAMETER :: routineN = 'pm_extend', routineP = moduleN//':'//routineN

      INTEGER                                            :: iatom, jatom, n
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: wij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: qblock, sblock
      TYPE(cp_dbcsr_iterator)                            :: dbcsr_iter

      CALL cp_dbcsr_iterator_start(dbcsr_iter, quambo)
      DO WHILE (cp_dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL cp_dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, qblock)
         CALL cp_dbcsr_get_block_p(matrix=smao, row=iatom, col=jatom, BLOCK=sblock, found=found)
         IF (found) THEN
            n = SIZE(qblock, 2)
            wij = ABS(SUM(qblock*sblock))/REAL(n, KIND=dp)
            IF (wij > 0.1_dp) THEN
               ecount(jatom, 1) = ecount(jatom, 1)+1
            ELSEIF (wij > 0.01_dp) THEN
               ecount(jatom, 2) = ecount(jatom, 2)+1
            ELSEIF (wij > 0.001_dp) THEN
               ecount(jatom, 3) = ecount(jatom, 3)+1
            END IF
         END IF
      END DO
      CALL cp_dbcsr_iterator_stop(dbcsr_iter)

   END SUBROUTINE pm_extend

! **************************************************************************************************
!> \brief ...
!> \param mao ...
!> \param smao ...
!> \param sovl ...
!> \param prmao ...
! **************************************************************************************************
   SUBROUTINE project_mao(mao, smao, sovl, prmao)
      TYPE(cp_dbcsr_type)                                :: mao, smao, sovl
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: prmao

      CHARACTER(len=*), PARAMETER :: routineN = 'project_mao', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iatom, jatom, n
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: wi
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: qblock, sblock, so
      TYPE(cp_dbcsr_iterator)                            :: dbcsr_iter

      CALL cp_dbcsr_iterator_start(dbcsr_iter, mao)
      DO WHILE (cp_dbcsr_iterator_blocks_left(dbcsr_iter))
         CALL cp_dbcsr_iterator_next_block(dbcsr_iter, iatom, jatom, qblock)
         CPASSERT(iatom == jatom)
         CALL cp_dbcsr_get_block_p(matrix=smao, row=iatom, col=jatom, BLOCK=sblock, found=found)
         IF (found) THEN
            CALL cp_dbcsr_get_block_p(matrix=sovl, row=iatom, col=jatom, BLOCK=so, found=found)
            n = SIZE(qblock, 2)
            DO i = 1, n
               wi = SUM(qblock(:, i)*sblock(:, i))
               prmao(iatom) = prmao(iatom)+wi/so(i, i)
            END DO
         END IF
      END DO
      CALL cp_dbcsr_iterator_stop(dbcsr_iter)

   END SUBROUTINE project_mao

! **************************************************************************************************
!> \brief Computes and prints the Cube Files for the minimal basis set
!> \param qs_env ...
!> \param print_section ...
!> \param minbas_coeff ...
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE post_minbas_cubes(qs_env, print_section, minbas_coeff, ispin)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: print_section
      TYPE(cp_fm_type), POINTER                          :: minbas_coeff
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(len=*), PARAMETER :: routineN = 'post_minbas_cubes', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename, title
      INTEGER                                            :: i, i_rep, ivec, iw, j, n_rep, natom
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes, first_bas, ilist, stride
      LOGICAL                                            :: explicit
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: wf_g, wf_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: minbas_section

      minbas_section => section_vals_get_subs_vals(print_section, "MINBAS_CUBE")
      CALL section_vals_get(minbas_section, explicit=explicit)
      IF (.NOT. explicit) RETURN

      CPASSERT(ASSOCIATED(minbas_coeff))

      logger => cp_get_default_logger()
      stride => section_get_ivals(print_section, "MINBAS_CUBE%STRIDE")

      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      subsys=subsys, cell=cell, particle_set=particle_set, pw_env=pw_env, dft_control=dft_control)
      CALL qs_subsys_get(subsys, particles=particles)

      CALL get_qs_env(qs_env=qs_env, natom=natom)
      ALLOCATE (blk_sizes(natom), first_bas(0:natom))
      CALL get_particle_set(particle_set, qs_kind_set, nmao=blk_sizes)
      first_bas(0) = 0
      DO i = 1, natom
         first_bas(i) = first_bas(i-1)+blk_sizes(i)
      END DO

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_r%pw, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_g%pw, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      ! loop over list of atoms
      CALL section_vals_val_get(minbas_section, "ATOM_LIST", n_rep_val=n_rep)
      IF (n_rep == 0) THEN
         DO i = 1, natom
            DO ivec = first_bas(i-1)+1, first_bas(i)
               WRITE (filename, '(a4,I5.5,a1,I1.1)') "MINBAS_", ivec, "_", ispin
               WRITE (title, *) "MINIMAL BASIS ", ivec, " atom ", i, " spin ", ispin
               iw = cp_print_key_unit_nr(logger, print_section, "MINBAS_CUBE", extension=".cube", &
                                         middle_name=TRIM(filename), file_position="REWIND", log_filename=.FALSE.)
               CALL calculate_wavefunction(minbas_coeff, ivec, wf_r, wf_g, atomic_kind_set, qs_kind_set, &
                                           cell, dft_control, particle_set, pw_env)
               CALL cp_pw_to_cube(wf_r%pw, iw, title, particles=particles, stride=stride)
               CALL cp_print_key_finished_output(iw, logger, print_section, "MINBAS_CUBE")
            END DO
         END DO
      ELSE
         DO i_rep = 1, n_rep
            CALL section_vals_val_get(minbas_section, "ATOM_LIST", i_rep_val=i_rep, i_vals=ilist)
            DO i = 1, SIZE(ilist, 1)
               j = ilist(i)
               DO ivec = first_bas(j-1)+1, first_bas(j)
                  WRITE (filename, '(a4,I5.5,a1,I1.1)') "MINBAS_", ivec, "_", ispin
                  WRITE (title, *) "MINIMAL BASIS ", ivec, " atom ", j, " spin ", ispin
                  iw = cp_print_key_unit_nr(logger, print_section, "MINBAS_CUBE", extension=".cube", &
                                            middle_name=TRIM(filename), file_position="REWIND", log_filename=.FALSE.)
                  CALL calculate_wavefunction(minbas_coeff, ivec, wf_r, wf_g, atomic_kind_set, qs_kind_set, &
                                              cell, dft_control, particle_set, pw_env)
                  CALL cp_pw_to_cube(wf_r%pw, iw, title, particles=particles, stride=stride)
                  CALL cp_print_key_finished_output(iw, logger, print_section, "MINBAS_CUBE")
               END DO
            END DO
         END DO
      END IF
      DEALLOCATE (blk_sizes, first_bas)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_g%pw)

   END SUBROUTINE post_minbas_cubes

END MODULE minbas_wfn_analysis
