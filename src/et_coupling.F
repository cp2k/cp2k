!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief calculates the electron transfer coupling elements
!>      Wu, Van Voorhis, JCP 125, 164105 (2006)
!> \author fschiff (01.2007)
! **************************************************************************************************
MODULE et_coupling
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: cdft_control_type,&
                                              dft_control_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_invert,&
                                              cp_fm_transpose
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pool_get_el_struct
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE cube_utils,                      ONLY: cube_info_type
   USE dbcsr_api,                       ONLY: dbcsr_deallocate_matrix_set,&
                                              dbcsr_p_type
   USE hirshfeld_methods,               ONLY: create_shape_function
   USE hirshfeld_types,                 ONLY: get_hirshfeld_info,&
                                              hirshfeld_type,&
                                              set_hirshfeld_info
   USE input_constants,                 ONLY: &
        becke_cutoff_element, becke_cutoff_global, cdft_combined_acceptor, cdft_combined_all, &
        cdft_combined_constraint, cdft_combined_donor, cdft_density_constraint, &
        cdft_magnetization_constraint, ref_charge_atomic, ref_charge_mulliken
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kahan_sum,                       ONLY: accurate_sum,&
                                              accurate_sum_arrays_product
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_sum
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_set
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_p_type,&
                                              pw_release
   USE qs_collocate_density,            ONLY: collocate_pgf_product_rspace
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integrate_potential,          ONLY: integrate_pgf_product_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_matrix_pools,                 ONLY: mpools_get
   USE qs_mo_types,                     ONLY: get_mo_set
   USE qs_modify_pab_block,             ONLY: FUNC_AB
   USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                              mpole_rho_atom,&
                                              rho0_mpole_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE realspace_grid_cube,             ONLY: cube_to_pw
   USE realspace_grid_types,            ONLY: &
        pw2rs, realspace_grid_desc_type, realspace_grid_type, rs2pw, rs_grid_create, &
        rs_grid_release, rs_grid_retain, rs_grid_zero, rs_pw_transfer
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling'
   INTEGER, PARAMETER, PRIVATE          :: divide = 1, multiply = 2
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

! *** Public subroutines ***

   PUBLIC :: calc_et_coupling, becke_restraint, hirshfeld_constraint

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE calc_et_coupling(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calc_et_coupling', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iw, j, k, nao, ncol_local, &
                                                            nmo, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: Sda, strength, Waa, Wbb, Wda
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: a, b, S_det
      REAL(KIND=dp), DIMENSION(2)                        :: eigenv
      REAL(KIND=dp), DIMENSION(2, 2)                     :: S_mat, tmp_mat, U, W_mat
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: rest_MO
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: mo_mo_fm_pools
      TYPE(cp_fm_struct_type), POINTER                   :: mo_mo_fmstruct
      TYPE(cp_fm_type), POINTER                          :: inverse_mat, SMO, Tinverse, tmp2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: et_coupling_section

      NULLIFY (tmp2, SMO, rest_MO, mo_mo_fmstruct, rest_MO, energy, matrix_s, dft_control, para_env)

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      et_coupling_section => section_vals_get_subs_vals(qs_env%input, &
                                                        "PROPERTIES%ET_COUPLING")

      CALL get_qs_env(qs_env, dft_control=dft_control, para_env=para_env)

      iw = cp_print_key_unit_nr(logger, et_coupling_section, "PROGRAM_RUN_INFO", &
                                extension=".log")

      ALLOCATE (rest_MO(2))
      ALLOCATE (a(dft_control%nspins))
      ALLOCATE (b(dft_control%nspins))
      ALLOCATE (S_det(dft_control%nspins))

      CALL mpools_get(qs_env%mpools, mo_mo_fm_pools=mo_mo_fm_pools)
      mo_mo_fmstruct => fm_pool_get_el_struct(mo_mo_fm_pools(1)%pool)
      DO i = 1, dft_control%nspins
         mo_mo_fmstruct => fm_pool_get_el_struct(mo_mo_fm_pools(i)%pool)

         CALL get_mo_set(mo_set=qs_env%mos(i)%mo_set, &
                         nao=nao, &
                         nmo=nmo)

         CALL cp_fm_create(matrix=tmp2, &
                           matrix_struct=qs_env%mos(i)%mo_set%mo_coeff%matrix_struct, &
                           name="ET_TMP"//TRIM(ADJUSTL(cp_to_string(2)))//"MATRIX")
         CALL cp_fm_create(matrix=inverse_mat, &
                           matrix_struct=mo_mo_fmstruct, &
                           name="INVERSE"//TRIM(ADJUSTL(cp_to_string(2)))//"MATRIX")
         CALL cp_fm_create(matrix=Tinverse, &
                           matrix_struct=mo_mo_fmstruct, &
                           name="T_INVERSE"//TRIM(ADJUSTL(cp_to_string(2)))//"MATRIX")
         CALL cp_fm_create(matrix=SMO, &
                           matrix_struct=mo_mo_fmstruct, &
                           name="ET_SMO"//TRIM(ADJUSTL(cp_to_string(1)))//"MATRIX")
         DO j = 1, 2
            NULLIFY (rest_MO(j)%matrix)
            CALL cp_fm_create(matrix=rest_MO(j)%matrix, &
                              matrix_struct=mo_mo_fmstruct, &
                              name="ET_rest_MO"//TRIM(ADJUSTL(cp_to_string(j)))//"MATRIX")
         END DO

!   calculate MO-overlap

         CALL get_qs_env(qs_env, matrix_s=matrix_s)
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, qs_env%et_coupling%et_mo_coeff(i)%matrix, &
                                      tmp2, nmo, 1.0_dp, 0.0_dp)
         CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, &
                      qs_env%mos(i)%mo_set%mo_coeff, &
                      tmp2, 0.0_dp, SMO)

!    calculate the MO-representation of the restraint matrix A

         CALL cp_dbcsr_sm_fm_multiply(qs_env%et_coupling%rest_mat(1)%matrix, &
                                      qs_env%et_coupling%et_mo_coeff(i)%matrix, &
                                      tmp2, nmo, 1.0_dp, 0.0_dp)

         CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, &
                      qs_env%mos(i)%mo_set%mo_coeff, &
                      tmp2, 0.0_dp, rest_MO(1)%matrix)

!    calculate the MO-representation of the restraint matrix D

         CALL cp_dbcsr_sm_fm_multiply(qs_env%et_coupling%rest_mat(2)%matrix, &
                                      qs_env%mos(i)%mo_set%mo_coeff, &
                                      tmp2, nmo, 1.0_dp, 0.0_dp)

         CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, &
                      qs_env%et_coupling%et_mo_coeff(i)%matrix, &
                      tmp2, 0.0_dp, rest_MO(2)%matrix)

         CALL cp_fm_invert(SMO, inverse_mat, S_det(i))

         CALL cp_fm_get_info(inverse_mat, nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices)
         b(i) = 0.0_dp

         DO j = 1, ncol_local
            DO k = 1, nrow_local
               b(i) = b(i)+rest_MO(2)%matrix%local_data(k, j)*inverse_mat%local_data(k, j)
            END DO
         END DO

         CALL cp_fm_transpose(inverse_mat, Tinverse)
         a(i) = 0.0_dp
         DO j = 1, ncol_local
            DO k = 1, nrow_local
               a(i) = a(i)+rest_MO(1)%matrix%local_data(k, j)*Tinverse%local_data(k, j)
            END DO
         END DO
         IF (dft_control%qs_control%becke_control%constraint_type == &
             cdft_magnetization_constraint .AND. i == 2) THEN
            a(i) = -a(i)
            b(i) = -b(i)
         END IF
         CALL mp_sum(a(i), para_env%group)

         CALL mp_sum(b(i), para_env%group)

         CALL cp_fm_release(tmp2)
         CALL cp_fm_release(rest_MO(1)%matrix)
         CALL cp_fm_release(rest_MO(2)%matrix)
         CALL cp_fm_release(SMO)
         CALL cp_fm_release(Tinverse)
         CALL cp_fm_release(inverse_mat)
      END DO
      DEALLOCATE (rest_MO)

!    solve eigenstates for the projector matrix

      IF (dft_control%nspins == 2) THEN
         Sda = S_det(1)*S_det(2)
         Wda = ((a(1)+a(2))+(b(1)+b(2)))*0.5_dp*Sda
      ELSE
         Sda = S_det(1)**2
         Wda = (a(1)+b(1))*Sda
      END IF

      IF (dft_control%qs_control%ddapc_restraint) THEN
         Waa = qs_env%et_coupling%order_p
         Wbb = dft_control%qs_control%ddapc_restraint_control(1)%ddapc_restraint_control%ddapc_order_p
         strength = dft_control%qs_control%ddapc_restraint_control(1)%ddapc_restraint_control%strength
      END IF
      IF (dft_control%qs_control%becke_restraint) THEN
         Waa = qs_env%et_coupling%order_p
         Wbb = dft_control%qs_control%becke_control%becke_order_p(1)
         strength = dft_control%qs_control%becke_control%strength(1)
      END IF

!!   construct S and W   !!!
      S_mat(1, 1) = 1.0_dp
      S_mat(2, 2) = 1.0_dp
      S_mat(2, 1) = Sda
      S_mat(1, 2) = Sda

      W_mat(1, 1) = Wbb
      W_mat(2, 2) = Waa
      W_mat(2, 1) = Wda
      W_mat(1, 2) = Wda

!!  solve WC=SCN
      CALL diamat_all(S_mat, eigenv, .TRUE.)
      ! U = S**(-1/2)
      U = 0.0_dp
      U(1, 1) = 1.0_dp/SQRT(eigenv(1))
      U(2, 2) = 1.0_dp/SQRT(eigenv(2))
      tmp_mat = MATMUL(U, TRANSPOSE(S_mat))
      U = MATMUL(S_mat, tmp_mat)
      tmp_mat = MATMUL(W_mat, U)
      W_mat = MATMUL(U, tmp_mat)
      CALL diamat_all(W_mat, eigenv, .TRUE.)
      tmp_mat = MATMUL(U, W_mat)

      CALL get_qs_env(qs_env, energy=energy)
      W_mat(1, 1) = energy%total
      W_mat(2, 2) = qs_env%et_coupling%energy
      a(1) = (energy%total+strength*Wbb)*Sda-strength*Wda
      a(2) = (qs_env%et_coupling%energy+qs_env%et_coupling%e1*Waa)*Sda-qs_env%et_coupling%e1*Wda
      W_mat(1, 2) = (a(1)+a(2))*0.5_dp
      W_mat(2, 1) = W_mat(1, 2)

      S_mat = MATMUL(W_mat, (tmp_mat))
      W_mat = MATMUL(TRANSPOSE(tmp_mat), S_mat)

      IF (iw > 0) THEN
         WRITE (iw, *)
         WRITE (iw, '(T3,A,T60,(3X,F12.6))') 'Strength of constraint A          :', qs_env%et_coupling%e1
         WRITE (iw, '(T3,A,T60,(3X,F12.6))') 'Strength of constraint B          :', strength
         WRITE (iw, '(T3,A,T60,(3X,F12.6))') 'Final target value of constraint A:', Waa
         WRITE (iw, '(T3,A,T60,(3X,F12.6))') 'Final target value of constraint B:', Wbb
         WRITE (iw, *)
         WRITE (iw, '(T3,A,T60,(3X,F12.6))') &
            'Diabatic electronic coupling matrix element(mHartree):', ABS(W_mat(1, 2)*1000.0_dp)

      ENDIF

      CALL dbcsr_deallocate_matrix_set(qs_env%et_coupling%rest_mat)

      CALL cp_print_key_finished_output(iw, logger, et_coupling_section, &
                                        "PROGRAM_RUN_INFO")
      CALL timestop(handle)
   END SUBROUTINE calc_et_coupling

! **************************************************************************************************
!> \brief calculates a Becke contraint
!> \param qs_env the qs_env where to build the constraint
!> \param becke_const the Becke real space potential
!> \param charge atomic Becke real space potentials defined on the constraint atoms
!> \param calc_pot if the potential needs to be recalculated or just integrated
!> \param calculate_forces logical if potential has to be calculated or only_energy
!> \par   History
!>        Created 01.2007 [fschiff]
!>        Extended functionality 12/15-12/16 [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE becke_restraint(qs_env, becke_const, charge, calc_pot, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type)                                    :: becke_const
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: charge
      LOGICAL                                            :: calc_pot, calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'becke_restraint', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=2)                                   :: element_symbol
      INTEGER :: atom_a, bounds(2), dir, handle, i, iatom, iex, ikind, ind(3), ip, ithread, ivar, &
         iw, j, jatom, k, katom, lb_index, natom, nelectron_total, ngrad, nkind, np(3), npme, &
         nskipped, nthread, numexp, nvar, tmp_index(2), ub_index, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: catom
      INTEGER, DIMENSION(2, 3)                           :: bo, bo_conf
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores, stride
      LOGICAL                                            :: build, in_memory, paw_atom, store_vectors
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: is_combined, is_constraint, skip_me
      REAL(kind=dp) :: alpha, chi, coef, dist1, dist2, dmyexp, dvol, dynamic_radius, eps_cavity, &
         ircov, jrcov, lb, my1, my1_homo, myexp, nelectron_frag, sign, sum_cell_f_all, &
         sum_cell_f_combined, sum_cell_f_constr, tc, th, tmp_const, ub, uij, zeff
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: cell_functions, coefficients, dE, &
                                                            distances, ds_dR_i, ds_dR_j, &
                                                            gapw_offset, strength, target_val
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: d_sum_const_dR, d_sum_Pm_dR, &
                                                            distance_vecs, dP_i_dRi, &
                                                            electronic_charge, position_vecs, R12
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dP_i_dRj, pair_dist_vecs
      REAL(kind=dp), DIMENSION(3)                        :: cell_v, dist_vec, dmy_dR_i, dmy_dR_j, &
                                                            dr, dr1_r2, dr_i_dR, dr_ij_dR, &
                                                            dr_j_dR, grid_p, r, r1, ra, shift
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cutoffs
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: gradients
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hirshfeld_type), POINTER                      :: cavity_env
      TYPE(mpole_rho_atom), DIMENSION(:), POINTER        :: mp_rho
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_frag
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(realspace_grid_type), POINTER                 :: rs_cavity
      TYPE(rho0_mpole_type), POINTER                     :: rho0_mpole
      TYPE(section_vals_type), POINTER                   :: becke_restraint_section

      store_vectors = .TRUE. ! Reuse vectors as much as possible
      NULLIFY (cores, stride, atom_list, cutoffs, cell, para_env, &
               dft_control, particle_set, rho_r, energy, rho, logger, &
               becke_restraint_section, qs_kind_set, gradients, particles, &
               subsys, pab, pw_env, rs_cavity, cavity_env, auxbas_pw_pool, &
               atomic_kind_set, mp_rho, rho0_mpole)
      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      particle_set=particle_set, &
                      rho=rho, &
                      natom=natom, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      qs_kind_set=qs_kind_set)

      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(qs_kind_set))

      becke_restraint_section => section_vals_get_subs_vals(qs_env%input, &
                                                            "DFT%QS%BECKE_RESTRAINT")
      iw = cp_print_key_unit_nr(logger, becke_restraint_section, "PROGRAM_RUN_INFO", &
                                extension=".cdftLog")

      CALL qs_rho_get(rho, rho_r=rho_r)

      IF (dft_control%qs_control%becke_restraint) THEN
         SELECT CASE (dft_control%qs_control%becke_control%constraint_type)
         CASE (cdft_density_constraint)
            nvar = 1
            ngrad = 1
         CASE (cdft_magnetization_constraint)
            nvar = 1
            ngrad = 1
            IF (dft_control%nspins == 1) &
               CALL cp_abort(__LOCATION__, &
                             "Becke magnetization density constraint requires UKS calculation.")
         CASE (cdft_combined_constraint)
            nvar = 2
            IF (dft_control%nspins == 1) &
               CALL cp_abort(__LOCATION__, &
                             "Becke combined charge+spin constraint requires UKS calculation.")
            SELECT CASE (dft_control%qs_control%becke_control%combined_type)
            CASE (cdft_combined_all)
               ngrad = 1
            CASE (cdft_combined_acceptor, cdft_combined_donor)
               ngrad = 2
            END SELECT
         END SELECT
         ALLOCATE (strength(nvar))
         ALLOCATE (target_val(nvar))
         ALLOCATE (dE(nvar))
         strength(:) = dft_control%qs_control%becke_control%strength(:)
         target_val(:) = dft_control%qs_control%becke_control%target(:)
         dE = 0.0_dp
         dvol = becke_const%pw%pw_grid%dvol
         eps_cavity = dft_control%qs_control%becke_control%eps_cavity

         IF (calc_pot) THEN

            becke_const%pw%cr3d = 0.0_dp

            ! Setup atomic radii for adjusting cell boundaries
            IF (dft_control%qs_control%becke_control%adjust) THEN
               IF (.NOT. ASSOCIATED(dft_control%qs_control%becke_control%radii)) THEN
                  CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
                  IF (.NOT. SIZE(atomic_kind_set) == SIZE(dft_control%qs_control%becke_control%radii_tmp)) &
                     CALL cp_abort(__LOCATION__, &
                                   "Length of keyword BECKE_RESTRAINT\ATOMIC_RADII does not "// &
                                   "match number of atomic kinds in the input coordinate file.")
                  ALLOCATE (dft_control%qs_control%becke_control%radii(SIZE(atomic_kind_set)))
                  dft_control%qs_control%becke_control%radii(:) = dft_control%qs_control%becke_control%radii_tmp(:)
                  DEALLOCATE (dft_control%qs_control%becke_control%radii_tmp)
               END IF
            END IF

            ! Setup cutoff scheme
            IF (.NOT. ASSOCIATED(dft_control%qs_control%becke_control%cutoffs)) THEN
               CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
               ALLOCATE (dft_control%qs_control%becke_control%cutoffs(natom))
               SELECT CASE (dft_control%qs_control%becke_control%cutoff_type)
               CASE (becke_cutoff_global)
                  dft_control%qs_control%becke_control%cutoffs(:) = dft_control%qs_control%becke_control%rglobal
               CASE (becke_cutoff_element)
                  IF (.NOT. SIZE(atomic_kind_set) == SIZE(dft_control%qs_control%becke_control%cutoffs_tmp)) &
                     CALL cp_abort(__LOCATION__, &
                                   "Length of keyword BECKE_RESTRAINT\ELEMENT_CUTOFFS does not "// &
                                   "match number of atomic kinds in the input coordinate file.")
                  DO ikind = 1, SIZE(atomic_kind_set)
                     CALL get_atomic_kind(atomic_kind_set(ikind), natom=katom, atom_list=atom_list)
                     DO iatom = 1, katom
                        atom_a = atom_list(iatom)
                        dft_control%qs_control%becke_control%cutoffs(atom_a) = &
                           dft_control%qs_control%becke_control%cutoffs_tmp(ikind)
                     END DO
                  END DO
                  DEALLOCATE (dft_control%qs_control%becke_control%cutoffs_tmp)
               END SELECT
            END IF
            cutoffs => dft_control%qs_control%becke_control%cutoffs

            IF (dft_control%qs_control%becke_control%atomic_charges) THEN
               DO iatom = 1, dft_control%qs_control%becke_control%natoms
                  charge(iatom)%pw%cr3d = 0.0_dp
               END DO
            END IF

            build = .FALSE.
            IF (dft_control%qs_control%becke_control%adjust .AND. &
                .NOT. ASSOCIATED(dft_control%qs_control%becke_control%aij)) THEN
               ALLOCATE (dft_control%qs_control%becke_control%aij(natom, natom))
               build = .TRUE.
            END IF
            ALLOCATE (catom(dft_control%qs_control%becke_control%natoms))

            IF (dft_control%qs_control%becke_control%dynamic_confine .OR. &
                dft_control%qs_control%becke_control%save_pot .OR. &
                dft_control%qs_control%becke_control%cavity_confine .OR. &
                dft_control%qs_control%becke_control%should_skip .OR. &
                dft_control%qs_control%becke_control%first_iteration) THEN
               ALLOCATE (is_constraint(natom))
               is_constraint = .FALSE.
            END IF
            ! This boolean is needed to prevent calculation of atom pairs ji when the pair ij has
            ! already been calculated (data for pair ji is set using symmetry)
            ! With gradient precomputation, symmetry exploited for both weight function and gradients
            ALLOCATE (skip_me(natom))
            in_memory = .FALSE.
            IF (dft_control%qs_control%becke_control%save_pot) THEN
               in_memory = dft_control%qs_control%becke_control%in_memory
            END IF
            IF (in_memory .OR. dft_control%qs_control%becke_control%first_iteration) &
               ALLOCATE (coefficients(natom))

            DO i = 1, dft_control%qs_control%becke_control%natoms
               catom(i) = dft_control%qs_control%becke_control%atoms(i)
               IF (dft_control%qs_control%becke_control%dynamic_confine .OR. &
                   dft_control%qs_control%becke_control%save_pot .OR. &
                   dft_control%qs_control%becke_control%cavity_confine .OR. &
                   dft_control%qs_control%becke_control%should_skip .OR. &
                   dft_control%qs_control%becke_control%first_iteration) &
                  is_constraint(catom(i)) = .TRUE.
               IF (in_memory .OR. dft_control%qs_control%becke_control%first_iteration) &
                  coefficients(catom(i)) = dft_control%qs_control%becke_control%coeff(i)
            ENDDO

            bo = becke_const%pw%pw_grid%bounds_local
            np = becke_const%pw%pw_grid%npts
            dvol = becke_const%pw%pw_grid%dvol
            dr = becke_const%pw%pw_grid%dr
            np = becke_const%pw%pw_grid%npts
            shift = -REAL(MODULO(np, 2), dp)*dr/2.0_dp

            ALLOCATE (cell_functions(natom))
            IF (store_vectors) THEN
               ALLOCATE (distances(natom))
               ALLOCATE (distance_vecs(3, natom))
               IF (in_memory) ALLOCATE (pair_dist_vecs(3, natom, natom))
               ALLOCATE (position_vecs(3, natom))
            END IF
            ! Setup static/dynamic confinement
            IF (dft_control%qs_control%becke_control%confine) THEN
               dir = dft_control%qs_control%becke_control%confine_dir
               lb = dft_control%qs_control%becke_control%confine_bounds(1)
               ub = dft_control%qs_control%becke_control%confine_bounds(2)
               ! Convert [0,L] to [-L/2, L/2]
               lb = lb-cell%hmat(dir, dir)/2._dp
               ub = ub-cell%hmat(dir, dir)/2._dp
            ELSE IF (dft_control%qs_control%becke_control%dynamic_confine) THEN
               dir = dft_control%qs_control%becke_control%confine_dir
               lb = HUGE(0.0_dp)
               ub = -HUGE(0.0_dp)
               dynamic_radius = dft_control%qs_control%becke_control%dynamic_radius
            END IF
            ALLOCATE (R12(natom, natom))
            IF (in_memory) THEN
               ALLOCATE (ds_dR_j(3))
               ALLOCATE (ds_dR_i(3))
               ALLOCATE (d_sum_Pm_dR(3, natom))
               ALLOCATE (d_sum_const_dR(3*ngrad, natom))
               ALLOCATE (dP_i_dRj(3, natom, natom))
               ALLOCATE (dP_i_dRi(3, natom))
               th = 1.0e-8_dp
            END IF
            !calculate distances from target atom (only implemented for a diatomic system)
            DO i = 1, 3
               cell_v(i) = cell%hmat(i, i)
            END DO
            DO iatom = 1, natom-1
               DO jatom = iatom+1, natom
                  r = particle_set(iatom)%r
                  r1 = particle_set(jatom)%r
                  DO i = 1, 3
                     r(i) = MODULO(r(i), cell%hmat(i, i))-cell%hmat(i, i)/2._dp
                     r1(i) = MODULO(r1(i), cell%hmat(i, i))-cell%hmat(i, i)/2._dp
                  END DO
                  dist_vec = (r-r1)-ANINT((r-r1)/cell_v)*cell_v
                  ! Store pbc corrected position and pairwise distance vectors for later reuse
                  IF (store_vectors) THEN
                     position_vecs(:, iatom) = r(:)
                     IF (iatom == 1 .AND. jatom == natom) position_vecs(:, jatom) = r1(:)
                     IF (in_memory) THEN
                        pair_dist_vecs(:, iatom, jatom) = dist_vec(:)
                        pair_dist_vecs(:, jatom, iatom) = -dist_vec(:)
                     END IF
                  END IF
                  R12(iatom, jatom) = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
                  R12(jatom, iatom) = R12(iatom, jatom)
                  IF (dft_control%qs_control%becke_control%dynamic_confine) THEN
                     IF (is_constraint(iatom) .AND. is_constraint(jatom)) THEN
                        ! ub = max(r(dir)+dynamic_radius), lb = min(r(dir)-dynamic_radius)
                        IF ((r(dir)+dynamic_radius) .GT. ub) ub = r(dir)+dynamic_radius
                        IF ((r1(dir)+dynamic_radius) .GT. ub) ub = r1(dir)+dynamic_radius
                        IF ((r(dir)-dynamic_radius) .LT. lb) lb = r(dir)-dynamic_radius
                        IF ((r1(dir)-dynamic_radius) .LT. lb) lb = r1(dir)-dynamic_radius
                     END IF
                  END IF
                  ! Set up heteronuclear cell partitioning using user defined radii
                  IF (build) THEN
                     CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                          kind_number=ikind)
                     ircov = dft_control%qs_control%becke_control%radii(ikind)
                     CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, &
                                          kind_number=ikind)
                     jrcov = dft_control%qs_control%becke_control%radii(ikind)
                     IF (ircov .NE. jrcov) THEN
                        chi = ircov/jrcov
                        uij = (chi-1.0_dp)/(chi+1.0_dp)
                        dft_control%qs_control%becke_control%aij(iatom, jatom) = uij/(uij**2-1.0_dp)
                        IF (dft_control%qs_control%becke_control%aij(iatom, jatom) .GT. 0.5_dp) THEN
                           dft_control%qs_control%becke_control%aij(iatom, jatom) = 0.5_dp
                        ELSE IF (dft_control%qs_control%becke_control%aij(iatom, jatom) .LT. -0.5_dp) THEN
                           dft_control%qs_control%becke_control%aij(iatom, jatom) = -0.5_dp
                        END IF
                     ELSE
                        dft_control%qs_control%becke_control%aij(iatom, jatom) = 0.0_dp
                     END IF
                     dft_control%qs_control%becke_control%aij(jatom, iatom) = &
                        -dft_control%qs_control%becke_control%aij(iatom, jatom) ! Note change of sign
                  END IF
               END DO
            END DO
            ! Dump some additional information about the calculation
            IF (dft_control%qs_control%becke_control%first_iteration) THEN
               IF (iw > 0) THEN
                  WRITE (iw, '(/,T3,A,T60)') &
                     '----------------------------- Becke atomic parameters -----------------------------'
                  IF (dft_control%qs_control%becke_control%adjust) THEN
                     WRITE (iw, '(/,T3,A,A)') &
                        'Atom  Element   Coefficient', '     Cutoff (angstrom)    CDFT Radius (angstrom)'
                     DO iatom = 1, natom
                        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                             element_symbol=element_symbol, &
                                             kind_number=ikind)
                        ircov = cp_unit_from_cp2k(dft_control%qs_control%becke_control%radii(ikind), &
                                                  "angstrom")
                        IF (is_constraint(iatom)) THEN
                           coef = coefficients(iatom)
                        ELSE
                           coef = 0.0_dp
                        END IF
                        WRITE (iw, "(i7,T15,A2,T23,F8.3,T39,F8.3,T61,F8.3,T75,F8.3)") &
                           iatom, element_symbol, coef, &
                           cp_unit_from_cp2k(dft_control%qs_control%becke_control%cutoffs(iatom), "angstrom"), &
                           ircov
                     END DO
                  ELSE
                     WRITE (iw, '(/,T3,A,A)') &
                        'Atom  Element   Coefficient', '     Cutoff (angstrom)'
                     DO iatom = 1, natom
                        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                             element_symbol=element_symbol)
                        IF (is_constraint(iatom)) THEN
                           coef = coefficients(iatom)
                        ELSE
                           coef = 0.0_dp
                        END IF
                        WRITE (iw, "(i7,T15,A2,T23,F8.3,T39,F8.3,T61,F8.3,T75)") &
                           iatom, element_symbol, coef, &
                           cp_unit_from_cp2k(dft_control%qs_control%becke_control%cutoffs(iatom), "angstrom")
                     END DO
                  END IF
                  WRITE (iw, '(/,T3,A)') &
                     '-----------------------------------------------------------------------------------'
               END IF
               dft_control%qs_control%becke_control%first_iteration = .FALSE.
            END IF
            ! Convert confinement bounds to RS grid indices
            IF (dft_control%qs_control%becke_control%dynamic_confine .OR. &
                dft_control%qs_control%becke_control%confine) THEN
               tmp_const = (ub-shift(dir))/dr(dir)
               IF (tmp_const .GE. 0.0_dp) THEN
                  ub_index = CEILING(tmp_const)
               ELSE
                  ub_index = FLOOR(tmp_const)
               END IF
               tmp_const = (lb-shift(dir))/dr(dir)
               IF (tmp_const .GE. 0.0_dp) THEN
                  lb_index = CEILING(tmp_const)
               ELSE
                  lb_index = FLOOR(tmp_const)
               END IF
               ! Store confinement bound indices for use in force calculation
               dft_control%qs_control%becke_control%confine_bounds_int(1) = lb_index
               dft_control%qs_control%becke_control%confine_bounds_int(2) = ub_index
            END IF
            ! Setup cavity confinement using spherical Gaussians
            IF (dft_control%qs_control%becke_control%cavity_confine) THEN
               cavity_env => dft_control%qs_control%becke_control%cavity_env
               CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, pw_env=pw_env)
               nkind = SIZE(qs_kind_set)
               ! Setup Gaussian shape function
               IF (.NOT. ASSOCIATED(cavity_env%kind_shape_fn)) &
                  CALL create_shape_function(cavity_env, qs_kind_set, atomic_kind_set, &
                                             radius=dft_control%qs_control%becke_control%rcavity, &
                                             radii_list=dft_control%qs_control%becke_control%radii)
               ! Form cavity by summing isolated Gaussian densities over constraint atoms
               NULLIFY (rs_cavity)
               CALL pw_env_get(pw_env, auxbas_rs_grid=rs_cavity, &
                               auxbas_pw_pool=auxbas_pw_pool)
               CALL rs_grid_retain(rs_cavity)
               CALL rs_grid_zero(rs_cavity)
               ALLOCATE (pab(1, 1))
               nthread = 1
               ithread = 0

               DO ikind = 1, SIZE(atomic_kind_set)
                  numexp = cavity_env%kind_shape_fn(ikind)%numexp
                  IF (numexp <= 0) CYCLE
                  CALL get_atomic_kind(atomic_kind_set(ikind), natom=katom, atom_list=atom_list)
                  ALLOCATE (cores(katom))
                  DO iex = 1, numexp
                     alpha = cavity_env%kind_shape_fn(ikind)%zet(iex)
                     coef = cavity_env%kind_shape_fn(ikind)%coef(iex)
                     npme = 0
                     cores = 0
                     DO iatom = 1, katom
                        IF (rs_cavity%desc%parallel .AND. .NOT. rs_cavity%desc%distributed) THEN
                           ! replicated realspace grid, split the atoms up between procs
                           IF (MODULO(iatom, rs_cavity%desc%group_size) == rs_cavity%desc%my_pos) THEN
                              npme = npme+1
                              cores(npme) = iatom
                           ENDIF
                        ELSE
                           npme = npme+1
                           cores(npme) = iatom
                        ENDIF
                     END DO
                     DO j = 1, npme
                        iatom = cores(j)
                        atom_a = atom_list(iatom)
                        pab(1, 1) = coef
                        IF (store_vectors) THEN
                           ra(:) = position_vecs(:, atom_a)+cell_v(:)/2._dp
                        ELSE
                           ra(:) = pbc(particle_set(atom_a)%r, cell)
                        END IF
                        IF (is_constraint(atom_a)) &
                           CALL collocate_pgf_product_rspace(0, alpha, 0, 0, 0.0_dp, 0, ra, &
                                                             (/0.0_dp, 0.0_dp, 0.0_dp/), 0.0_dp, 1.0_dp, &
                                                             pab, 0, 0, rs_cavity, cell, pw_env%cube_info(1), &
                                                             dft_control%qs_control%eps_rho_rspace, &
                                                             ga_gb_function=FUNC_AB, ithread=ithread, &
                                                             use_subpatch=.TRUE., subpatch_pattern=0_int_8, &
                                                             lmax_global=0)
                     END DO
                  END DO
                  DEALLOCATE (cores)
               END DO
               DEALLOCATE (pab)

               CALL pw_pool_create_pw(auxbas_pw_pool, dft_control%qs_control%becke_control%cavity%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL rs_pw_transfer(rs_cavity, dft_control%qs_control%becke_control%cavity%pw, rs2pw)
               CALL rs_grid_release(rs_cavity)
               ! Ignore grid points where the Gaussian density falls below eps_cavity
               IF ((in_memory .OR. dft_control%qs_control%becke_control%save_pot) .AND. &
                   .NOT. dft_control%qs_control%becke_control%dynamic_confine) THEN
                  CALL hfun_zero(dft_control%qs_control%becke_control%cavity%pw%cr3d, &
                                 eps_cavity, just_bounds=.TRUE., bounds=bounds)
                  dft_control%qs_control%becke_control%confine_bounds_int(1) = bounds(1)
                  dft_control%qs_control%becke_control%confine_bounds_int(2) = bounds(2)
               END IF
               ! Optional printing of cavity (meant for testing, so options currently hardcoded...)
               IF (dft_control%qs_control%becke_control%print_cavity) THEN
                  CALL hfun_zero(dft_control%qs_control%becke_control%cavity%pw%cr3d, &
                                 eps_cavity, just_bounds=.FALSE.)
                  NULLIFY (stride, subsys, particles)
                  ALLOCATE (stride(3))
                  stride = (/2, 2, 2/)
                  ! Note PROGRAM_RUN_INFO section neeeds to be active!
                  unit_nr = cp_print_key_unit_nr(logger, becke_restraint_section, "PROGRAM_RUN_INFO", &
                                                 middle_name="BECKE_CAVITY", &
                                                 extension=".cube", file_position="REWIND", &
                                                 log_filename=.FALSE.)
                  IF (para_env%mepos == para_env%source .AND. unit_nr .LT. 1) &
                     CALL cp_abort(__LOCATION__, &
                                   "Please turn on PROGRAM_RUN_INFO to print cavity")
                  CALL get_qs_env(qs_env, subsys=subsys)
                  CALL qs_subsys_get(subsys, particles=particles)
                  CALL cp_pw_to_cube(dft_control%qs_control%becke_control%cavity%pw, &
                                     unit_nr, "CAVITY", particles=particles, &
                                     stride=stride)
                  DEALLOCATE (stride)
               END IF
            END IF
            ! If requested, allocate storage for gradients
            IF (in_memory) THEN
               bo_conf = bo
               ! With confinement active, we dont need to store gradients outside
               ! the confinement bounds since they vanish for all particles
               IF (dft_control%qs_control%becke_control%dynamic_confine .OR. &
                   dft_control%qs_control%becke_control%confine) THEN
                  ! Ensure bounds are within local grid bounds
                  IF (ub_index .GT. bo(2, dir)) THEN
                     tmp_index(2) = bo(2, dir)-1
                  ELSE
                     tmp_index(2) = ub_index-1
                  END IF
                  IF (lb_index .LT. bo(1, dir)) THEN
                     tmp_index(1) = bo(1, dir)+1
                  ELSE
                     tmp_index(1) = lb_index+1
                  END IF
                  bo_conf(2, dir) = tmp_index(2)
                  bo_conf(1, dir) = tmp_index(1)
               ELSE IF (dft_control%qs_control%becke_control%cavity_confine) THEN
                  IF (bounds(2) .LT. bo(2, 3)) bo_conf(2, 3) = bounds(2)-1
                  IF (bounds(1) .GT. bo(1, 3)) bo_conf(1, 3) = bounds(1)+1
               END IF
               ALLOCATE (dft_control%qs_control%becke_control%gradients(3*natom*ngrad, bo_conf(1, 1):bo_conf(2, 1), &
                                                                        bo_conf(1, 2):bo_conf(2, 2), &
                                                                        bo_conf(1, 3):bo_conf(2, 3)))
               gradients => dft_control%qs_control%becke_control%gradients
               gradients = 0.0_dp
               ! If requested, store skipped grid points that were not handled by confinement
               IF (dft_control%qs_control%becke_control%should_skip .AND. .NOT. &
                   dft_control%qs_control%becke_control%cavity_confine) THEN
                  ! With cavity confine, we use the cavity to store these values and dont need this array
                  ALLOCATE (dft_control%qs_control%becke_control%skip_list(bo_conf(1, 1):bo_conf(2, 1), &
                                                                           bo_conf(1, 2):bo_conf(2, 2)+1, &
                                                                           bo_conf(1, 3):bo_conf(2, 3)))
                  dft_control%qs_control%becke_control%skip_list = .FALSE.
               END IF
            END IF
            ! Combined weight: allocate storage for the spin weight function
            IF (dft_control%qs_control%becke_control%constraint_type == cdft_combined_constraint) THEN
               SELECT CASE (dft_control%qs_control%becke_control%combined_type)
               CASE (cdft_combined_all)
                  ! Pass. Charge and spin constraints share same weight
               CASE (cdft_combined_acceptor, cdft_combined_donor)
                  ALLOCATE (is_combined(natom))
                  is_combined = .FALSE.
                  NULLIFY (pw_env, auxbas_pw_pool)
                  CALL get_qs_env(qs_env, pw_env=pw_env)
                  CALL pw_env_get(pw_env, auxbas_rs_grid=rs_cavity, &
                                  auxbas_pw_pool=auxbas_pw_pool)
                  CALL pw_pool_create_pw(auxbas_pw_pool, dft_control%qs_control%becke_control%combined_weight%pw, &
                                         use_data=REALDATA3D, in_space=REALSPACE)
                  dft_control%qs_control%becke_control%combined_weight%pw%cr3d = 0.0_dp
                  DO iatom = 1, dft_control%qs_control%becke_control%natoms
                     IF (dft_control%qs_control%becke_control%combined_type == cdft_combined_acceptor .AND. &
                         dft_control%qs_control%becke_control%coeff(iatom) == -1.0_dp) THEN
                        is_combined(catom(iatom)) = .TRUE.
                     ELSE IF (dft_control%qs_control%becke_control%combined_type == cdft_combined_donor .AND. &
                              dft_control%qs_control%becke_control%coeff(iatom) == +1.0_dp) THEN
                        is_combined(catom(iatom)) = .TRUE.
                     END IF
                  END DO
               END SELECT
            END IF
            ! Calculate weight and possibly the gradient
            DO k = bo(1, 1), bo(2, 1)
               DO j = bo(1, 2), bo(2, 2)
                  DO i = bo(1, 3), bo(2, 3)
                     ! If the grid point is too far from all constraint atoms and any confinement method is active,
                     ! we can skip this grid point as it does not contribute to the weight or gradients
                     IF (dft_control%qs_control%becke_control%cavity_confine) THEN
                        IF (dft_control%qs_control%becke_control%cavity%pw%cr3d(k, j, i) < eps_cavity) CYCLE
                     END IF
                     ind = (/k, j, i/)
                     IF (dft_control%qs_control%becke_control%confine .OR. &
                         dft_control%qs_control%becke_control%dynamic_confine) THEN
                        IF (ind(dir) .LE. lb_index .OR. ind(dir) .GE. ub_index) CYCLE
                     END IF
                     grid_p(1) = k*dr(1)+shift(1)
                     grid_p(2) = j*dr(2)+shift(2)
                     grid_p(3) = i*dr(3)+shift(3)
                     nskipped = 0
                     cell_functions = 1.0_dp
                     skip_me = .FALSE.
                     IF (store_vectors) distances = 0.0_dp
                     IF (in_memory) THEN
                        d_sum_Pm_dR = 0.0_dp
                        d_sum_const_dR = 0.0_dp
                        dP_i_dRi = 0.0_dp
                     END IF
                     ! Iterate over all atoms in the system
                     DO iatom = 1, natom
                        IF (skip_me(iatom)) THEN
                           cell_functions(iatom) = 0.0_dp
                           IF (dft_control%qs_control%becke_control%should_skip) THEN
                              IF (is_constraint(iatom)) nskipped = nskipped+1
                              IF (nskipped == dft_control%qs_control%becke_control%natoms) THEN
                                 IF (in_memory) THEN
                                    IF (dft_control%qs_control%becke_control%cavity_confine) THEN
                                       dft_control%qs_control%becke_control%cavity%pw%cr3d(k, j, i) = 0.0_dp
                                    ELSE
                                       dft_control%qs_control%becke_control%skip_list(k, j, i) = .TRUE.
                                    END IF
                                 END IF
                                 EXIT
                              END IF
                           END IF
                           CYCLE
                        END IF
                        IF (store_vectors) THEN
                           IF (distances(iatom) .EQ. 0.0_dp) THEN
                              r = position_vecs(:, iatom)
                              dist_vec = (r-grid_p)-ANINT((r-grid_p)/cell_v)*cell_v
                              dist1 = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
                              distance_vecs(:, iatom) = dist_vec
                              distances(iatom) = dist1
                           ELSE
                              dist_vec = distance_vecs(:, iatom)
                              dist1 = distances(iatom)
                           END IF
                        ELSE
                           r = particle_set(iatom)%r
                           DO ip = 1, 3
                              r(ip) = MODULO(r(ip), cell%hmat(ip, ip))-cell%hmat(ip, ip)/2._dp
                           END DO
                           dist_vec = (r-grid_p)-ANINT((r-grid_p)/cell_v)*cell_v
                           dist1 = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
                        END IF
                        IF (dist1 .LE. cutoffs(iatom)) THEN
                           IF (in_memory) THEN
                              IF (dist1 .LE. th) dist1 = th
                              dr_i_dR(:) = dist_vec(:)/dist1
                           END IF
                           DO jatom = 1, natom
                              IF (jatom .NE. iatom) THEN
                                 ! Using pairwise symmetry, execute block only for such j<i
                                 ! that have previously not been looped over
                                 ! Note that if skip_me(jatom) = .TRUE., this means that the outer
                                 ! loop over iatom skipped this index when iatom=jatom, but we still
                                 ! need to compute the pair for iatom>jatom
                                 IF (jatom < iatom) THEN
                                    IF (.NOT. skip_me(jatom)) CYCLE
                                 END IF
                                 IF (store_vectors) THEN
                                    IF (distances(jatom) .EQ. 0.0_dp) THEN
                                       r1 = position_vecs(:, jatom)
                                       dist_vec = (r1-grid_p)-ANINT((r1-grid_p)/cell_v)*cell_v
                                       dist2 = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
                                       distance_vecs(:, jatom) = dist_vec
                                       distances(jatom) = dist2
                                    ELSE
                                       dist_vec = distance_vecs(:, jatom)
                                       dist2 = distances(jatom)
                                    END IF
                                 ELSE
                                    r1 = particle_set(jatom)%r
                                    DO ip = 1, 3
                                       r1(ip) = MODULO(r1(ip), cell%hmat(ip, ip))-cell%hmat(ip, ip)/2._dp
                                    END DO
                                    dist_vec = (r1-grid_p)-ANINT((r1-grid_p)/cell_v)*cell_v
                                    dist2 = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
                                 END IF
                                 IF (in_memory) THEN
                                    IF (store_vectors) THEN
                                       dr1_r2 = pair_dist_vecs(:, iatom, jatom)
                                    ELSE
                                       dr1_r2 = (r-r1)-ANINT((r-r1)/cell_v)*cell_v
                                    END IF
                                    IF (dist2 .LE. th) dist2 = th
                                    tmp_const = (R12(iatom, jatom)**3)
                                    dr_ij_dR(:) = dr1_r2(:)/tmp_const

                                    !derivative w.r.t. Rj
                                    dr_j_dR = dist_vec(:)/dist2

                                    dmy_dR_j(:) = -(dr_j_dR(:)/R12(iatom, jatom)-(dist1-dist2)*dr_ij_dR(:))

                                    !derivative w.r.t. Ri
                                    dmy_dR_i(:) = dr_i_dR(:)/R12(iatom, jatom)-(dist1-dist2)*dr_ij_dR(:)
                                 END IF

                                 my1 = (dist1-dist2)/R12(iatom, jatom)

                                 IF (dft_control%qs_control%becke_control%adjust) THEN
                                    my1_homo = my1 ! Homonuclear quantity needed for gradient
                                    my1 = my1+ &
                                          dft_control%qs_control%becke_control%aij(iatom, jatom)*(1.0_dp-my1**2)
                                 END IF
                                 myexp = 1.5_dp*my1-0.5_dp*my1**3

                                 IF (in_memory) THEN
                                    dmyexp = 1.5_dp-1.5_dp*my1**2
                                    tmp_const = (1.5_dp**2)*dmyexp*(1-myexp**2)* &
                                                (1.0_dp-((1.5_dp*myexp-0.5_dp*(myexp**3))**2))

                                    ds_dR_i(:) = -0.5_dp*tmp_const*dmy_dR_i(:) ! d s(myij)/d R_i
                                    ds_dR_j(:) = -0.5_dp*tmp_const*dmy_dR_j(:) ! d s(myij)/d R_j

                                    IF (dft_control%qs_control%becke_control%adjust) THEN
                                       tmp_const = 1.0_dp-2.0_dp*my1_homo* &
                                                   dft_control%qs_control%becke_control%aij(iatom, jatom)
                                       ds_dR_i(:) = ds_dR_i(:)*tmp_const
                                       ! tmp_const is same for both since aij=-aji and myij=-myji
                                       ds_dR_j(:) = ds_dR_j(:)*tmp_const
                                    END IF
                                 END IF

                                 myexp = 1.5_dp*myexp-0.5_dp*myexp**3
                                 myexp = 1.5_dp*myexp-0.5_dp*myexp**3
                                 tmp_const = 0.5_dp*(1.0_dp-myexp) ! s(myij)
                                 cell_functions(iatom) = cell_functions(iatom)*tmp_const
                                 IF (in_memory) THEN
                                    IF (ABS(tmp_const) .LE. th) tmp_const = tmp_const+th
                                    ! P_i independent part of dP_i/dR_i
                                    dP_i_dRi(:, iatom) = dP_i_dRi(:, iatom)+ds_dR_i(:)/tmp_const
                                    ! P_i independent part of dP_i/dR_j
                                    dP_i_dRj(:, iatom, jatom) = ds_dR_j(:)/tmp_const
                                 END IF

                                 IF (dist2 .LE. cutoffs(jatom)) THEN
                                    tmp_const = 0.5_dp*(1.0_dp+myexp) ! s(myji)
                                    cell_functions(jatom) = cell_functions(jatom)*tmp_const
                                    IF (in_memory) THEN
                                       IF (ABS(tmp_const) .LE. th) tmp_const = tmp_const+th
                                       ! P_j independent part of dP_j/dR_i
                                       ! d s(myji)/d R_i = -d s(myij)/d R_i
                                       dP_i_dRj(:, jatom, iatom) = -ds_dR_i(:)/tmp_const
                                       ! P_j independent part of dP_j/dR_j
                                       ! d s(myji)/d R_j = -d s(myij)/d R_j
                                       dP_i_dRi(:, jatom) = dP_i_dRi(:, jatom)-ds_dR_j(:)/tmp_const
                                    END IF
                                 ELSE
                                    skip_me(jatom) = .TRUE.
                                 END IF
                              END IF
                           END DO ! jatom
                           IF (in_memory) THEN
                              dP_i_dRi(:, iatom) = cell_functions(iatom)*dP_i_dRi(:, iatom)
                              d_sum_Pm_dR(:, iatom) = d_sum_Pm_dR(:, iatom)+dP_i_dRi(:, iatom)
                              IF (is_constraint(iatom)) &
                                 d_sum_const_dR(1:3, iatom) = d_sum_const_dR(1:3, iatom)+dP_i_dRi(:, iatom)* &
                                                              coefficients(iatom)
                              IF (ngrad == 2) THEN
                                 IF (is_combined(iatom)) d_sum_const_dR(4:6, iatom) = d_sum_const_dR(4:6, iatom) &
                                                                                      +dP_i_dRi(:, iatom)
                              END IF
                              DO jatom = 1, natom
                                 IF (jatom .NE. iatom) THEN
                                    IF (jatom < iatom) THEN
                                       IF (.NOT. skip_me(jatom)) THEN
                                          dP_i_dRj(:, iatom, jatom) = cell_functions(iatom)*dP_i_dRj(:, iatom, jatom)
                                          d_sum_Pm_dR(:, jatom) = d_sum_Pm_dR(:, jatom)+dP_i_dRj(:, iatom, jatom)
                                          IF (is_constraint(iatom)) &
                                             d_sum_const_dR(1:3, jatom) = d_sum_const_dR(1:3, jatom)+ &
                                                                          dP_i_dRj(:, iatom, jatom)* &
                                                                          coefficients(iatom)
                                          IF (ngrad == 2) THEN
                                             IF (is_combined(iatom)) d_sum_const_dR(4:6, jatom) = &
                                                d_sum_const_dR(4:6, jatom) &
                                                +dP_i_dRj(:, iatom, jatom)
                                          END IF
                                          CYCLE
                                       END IF
                                    END IF
                                    dP_i_dRj(:, iatom, jatom) = cell_functions(iatom)*dP_i_dRj(:, iatom, jatom)
                                    d_sum_Pm_dR(:, jatom) = d_sum_Pm_dR(:, jatom)+dP_i_dRj(:, iatom, jatom)
                                    IF (is_constraint(iatom)) &
                                       d_sum_const_dR(1:3, jatom) = d_sum_const_dR(1:3, jatom)+ &
                                                                    dP_i_dRj(:, iatom, jatom)* &
                                                                    coefficients(iatom)
                                    IF (ngrad == 2) THEN
                                       IF (is_combined(iatom)) d_sum_const_dR(4:6, jatom) = d_sum_const_dR(4:6, jatom) &
                                                                                            +dP_i_dRj(:, iatom, jatom)
                                    END IF
                                 END IF
                              END DO
                           END IF
                        ELSE
                           cell_functions(iatom) = 0.0_dp
                           skip_me(iatom) = .TRUE.
                           IF (dft_control%qs_control%becke_control%should_skip) THEN
                              IF (is_constraint(iatom)) nskipped = nskipped+1
                              IF (nskipped == dft_control%qs_control%becke_control%natoms) THEN
                                 IF (in_memory) THEN
                                    IF (dft_control%qs_control%becke_control%cavity_confine) THEN
                                       dft_control%qs_control%becke_control%cavity%pw%cr3d(k, j, i) = 0.0_dp
                                    ELSE
                                       dft_control%qs_control%becke_control%skip_list(k, j, i) = .TRUE.
                                    END IF
                                 END IF
                                 EXIT
                              END IF
                           END IF
                        END IF
                     END DO !iatom
                     IF (nskipped == dft_control%qs_control%becke_control%natoms) CYCLE
                     sum_cell_f_constr = 0.0_dp
                     sum_cell_f_combined = 0.0_dp
                     DO ip = 1, dft_control%qs_control%becke_control%natoms
                        sum_cell_f_constr = sum_cell_f_constr+cell_functions(catom(ip))* &
                                            dft_control%qs_control%becke_control%coeff(ip)
                        IF (ngrad == 2) THEN
                           ! The coeff is implicitly 1.0
                           IF (is_combined(catom(ip))) &
                              sum_cell_f_combined = sum_cell_f_combined+cell_functions(catom(ip))
                        END IF
                     END DO

                     sum_cell_f_all = 0.0_dp
                     DO ip = 1, natom
                        sum_cell_f_all = sum_cell_f_all+cell_functions(ip)
                     END DO
                     ! Gradients at (k,j,i)
                     IF (in_memory .AND. ABS(sum_cell_f_all) .GT. 0.0_dp) THEN
                        DO iatom = 1, natom
                           gradients(3*(iatom-1)+1:3*(iatom-1)+3, k, j, i) = &
                              d_sum_const_dR(1:3, iatom)/sum_cell_f_all-sum_cell_f_constr* &
                              d_sum_Pm_dR(1:3, iatom)/(sum_cell_f_all**2)
                           IF (ngrad == 2) THEN
                              gradients(3*(iatom-1)+1+natom:3*(iatom-1)+3+natom, k, j, i) = &
                                 d_sum_const_dR(4:6, iatom)/sum_cell_f_all-sum_cell_f_constr* &
                                 d_sum_Pm_dR(1:3, iatom)/(sum_cell_f_all**2)
                           END IF
                        END DO
                     END IF
                     ! Weight function at (k,j,i)
                     IF (ABS(sum_cell_f_all) .GT. 0.000001) THEN
                        becke_const%pw%cr3d(k, j, i) = sum_cell_f_constr/sum_cell_f_all
                        IF (dft_control%qs_control%becke_control%atomic_charges) THEN
                           DO iatom = 1, dft_control%qs_control%becke_control%natoms
                              charge(iatom)%pw%cr3d(k, j, i) = cell_functions(catom(iatom))/sum_cell_f_all
                           END DO
                        END IF
                        IF (ngrad == 2) THEN
                           dft_control%qs_control%becke_control%combined_weight%pw%cr3d(k, j, i) = sum_cell_f_combined/ &
                                                                                                   sum_cell_f_all
                        END IF
                     END IF
                  END DO
               END DO
            END DO
            ! Release storage
            IF (ALLOCATED(coefficients)) &
               DEALLOCATE (coefficients)
            IF (in_memory) THEN
               DEALLOCATE (ds_dR_j)
               DEALLOCATE (ds_dR_i)
               DEALLOCATE (d_sum_Pm_dR)
               DEALLOCATE (d_sum_const_dR)
               DEALLOCATE (dP_i_dRj)
               DEALLOCATE (dP_i_dRi)
               NULLIFY (gradients)
               IF (store_vectors) THEN
                  DEALLOCATE (pair_dist_vecs)
               END IF
            END IF
            NULLIFY (cutoffs)
            IF (ALLOCATED(is_constraint)) &
               DEALLOCATE (is_constraint)
            DEALLOCATE (catom)
            DEALLOCATE (R12)
            DEALLOCATE (cell_functions)
            DEALLOCATE (skip_me)
            IF (store_vectors) THEN
               DEALLOCATE (distances)
               DEALLOCATE (distance_vecs)
               DEALLOCATE (position_vecs)
            END IF
            IF (ALLOCATED(is_combined)) DEALLOCATE (is_combined)
         END IF
         bo = becke_const%pw%pw_grid%bounds_local
         IF (dft_control%qs_control%becke_control%atomic_charges) THEN
            ALLOCATE (electronic_charge(dft_control%qs_control%becke_control%natoms, dft_control%nspins))
            electronic_charge = 0.0_dp
         END IF
         ! Calculate value of constraint i.e. int ( rho(r) w(r) dr)
         DO i = 1, dft_control%nspins
            DO ivar = 1, nvar
               ! With external control, we can use cavity_mat as a mask to kahan sum
               IF (dft_control%qs_control%becke_control%external_control .AND. &
                   dft_control%qs_control%becke_control%cavity_confine) THEN
                  IF (ivar == 2) THEN
                     sign = 1.0_dp
                     IF (i == 2) sign = -1.0_dp
                     SELECT CASE (dft_control%qs_control%becke_control%combined_type)
                     CASE (cdft_combined_all)
                        dE(ivar) = dE(ivar)+ &
                                   sign*accurate_sum_arrays_product(becke_const%pw%cr3d, rho_r(i)%pw%cr3d, &
                                                                    dft_control%qs_control%becke_control%cavity_mat, &
                                                                    eps_cavity)*dvol
                     CASE (cdft_combined_donor, cdft_combined_acceptor)
                        ! Transfer spin constraint weight truncated real-3D matrix -> pw%cr3d
                        IF (ASSOCIATED(dft_control%qs_control%becke_control%combined_mat)) THEN
                           NULLIFY (pw_env, auxbas_pw_pool)
                           CALL get_qs_env(qs_env, pw_env=pw_env)
                           CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
                           CALL pw_pool_create_pw(auxbas_pw_pool, &
                                                  dft_control%qs_control%becke_control%combined_weight%pw, &
                                                  use_data=REALDATA3D, in_space=REALSPACE)
                           dft_control%qs_control%becke_control%combined_weight%pw%cr3d = 0.0_dp
                           CALL pw_set(dft_control%qs_control%becke_control%combined_weight%pw, &
                                       dft_control%qs_control%becke_control%combined_mat)
                           DEALLOCATE (dft_control%qs_control%becke_control%combined_mat)
                        END IF
                        dE(ivar) = dE(ivar)+sign* &
                                   accurate_sum_arrays_product(dft_control%qs_control%becke_control% &
                                                               combined_weight%pw%cr3d, &
                                                               rho_r(i)%pw%cr3d, &
                                                               dft_control%qs_control%becke_control%cavity_mat, &
                                                               eps_cavity)*dvol
                     END SELECT
                  ELSE
                     sign = 1.0_dp
                     ! Note the negative sign for beta spin if spin constraint is active
                     IF (dft_control%qs_control%becke_control%constraint_type == &
                         cdft_magnetization_constraint .AND. i == 2) sign = -1.0_dp
                     dE(ivar) = dE(ivar)+sign* &
                                accurate_sum_arrays_product(becke_const%pw%cr3d, rho_r(i)%pw%cr3d, &
                                                            dft_control%qs_control%becke_control%cavity_mat, &
                                                            eps_cavity)*dvol
                  END IF
               ELSE
                  IF (ivar == 2) THEN
                     sign = 1.0_dp
                     IF (i == 2) sign = -1.0_dp
                     SELECT CASE (dft_control%qs_control%becke_control%combined_type)
                     CASE (cdft_combined_all)
                        dE(ivar) = dE(ivar)+sign*accurate_sum(becke_const%pw%cr3d*rho_r(i)%pw%cr3d)*dvol
                     CASE (cdft_combined_donor, cdft_combined_acceptor)
                        dE(ivar) = dE(ivar)+sign* &
                                   accurate_sum(dft_control%qs_control%becke_control%combined_weight%pw%cr3d* &
                                                rho_r(i)%pw%cr3d)*dvol
                     END SELECT
                  ELSE
                     sign = 1.0_dp
                     ! Note the negative sign for beta spin if spin constraint is active
                     IF (dft_control%qs_control%becke_control%constraint_type == &
                         cdft_magnetization_constraint .AND. i == 2) sign = -1.0_dp
                     dE(ivar) = dE(ivar)+sign*accurate_sum(becke_const%pw%cr3d*rho_r(i)%pw%cr3d)*dvol
                  END IF
               END IF
            END DO
            IF (dft_control%qs_control%becke_control%atomic_charges) THEN
               DO iatom = 1, dft_control%qs_control%becke_control%natoms
                  electronic_charge(iatom, i) = accurate_sum(charge(iatom)%pw%cr3d*rho_r(i)%pw%cr3d)*dvol
               END DO
            END IF
         END DO
         CALL get_qs_env(qs_env, energy=energy)
         CALL mp_sum(dE, para_env%group)
         IF (dft_control%qs_control%becke_control%atomic_charges) THEN
            CALL mp_sum(electronic_charge, para_env%group)
         END IF
         ! Use fragment densities as reference value (= Becke deformation density)
         IF (dft_control%qs_control%becke_control%fragment_density) THEN
            ! Fragment densities are meaningful only for some calculation types
            IF (.NOT. qs_env%single_point_run) &
               CALL cp_abort(__LOCATION__, &
                             "Becke deformation density only compatible with single "// &
                             "point calculations (run_type ENERGY or ENERGY_FORCE).")
            IF (dft_control%qs_control%gapw) &
               CALL cp_abort(__LOCATION__, &
                             "Becke deformation density not compatible with GAPW.")
            SELECT CASE (dft_control%qs_control%becke_control%constraint_type)
            CASE (cdft_density_constraint)
               ! Pass
            CASE (cdft_magnetization_constraint)
               CALL cp_abort(__LOCATION__, &
                             "Becke deformation density not yet compatible with magnetization density constraint.")
            CASE (cdft_combined_constraint)
               CALL cp_abort(__LOCATION__, &
                             "Becke deformation density not yet compatible with combined charge+spin constraint.")
            END SELECT
            IF (.NOT. dft_control%qs_control%becke_control%fragments_integrated) THEN
               ALLOCATE (dft_control%qs_control%becke_control%fragments(2))
               CALL get_qs_env(qs_env, pw_env=pw_env)
               CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
               ! Read fragment reference densities (total density => no spin density data)
               ! This step can be very slow if the rs_grid is dense
               CALL pw_pool_create_pw(auxbas_pw_pool, dft_control%qs_control%becke_control%fragments(1)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL cube_to_pw(dft_control%qs_control%becke_control%fragments(1)%pw, &
                               dft_control%qs_control%becke_control%fragment_a_fname, 1.0_dp)
               CALL pw_pool_create_pw(auxbas_pw_pool, dft_control%qs_control%becke_control%fragments(2)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL cube_to_pw(dft_control%qs_control%becke_control%fragments(2)%pw, &
                               dft_control%qs_control%becke_control%fragment_b_fname, 1.0_dp)
               ! Sum up fragments
               CALL pw_pool_create_pw(auxbas_pw_pool, rho_frag%pw, use_data=REALDATA3D, &
                                      in_space=REALSPACE)
               CALL pw_copy(dft_control%qs_control%becke_control%fragments(1)%pw, rho_frag%pw)
               CALL pw_axpy(dft_control%qs_control%becke_control%fragments(2)%pw, rho_frag%pw, 1.0_dp)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, dft_control%qs_control%becke_control%fragments(1)%pw)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, dft_control%qs_control%becke_control%fragments(2)%pw)
               DEALLOCATE (dft_control%qs_control%becke_control%fragments)
               ! Check that the number of electrons is consistent
               CALL get_qs_env(qs_env, subsys=subsys)
               CALL qs_subsys_get(subsys, nelectron_total=nelectron_total)
               nelectron_frag = pw_integrate_function(rho_frag%pw)
               IF (NINT(nelectron_frag) /= nelectron_total) &
                  CALL cp_abort(__LOCATION__, &
                                "The number of electrons in the reference and interacting "// &
                                " configurations does not match. Check your fragment cube files.")
               ! Update constraint target value i.e. perform integration w*rho_frag*dr
               ! Note, we assume the weight w is defined only on one of the fragments
               dft_control%qs_control%becke_control%target = 0.0_dp
               IF (ANY(dft_control%qs_control%becke_control%coeff .NE. 1.0_dp)) &
                  CALL cp_abort(__LOCATION__, &
                                "Absolute weight required for Becke deformation density.")
               IF (dft_control%qs_control%becke_control%external_control .AND. &
                   dft_control%qs_control%becke_control%cavity_confine) THEN
                  dft_control%qs_control%becke_control%target(1) = &
                     dft_control%qs_control%becke_control%target(1)+ &
                     accurate_sum_arrays_product(becke_const%pw%cr3d, rho_frag%pw%cr3d, &
                                                 dft_control%qs_control%becke_control%cavity_mat, eps_cavity)*dvol
               ELSE
                  dft_control%qs_control%becke_control%target(1) = &
                     dft_control%qs_control%becke_control%target(1)+ &
                     accurate_sum(becke_const%pw%cr3d*rho_frag%pw%cr3d)*dvol
               END IF
               CALL mp_sum(dft_control%qs_control%becke_control%target, para_env%group)
               IF (debug_this_module) THEN
                  ALLOCATE (stride(3))
                  stride = (/1, 1, 1/)
                  unit_nr = cp_print_key_unit_nr(logger, becke_restraint_section, "PROGRAM_RUN_INFO", &
                                                 middle_name="BECKE_WEIGHT", &
                                                 extension=".cube", file_position="REWIND", &
                                                 log_filename=.FALSE.)
                  CALL get_qs_env(qs_env, subsys=subsys)
                  CALL qs_subsys_get(subsys, particles=particles)
                  CALL cp_pw_to_cube(becke_const%pw, &
                                     unit_nr, "BECKE_WEIGHT", particles=particles, &
                                     stride=stride)
                  unit_nr = cp_print_key_unit_nr(logger, becke_restraint_section, "PROGRAM_RUN_INFO", &
                                                 middle_name="FRAG_DENS", &
                                                 extension=".cube", file_position="REWIND", &
                                                 log_filename=.FALSE.)
                  CALL cp_pw_to_cube(rho_frag%pw, &
                                     unit_nr, "FRAG_DENS", particles=particles, &
                                     stride=stride)
                  DEALLOCATE (stride)
               END IF
               ! Calculate reference atomic charges int( w_i * rho_frag * dr )
               IF (dft_control%qs_control%becke_control%atomic_charges) THEN
                  ALLOCATE (dft_control%qs_control%becke_control% &
                            charges_fragment(dft_control%qs_control%becke_control%natoms))
                  DO iatom = 1, dft_control%qs_control%becke_control%natoms
                     dft_control%qs_control%becke_control%charges_fragment(iatom) = &
                        accurate_sum(charge(iatom)%pw%cr3d*rho_frag%pw%cr3d)*dvol
                  END DO
                  CALL mp_sum(dft_control%qs_control%becke_control%charges_fragment, para_env%group)
               END IF
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_frag%pw)
               target_val(:) = dft_control%qs_control%becke_control%target(:)
               IF (dft_control%qs_control%cdft) dft_control%qs_control%cdft_control%target = target_val
               dft_control%qs_control%becke_control%fragments_integrated = .TRUE.
            END IF
         END IF
         IF (dft_control%qs_control%gapw) THEN
            ! GAPW: add core charges (rho_hard - rho_soft)
            IF (dft_control%nspins == 2 .AND. nvar == 2) THEN
               ALLOCATE (gapw_offset(dft_control%nspins+nvar))
            ELSE
               ALLOCATE (gapw_offset(dft_control%nspins))
            END IF
            gapw_offset = 0.0_dp
            CALL get_qs_env(qs_env, rho0_mpole=rho0_mpole)
            CALL get_rho0_mpole(rho0_mpole, mp_rho=mp_rho)
            DO iatom = 1, dft_control%qs_control%becke_control%natoms
               jatom = dft_control%qs_control%becke_control%atoms(iatom)
               CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), paw_atom=paw_atom)
               IF (paw_atom) THEN
                  DO i = 1, dft_control%nspins
                     gapw_offset(i) = gapw_offset(i)+ &
                                      dft_control%qs_control%becke_control%coeff(iatom)*mp_rho(jatom)%q0(i)
                     IF (nvar == 2 .AND. is_combined(jatom)) THEN
                        sign = 1.0_dp
                        IF (i == 2) sign = -1.0_dp
                        gapw_offset(i+nvar) = gapw_offset(i+nvar)+sign*mp_rho(jatom)%q0(i)
                     END IF
                     IF (dft_control%qs_control%becke_control%atomic_charges) &
                        electronic_charge(iatom, i) = electronic_charge(iatom, i)+mp_rho(jatom)%q0(i)
                  END DO
               END IF
            END DO
            DO i = 1, dft_control%nspins
               DO ivar = 1, nvar
                  IF (dft_control%qs_control%becke_control%constraint_type == &
                      cdft_magnetization_constraint .AND. i == 2) THEN
                     dE(ivar) = dE(ivar)-gapw_offset(i)
                  ELSE IF (ivar == 2) THEN
                     dE(ivar) = dE(ivar)+gapw_offset(i+nvar)
                  ELSE
                     dE(ivar) = dE(ivar)+gapw_offset(i)
                  END IF
               END DO
            END DO
            DEALLOCATE (gapw_offset)
         END IF

         dft_control%qs_control%becke_control%becke_order_p(:) = dE(:)
         energy%becke = 0.0_dp
         DO ivar = 1, nvar
            energy%becke = energy%becke+(dE(ivar)-target_val(ivar))*strength(ivar)
         END DO

         IF (dft_control%qs_control%cdft) &
            dft_control%qs_control%cdft_control%value(:) = dE(:)

         IF (iw > 0) THEN
            WRITE (iw, '(/,T3,A,T60)') &
               '--------------------- Becke constraint information ---------------------'
            SELECT CASE (dft_control%qs_control%becke_control%constraint_type)
            CASE (cdft_density_constraint, cdft_magnetization_constraint)
               WRITE (iw, '(T3,A,T54,(3X,F18.12))') &
                  'Target value of constraint         :', target_val(1)
               WRITE (iw, '(T3,A,T54,(3X,F18.12))') &
                  'Current value of constraint        :', dE(1)
               WRITE (iw, '(T3,A,T59,(3X,ES13.3))') &
                  'Deviation from target              :', ABS(dE(1)-target_val(1))
               WRITE (iw, '(T3,A,T54,(3X,F18.12))') &
                  'Strength of constraint             :', strength(1)
            CASE (cdft_combined_constraint)
               WRITE (iw, '(T3,A,T54,(3X,F18.12))') &
                  'Target value of charge constraint :', target_val(1)
               WRITE (iw, '(T3,A,T54,(3X,F18.12))') &
                  'Target value of spin constraint   :', target_val(2)
               WRITE (iw, '(T3,A,T54,(3X,F18.12))') &
                  'Current value of charge constraint:', dE(1)
               WRITE (iw, '(T3,A,T54,(3X,F18.12))') &
                  'Current value of spin constraint  :', dE(2)
               WRITE (iw, '(T3,A,T59,(3X,ES13.3))') &
                  'Deviation from target (charge)    :', ABS(dE(1)-target_val(1))
               WRITE (iw, '(T3,A,T59,(3X,ES13.3))') &
                  'Deviation from target (spin)      :', ABS(dE(2)-target_val(2))
               WRITE (iw, '(T3,A,T54,(3X,F18.12))') &
                  'Strength of constraint (charge)   :', strength(1)
               WRITE (iw, '(T3,A,T54,(3X,F18.12))') &
                  'Strength of constraint (spin)     :', strength(2)
            END SELECT
            WRITE (iw, '(T3,A)') &
               '------------------------------------------------------------------------'
            IF (dft_control%qs_control%becke_control%atomic_charges) THEN
               WRITE (iw, '(/,T3,A,T60)') &
                  '-------------------------------- Becke atomic charges -------------------------------'
               IF (.NOT. dft_control%qs_control%becke_control%fragment_density) THEN
                  IF (dft_control%nspins == 1) THEN
                     WRITE (iw, '(/,T3,A,A)') &
                        '#Atom  Element   Coefficient', '     Core charge    Population (total)'// &
                        '    Net charge'
                     tc = 0.0_dp
                     DO iatom = 1, dft_control%qs_control%becke_control%natoms
                        jatom = dft_control%qs_control%becke_control%atoms(iatom)
                        CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, &
                                             element_symbol=element_symbol, &
                                             kind_number=ikind)
                        CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
                        WRITE (iw, "(i7,T15,A2,T23,F8.3,T39,F8.3,T61,F8.3,T75,F8.3)") &
                           jatom, element_symbol, &
                           dft_control%qs_control%becke_control%coeff(iatom), &
                           zeff, electronic_charge(iatom, 1), &
                           (zeff-electronic_charge(iatom, 1))
                        tc = tc+(zeff-electronic_charge(iatom, 1))
                     END DO
                     WRITE (iw, '(/,T3,A,T75,F8.3,/)') "Total Charge: ", tc
                  ELSE
                     WRITE (iw, '(/,T3,A,A)') &
                        '#Atom  Element   Coefficient', '     Core charge    Population (alpha, beta)'// &
                        '    Net charge'
                     tc = 0.0_dp
                     DO iatom = 1, dft_control%qs_control%becke_control%natoms
                        jatom = dft_control%qs_control%becke_control%atoms(iatom)
                        CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, &
                                             element_symbol=element_symbol, &
                                             kind_number=ikind)
                        CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
                        WRITE (iw, "(i7,T15,A2,T23,F8.3,T39,F8.3,T53,F8.3,T63,F8.3,T81,F8.3)") &
                           jatom, element_symbol, &
                           dft_control%qs_control%becke_control%coeff(iatom), &
                           zeff, electronic_charge(iatom, 1), electronic_charge(iatom, 2), &
                           (zeff-electronic_charge(iatom, 1)-electronic_charge(iatom, 2))
                        tc = tc+(zeff-electronic_charge(iatom, 1)-electronic_charge(iatom, 2))
                     END DO
                     WRITE (iw, '(/,T3,A,T81,F8.3,/)') "Total Charge: ", tc
                  END IF
               ELSE
                  WRITE (iw, '(/,T3,A,A)') &
                     '#Atom  Element   Coefficient', '     Fragment charge    Population (total)'// &
                     '    Net charge'
                  tc = 0.0_dp
                  DO iatom = 1, dft_control%qs_control%becke_control%natoms
                     jatom = dft_control%qs_control%becke_control%atoms(iatom)
                     CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, &
                                          element_symbol=element_symbol, &
                                          kind_number=ikind)
                     IF (dft_control%nspins == 1) THEN
                        WRITE (iw, "(i7,T15,A2,T23,F8.3,T39,F8.3,T61,F8.3,T75,F8.3)") &
                           jatom, element_symbol, &
                           dft_control%qs_control%becke_control%coeff(iatom), &
                           dft_control%qs_control%becke_control%charges_fragment(iatom), &
                           electronic_charge(iatom, 1), &
                           (electronic_charge(iatom, 1)- &
                            dft_control%qs_control%becke_control%charges_fragment(iatom))
                        tc = tc+(electronic_charge(iatom, 1)- &
                                 dft_control%qs_control%becke_control%charges_fragment(iatom))
                     ELSE
                        WRITE (iw, "(i7,T15,A2,T23,F8.3,T39,F8.3,T61,F8.3,T75,F8.3)") &
                           jatom, element_symbol, &
                           dft_control%qs_control%becke_control%coeff(iatom), &
                           dft_control%qs_control%becke_control%charges_fragment(iatom), &
                           electronic_charge(iatom, 1)+electronic_charge(iatom, 2), &
                           (electronic_charge(iatom, 1)+electronic_charge(iatom, 2)- &
                            dft_control%qs_control%becke_control%charges_fragment(iatom))
                        tc = tc+(electronic_charge(iatom, 1)+electronic_charge(iatom, 2)- &
                                 dft_control%qs_control%becke_control%charges_fragment(iatom))
                     END IF
                  END DO
                  WRITE (iw, '(/,T3,A,T75,F8.3,/)') "Total Charge: ", tc
               END IF
            END IF
         END IF
         DEALLOCATE (dE, strength, target_val)

         IF (dft_control%qs_control%becke_control%atomic_charges) &
            DEALLOCATE (electronic_charge)

         CALL cp_print_key_finished_output(iw, logger, becke_restraint_section, &
                                           "PROGRAM_RUN_INFO")

         IF (calculate_forces) CALL becke_force(qs_env, becke_const)

      END IF
      CALL timestop(handle)

   END SUBROUTINE becke_restraint

! **************************************************************************************************
!> \brief calculates Becke constraint gradients
!> \param qs_env the qs_env where to build the gradients
!> \param becke_const the Becke real space potential
!> \par   History
!>        Created 01.2007 [fschiff]
!>        Extended functionality 12/15-12/16 [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE becke_force(qs_env, becke_const)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type)                                    :: becke_const

      CHARACTER(len=*), PARAMETER :: routineN = 'becke_force', routineP = moduleN//':'//routineN

      INTEGER :: dir, handle, i, iatom, ikind, ind(3), ip, ispin, ivar, j, jatom, k, lb_index, &
         natom, ngrad, np(3), nskipped, nvar, tmp_index(2), ub_index
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, catom, kind_of
      INTEGER, DIMENSION(2, 3)                           :: bo, bo_conf
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: is_combined, is_constraint
      REAL(kind=dp) :: dist1, dist2, dpolyn, dvol, eps_cavity, f_3, my1, my1_homo, polyn, sign, &
         sum_cell_f_combined, sum_cell_f_constr, sum_Pm, th, tmp_const
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: coefficients, P_i, strength
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: d_sum_const_dR, d_sum_Pm_dR, dP_i_dRi, &
                                                            dw_dR, integrated, R12, s_my_ij
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dP_i_dRj, ds_dR
      REAL(kind=dp), DIMENSION(3)                        :: cell_v, dist_vec_i, dist_vec_j, &
                                                            dmy_dR_i, dmy_dR_j, dr, dr1_r2, &
                                                            dr_i_dR, dr_ij_dR, dr_j_dR, grid_p, r, &
                                                            r1, shift
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cutoffs
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: gradients
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)
      NULLIFY (atomic_kind_set, cell, para_env, dft_control, particle_set, &
               rho, rho_r, force, cutoffs)

      CALL get_qs_env(qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      natom=natom, &
                      particle_set=particle_set, &
                      cell=cell, &
                      rho=rho, &
                      force=force, &
                      dft_control=dft_control, &
                      para_env=para_env)
      CALL qs_rho_get(rho, rho_r=rho_r)

      th = 1.0e-8_dp
      IF (dft_control%qs_control%becke_restraint) THEN
         SELECT CASE (dft_control%qs_control%becke_control%constraint_type)
         CASE (cdft_density_constraint, cdft_magnetization_constraint)
            nvar = 1
            ngrad = 1
         CASE (cdft_combined_constraint)
            nvar = 2
            SELECT CASE (dft_control%qs_control%becke_control%combined_type)
            CASE (cdft_combined_all)
               ngrad = 1
            CASE (cdft_combined_acceptor, cdft_combined_donor)
               ngrad = 2
            END SELECT
         END SELECT
         ALLOCATE (strength(nvar))
         strength(:) = dft_control%qs_control%becke_control%strength(:)
         cutoffs => dft_control%qs_control%becke_control%cutoffs
         eps_cavity = dft_control%qs_control%becke_control%eps_cavity
         ALLOCATE (atom_of_kind(natom))
         ALLOCATE (kind_of(natom))

         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  atom_of_kind=atom_of_kind, &
                                  kind_of=kind_of)

         bo = becke_const%pw%pw_grid%bounds_local
         np = becke_const%pw%pw_grid%npts
         dvol = becke_const%pw%pw_grid%dvol
         dr = becke_const%pw%pw_grid%dr
         np = becke_const%pw%pw_grid%npts
         shift = -REAL(MODULO(np, 2), dp)*dr/2.0_dp
         ALLOCATE (integrated(3*ngrad, natom))

         IF (dft_control%qs_control%becke_control%confine .OR. &
             dft_control%qs_control%becke_control%dynamic_confine) THEN
            dir = dft_control%qs_control%becke_control%confine_dir
            lb_index = dft_control%qs_control%becke_control%confine_bounds_int(1)
            ub_index = dft_control%qs_control%becke_control%confine_bounds_int(2)
         END IF

         IF (dft_control%qs_control%becke_control%in_memory) THEN
            gradients => dft_control%qs_control%becke_control%gradients
            bo_conf = bo
            IF (dft_control%qs_control%becke_control%dynamic_confine .OR. &
                dft_control%qs_control%becke_control%confine) THEN
               IF (ub_index .GT. bo(2, dir)) THEN
                  tmp_index(2) = bo(2, dir)-1
               ELSE
                  tmp_index(2) = ub_index-1
               END IF
               IF (lb_index .LT. bo(1, dir)) THEN
                  tmp_index(1) = bo(1, dir)+1
               ELSE
                  tmp_index(1) = lb_index+1
               END IF
               bo_conf(2, dir) = tmp_index(2)
               bo_conf(1, dir) = tmp_index(1)
            ELSE IF (dft_control%qs_control%becke_control%cavity_confine) THEN
               IF (dft_control%qs_control%becke_control%confine_bounds_int(2) .LT. bo(2, 3)) &
                  bo_conf(2, 3) = dft_control%qs_control%becke_control%confine_bounds_int(2)-1
               IF (dft_control%qs_control%becke_control%confine_bounds_int(1) .GT. bo(1, 3)) &
                  bo_conf(1, 3) = dft_control%qs_control%becke_control%confine_bounds_int(1)+1
            END IF
         END IF

         integrated = 0.0_dp

         IF (dft_control%qs_control%becke_control%in_memory) THEN
            ! Gradients have been precomputed, so we can just integrate
            IF (.NOT. ASSOCIATED(dft_control%qs_control%becke_control%cavity_mat)) THEN
               ! No external control
               DO k = bo(1, 1), bo(2, 1)
                  DO j = bo(1, 2), bo(2, 2)
                     DO i = bo(1, 3), bo(2, 3)
                        ! First check if this grid point should be skipped
                        IF (dft_control%qs_control%becke_control%cavity_confine) THEN
                           IF (dft_control%qs_control%becke_control%cavity%pw%cr3d(k, j, i) < eps_cavity) CYCLE
                        END IF
                        ind = (/k, j, i/)
                        IF (dft_control%qs_control%becke_control%confine .OR. &
                            dft_control%qs_control%becke_control%dynamic_confine) THEN
                           IF (ind(dir) .LE. lb_index .OR. ind(dir) .GE. ub_index) CYCLE
                        END IF
                        IF (dft_control%qs_control%becke_control%should_skip .AND. .NOT. &
                            dft_control%qs_control%becke_control%cavity_confine) THEN
                           IF (dft_control%qs_control%becke_control%skip_list(k, j, i)) CYCLE
                        END IF
                        DO iatom = 1, natom
                           DO ispin = 1, dft_control%nspins
                              DO ivar = 1, nvar
                                 sign = 1.0_dp
                                 IF (ivar == 2) THEN
                                    IF (ispin == 2) sign = -1.0_dp
                                    integrated(1+3*(ngrad-1):3+3*(ngrad-1), iatom) = &
                                       integrated(1+3*(ngrad-1):3+3*(ngrad-1), iatom)+sign* &
                                       gradients(3*(iatom-1)+1+natom*(ngrad-1):3*(iatom-1)+3+natom*(ngrad-1), k, j, i)* &
                                       rho_r(ispin)%pw%cr3d(k, j, i)*dvol
                                 ELSE
                                    IF (dft_control%qs_control%becke_control%constraint_type == &
                                        cdft_magnetization_constraint .AND. ispin == 2) sign = -1.0_dp
                                    integrated(1:3, iatom) = integrated(1:3, iatom)+sign* &
                                                             gradients(3*(iatom-1)+1:3*(iatom-1)+3, k, j, i)* &
                                                             rho_r(ispin)%pw%cr3d(k, j, i)*dvol
                                 END IF
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            ELSE
               DO k = LBOUND(dft_control%qs_control%becke_control%cavity_mat, 1), &
                      UBOUND(dft_control%qs_control%becke_control%cavity_mat, 1)
                  DO j = LBOUND(dft_control%qs_control%becke_control%cavity_mat, 2), &
                         UBOUND(dft_control%qs_control%becke_control%cavity_mat, 2)
                     DO i = LBOUND(dft_control%qs_control%becke_control%cavity_mat, 3), &
                            UBOUND(dft_control%qs_control%becke_control%cavity_mat, 3)
                        ! First check if this grid point should be skipped
                        IF (dft_control%qs_control%becke_control%cavity_mat(k, j, i) < eps_cavity) CYCLE
                        DO iatom = 1, natom
                           DO ispin = 1, dft_control%nspins
                              DO ivar = 1, nvar
                                 sign = 1.0_dp
                                 IF (ivar == 2) THEN
                                    IF (ispin == 2) sign = -1.0_dp
                                    integrated(1+3*(ngrad-1):3+3*(ngrad-1), iatom) = &
                                       integrated(1+3*(ngrad-1):3+3*(ngrad-1), iatom)+sign* &
                                       gradients(3*(iatom-1)+1+natom*(ngrad-1):3*(iatom-1)+3+natom*(ngrad-1), k, j, i)* &
                                       rho_r(ispin)%pw%cr3d(k, j, i)*dvol
                                 ELSE
                                    IF (dft_control%qs_control%becke_control%constraint_type == &
                                        cdft_magnetization_constraint .AND. ispin == 2) sign = -1.0_dp
                                    integrated(1:3, iatom) = integrated(1:3, iatom)+sign* &
                                                             gradients(3*(iatom-1)+1:3*(iatom-1)+3, k, j, i)* &
                                                             rho_r(ispin)%pw%cr3d(k, j, i)*dvol
                                 END IF
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END IF
            NULLIFY (gradients)
            IF (.NOT. dft_control%qs_control%becke_control%transfer_pot) &
               DEALLOCATE (dft_control%qs_control%becke_control%gradients)
            IF (dft_control%qs_control%becke_control%should_skip .AND. .NOT. &
                dft_control%qs_control%becke_control%cavity_confine) &
               DEALLOCATE (dft_control%qs_control%becke_control%skip_list)
         ELSE
            ! Construct gradients and then integrate
            ALLOCATE (catom(dft_control%qs_control%becke_control%natoms))
            ALLOCATE (is_constraint(natom))
            ALLOCATE (coefficients(natom))

            is_constraint = .FALSE.
            DO i = 1, dft_control%qs_control%becke_control%natoms
               catom(i) = dft_control%qs_control%becke_control%atoms(i)
               is_constraint(dft_control%qs_control%becke_control%atoms(i)) = .TRUE.
               coefficients(catom(i)) = dft_control%qs_control%becke_control%coeff(i)
            ENDDO
            IF (dft_control%qs_control%becke_control%constraint_type == cdft_combined_constraint) THEN
               SELECT CASE (dft_control%qs_control%becke_control%combined_type)
               CASE (cdft_combined_all)
                  ! Do nothing
               CASE (cdft_combined_acceptor, cdft_combined_donor)
                  ALLOCATE (is_combined(natom))
                  is_combined = .FALSE.
                  DO i = 1, dft_control%qs_control%becke_control%natoms
                     IF (dft_control%qs_control%becke_control%combined_type == cdft_combined_acceptor .AND. &
                         dft_control%qs_control%becke_control%coeff(i) == -1.0_dp) THEN
                        is_combined(catom(i)) = .TRUE.
                     ELSE IF (dft_control%qs_control%becke_control%combined_type == cdft_combined_donor .AND. &
                              dft_control%qs_control%becke_control%coeff(i) == +1.0_dp) THEN
                        is_combined(catom(i)) = .TRUE.
                     END IF
                  END DO
               END SELECT
            END IF
            !   calculate distances from target atom (only implemented for a diatomic system)
            DO i = 1, 3
               cell_v(i) = cell%hmat(i, i)
            END DO

            ALLOCATE (R12(natom, natom))
            ALLOCATE (s_my_ij(natom, natom))
            ALLOCATE (ds_dR(3, natom, natom))
            ALLOCATE (P_i(natom))
            ALLOCATE (dw_dR(3*ngrad, natom))
            ALLOCATE (d_sum_Pm_dR(3, natom))
            ALLOCATE (d_sum_const_dR(3*ngrad, natom))
            ALLOCATE (dP_i_dRj(3, natom, natom))
            ALLOCATE (dP_i_dRi(3, natom))

            DO iatom = 1, natom-1
               DO jatom = iatom+1, natom
                  r = particle_set(iatom)%r
                  r1 = particle_set(jatom)%r
                  DO i = 1, 3
                     r(i) = MODULO(r(i), cell%hmat(i, i))-cell%hmat(i, i)/2._dp
                     r1(i) = MODULO(r1(i), cell%hmat(i, i))-cell%hmat(i, i)/2._dp
                  END DO
                  dist_vec_i = (r-r1)-ANINT((r-r1)/cell_v)*cell_v
                  R12(iatom, jatom) = SQRT(DOT_PRODUCT(dist_vec_i, dist_vec_i))
                  R12(jatom, iatom) = R12(iatom, jatom)
               END DO
            END DO

            DO k = bo(1, 1), bo(2, 1)
               DO j = bo(1, 2), bo(2, 2)
                  DO i = bo(1, 3), bo(2, 3)
                     IF (dft_control%qs_control%becke_control%cavity_confine) THEN
                        IF (dft_control%qs_control%becke_control%cavity%pw%cr3d(k, j, i) < eps_cavity) CYCLE
                     END IF
                     grid_p(1) = k*dr(1)+shift(1)
                     grid_p(2) = j*dr(2)+shift(2)
                     grid_p(3) = i*dr(3)+shift(3)
                     ind = (/k, j, i/)
                     IF (dft_control%qs_control%becke_control%confine .OR. &
                         dft_control%qs_control%becke_control%dynamic_confine) THEN
                        IF (ind(dir) .LE. lb_index .OR. ind(dir) .GE. ub_index) CYCLE
                     END IF
                     nskipped = 0
                     d_sum_Pm_dR = 0.0_dp
                     d_sum_const_dR = 0.0_dp
                     P_i = 1.0_dp
                     dP_i_dRi = 0.0_dp
                     DO iatom = 1, natom
                        r = particle_set(iatom)%r
                        DO ip = 1, 3
                           r(ip) = MODULO(r(ip), cell%hmat(ip, ip))-cell%hmat(ip, ip)/2._dp
                        END DO
                        dist_vec_i = (r-grid_p)-ANINT((r-grid_p)/cell_v)*cell_v
                        dist1 = SQRT(DOT_PRODUCT(dist_vec_i, dist_vec_i))
                        IF (dist1 .LE. cutoffs(iatom)) THEN
                           IF (dist1 .LE. th) dist1 = th
                           dr_i_dR(:) = dist_vec_i(:)/dist1
                           DO jatom = 1, natom
                              IF (jatom .NE. iatom) THEN
                                 r1 = particle_set(jatom)%r
                                 DO ip = 1, 3
                                    r1(ip) = MODULO(r1(ip), cell%hmat(ip, ip))-cell%hmat(ip, ip)/2._dp
                                 END DO
                                 dist_vec_j = (r1-grid_p)-ANINT((r1-grid_p)/cell_v)*cell_v
                                 dr1_r2 = (r-r1)-ANINT((r-r1)/cell_v)*cell_v
                                 dist2 = SQRT(DOT_PRODUCT(dist_vec_j, dist_vec_j))
                                 IF (dist2 .LE. th) dist2 = th
                                 tmp_const = (R12(iatom, jatom)**3)
                                 dr_ij_dR(:) = dr1_r2(:)/tmp_const

                                 !derivativ w.r.t. Rj
                                 dr_j_dR = dist_vec_j(:)/dist2

                                 dmy_dR_j(:) = -(dr_j_dR(:)/R12(iatom, jatom)-(dist1-dist2)*dr_ij_dR(:))

                                 !derivativ w.r.t. Ri
                                 dmy_dR_i(:) = dr_i_dR(:)/R12(iatom, jatom)-(dist1-dist2)*dr_ij_dR(:)

                                 my1 = (dist1-dist2)/R12(iatom, jatom)

                                 IF (dft_control%qs_control%becke_control%adjust) THEN
                                    my1_homo = my1
                                    my1 = my1+ &
                                          dft_control%qs_control%becke_control%aij(iatom, jatom)*(1.0_dp-my1**2)
                                 END IF

                                 polyn = 1.5_dp*my1-0.5_dp*my1**3

                                 dpolyn = 1.5_dp-1.5_dp*my1**2

                                 tmp_const = (1.5_dp**2)*(dpolyn-dpolyn*(polyn**2))- &
                                             (1.5_dp)*((1.5_dp*polyn-0.5*(polyn**3))**2)* &
                                             (1.5_dp*dpolyn-1.5_dp*dpolyn*(polyn**2))

                                 ds_dR(:, iatom, jatom) = -0.5_dp*tmp_const*dmy_dR_i(:)
                                 ds_dR(:, jatom, iatom) = -0.5_dp*tmp_const*dmy_dR_j(:)
                                 IF (dft_control%qs_control%becke_control%adjust) THEN
                                    tmp_const = 1.0_dp-2.0_dp*my1_homo* &
                                                dft_control%qs_control%becke_control%aij(iatom, jatom)
                                    ds_dR(:, iatom, jatom) = ds_dR(:, iatom, jatom)*tmp_const
                                    ds_dR(:, jatom, iatom) = ds_dR(:, jatom, iatom)*tmp_const
                                 END IF
                                 f_3 = 1.5_dp*polyn-0.5_dp*polyn**3
                                 f_3 = 1.5_dp*f_3-0.5_dp*f_3**3
                                 s_my_ij(iatom, jatom) = 0.5_dp*(1-f_3)
                                 P_i(iatom) = P_i(iatom)*s_my_ij(iatom, jatom)
                                 IF (ABS(s_my_ij(iatom, jatom)) .LE. th) &
                                    s_my_ij(iatom, jatom) = s_my_ij(iatom, jatom)+th
                              END IF
                           END DO

                           DO jatom = 1, natom
                              IF (iatom .NE. jatom) THEN
                                 IF (ABS(s_my_ij(iatom, jatom)) .LE. th) &
                                    s_my_ij(iatom, jatom) = s_my_ij(iatom, jatom)+th
                                 dP_i_dRi(:, iatom) = dP_i_dRi(:, iatom)+ &
                                                      P_i(iatom)/s_my_ij(iatom, jatom)*ds_dR(:, iatom, jatom)
                                 dP_i_dRj(:, iatom, jatom) = P_i(iatom)/s_my_ij(iatom, jatom)*ds_dR(:, jatom, iatom)
                              END IF
                           END DO

                           d_sum_Pm_dR(:, iatom) = d_sum_Pm_dR(:, iatom)+dP_i_dRi(:, iatom)

                           IF (is_constraint(iatom)) &
                              d_sum_const_dR(1:3, iatom) = d_sum_const_dR(1:3, iatom)+dP_i_dRi(:, iatom)* &
                                                           coefficients(iatom)
                           IF (ngrad == 2) THEN
                              IF (is_combined(iatom)) &
                                 d_sum_const_dR(4:6, iatom) = d_sum_const_dR(4:6, iatom)+dP_i_dRi(:, iatom)
                           END IF
                           DO jatom = 1, natom
                              IF (iatom .NE. jatom) THEN
                                 d_sum_Pm_dR(:, jatom) = d_sum_Pm_dR(:, jatom)+dP_i_dRj(:, iatom, jatom)

                                 IF (is_constraint(iatom)) &
                                    d_sum_const_dR(:, jatom) = d_sum_const_dR(:, jatom)+dP_i_dRj(:, iatom, jatom)* &
                                                               coefficients(iatom)
                                 IF (ngrad == 2) THEN
                                    IF (is_combined(iatom)) &
                                       d_sum_const_dR(4:6, jatom) = d_sum_const_dR(4:6, jatom)+ &
                                                                    dP_i_dRj(:, iatom, jatom)
                                 END IF
                              END IF
                           END DO

                        ELSE
                           P_i(iatom) = 0.0_dp
                           IF (dft_control%qs_control%becke_control%should_skip) THEN
                              IF (is_constraint(iatom)) nskipped = nskipped+1
                              IF (nskipped == dft_control%qs_control%becke_control%natoms) EXIT
                           END IF
                        END IF
                     END DO !iatom
                     IF (nskipped == dft_control%qs_control%becke_control%natoms) CYCLE

                     sum_Pm = 0.0_dp
                     DO ip = 1, natom
                        sum_Pm = sum_Pm+P_i(ip)
                     END DO

                     sum_cell_f_constr = 0.0_dp
                     sum_cell_f_combined = 0.0_dp
                     DO ip = 1, dft_control%qs_control%becke_control%natoms
                        sum_cell_f_constr = sum_cell_f_constr+P_i(catom(ip))* &
                                            dft_control%qs_control%becke_control%coeff(ip)
                        IF (ngrad == 2) THEN
                           IF (is_combined(catom(ip))) &
                              sum_cell_f_combined = sum_cell_f_combined+P_i(catom(ip))
                        END IF
                     END DO

                     DO iatom = 1, natom
                        IF (ABS(sum_Pm) .EQ. 0.0_dp) THEN
                           dw_dR(:, iatom) = 0.0_dp
                        ELSE
                           dw_dR(1:3, iatom) = d_sum_const_dR(1:3, iatom)/sum_Pm- &
                                               sum_cell_f_constr*d_sum_Pm_dR(1:3, iatom)/(sum_Pm**2)
                           IF (ngrad == 2) THEN
                              dw_dR(4:6, iatom) = d_sum_const_dR(4:6, iatom)/sum_Pm- &
                                                  sum_cell_f_combined*d_sum_Pm_dR(1:3, iatom)/(sum_Pm**2)
                           END IF
                        END IF
                     END DO

                     DO iatom = 1, natom
                        DO ispin = 1, dft_control%nspins
                           DO ivar = 1, nvar
                              sign = 1.0_dp
                              IF (ivar == 2) THEN
                                 IF (ispin == 2) sign = -1.0_dp
                                 integrated(1+3*(ngrad-1):3+3*(ngrad-1), iatom) = &
                                    integrated(1+3*(ngrad-1):3+3*(ngrad-1), iatom)+ &
                                    sign*dw_dR(1+3*(ngrad-1):3+3*(ngrad-1), iatom)*rho_r(ispin)%pw%cr3d(k, j, i)*dvol
                              ELSE
                                 IF (dft_control%qs_control%becke_control%constraint_type == &
                                     cdft_magnetization_constraint .AND. ispin == 2) sign = -1.0_dp
                                 integrated(1:3, iatom) = integrated(1:3, iatom)+ &
                                                          sign*dw_dR(1:3, iatom)*rho_r(ispin)%pw%cr3d(k, j, i)*dvol
                              END IF
                           END DO ! ivar
                        END DO ! ispin
                     END DO ! iatom
                  END DO ! i
               END DO ! j
            END DO ! k
            DEALLOCATE (catom)
            DEALLOCATE (R12)
            DEALLOCATE (s_my_ij)
            DEALLOCATE (ds_dR)
            DEALLOCATE (P_i)
            DEALLOCATE (dw_dR)
            DEALLOCATE (is_constraint)
            DEALLOCATE (d_sum_Pm_dR)
            DEALLOCATE (d_sum_const_dR)
            DEALLOCATE (dP_i_dRi)
            DEALLOCATE (dP_i_dRj)
            DEALLOCATE (coefficients)
            IF (ALLOCATED(is_combined)) DEALLOCATE (is_combined)
         END IF ! Gradient integration

         CALL mp_sum(integrated, para_env%group)

         ! Update force only on master process. Otherwise force due to constraint becomes multiplied
         ! by the number of processes when final force%rho_elec is constructed in qs_force
         ! by mp_summing [the final integrated(:,:) is distributed on all processors]
         IF (para_env%mepos == para_env%source) THEN
            DO iatom = 1, natom
               ikind = kind_of(iatom)
               i = atom_of_kind(iatom)
               IF (ngrad == 2) THEN
                  ! With combined constraint donor/acceptor, gradient due to charge constraint is in integrated(1:3,:)
                  ! and due to spin constraint in integrated(4:6,:)
                  force(ikind)%rho_elec(:, i) = force(ikind)%rho_elec(:, i) &
                                                +integrated(1:3, iatom)*strength(1)
                  force(ikind)%rho_elec(:, i) = &
                     force(ikind)%rho_elec(:, i)+ &
                     integrated(4:6, iatom)*strength(2)
               ELSE
                  ! With combined constraint all, sum of gradients from both constraints are in integrated(1:3,:)
                  force(ikind)%rho_elec(:, i) = force(ikind)%rho_elec(:, i) &
                                                +integrated(1:3, iatom)*SUM(strength)
               END IF
            END DO
         END IF
         DEALLOCATE (strength)
         DEALLOCATE (integrated)
         DEALLOCATE (atom_of_kind)
         DEALLOCATE (kind_of)
      END IF

      CALL timestop(handle)

   END SUBROUTINE becke_force

! **************************************************************************************************
!> \brief calculates a Gaussian Hirshfeld constraint
!> \param qs_env the qs_env where to build the constraint
!> \param cdft_control the container for constraint related structures
!> \param calc_pot if the constraint potential should be rebuilt or just integrated
!> \param calculate_forces logical if potential has to be calculated or only_energy
!> \author  Nico Holmberg (01.2016)
! **************************************************************************************************
   SUBROUTINE hirshfeld_constraint(qs_env, cdft_control, calc_pot, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cdft_control_type)                            :: cdft_control
      LOGICAL                                            :: calc_pot, calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'hirshfeld_constraint', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, handle, i, iat, iatom, iex, &
                                                            ikind, ithread, iw, j, natom, nkind, &
                                                            npme, nthread, numexp, unit_nr
      INTEGER(KIND=int_8)                                :: subpatch_pattern
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: catom
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores, stride
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: is_constraint
      REAL(kind=dp)                                      :: alpha, coef, dE, dvol, eps_rho_rspace, &
                                                            ra(3), tnfun, tnfun2, zeff
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: coefficients
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cube_info_type)                               :: cube_info
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hirshfeld_type), POINTER                      :: hirshfeld_env
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: tmp
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(pw_p_type), POINTER                           :: fnorm
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(realspace_grid_desc_type), POINTER            :: auxbas_rs_desc
      TYPE(realspace_grid_type), POINTER                 :: rs_rho_all, rs_rho_constr
      TYPE(section_vals_type), POINTER                   :: print_section

      NULLIFY (atom_list, cores, stride, pab, atomic_kind, atomic_kind_set, cell, &
               logger, para_env, dft_control, hirshfeld_env, particles, particle_set, &
               pw_env, fnorm, rho_r, auxbas_pw_pool, energy, rho, qs_kind_set, &
               subsys, auxbas_rs_desc, rs_rho_all, rs_rho_constr, print_section)
      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      particle_set=particle_set, &
                      rho=rho, &
                      natom=natom, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      qs_kind_set=qs_kind_set, &
                      subsys=subsys)

      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(qs_kind_set))

      print_section => section_vals_get_subs_vals(qs_env%input, &
                                                  "DFT%QS%CDFT%HIRSHFELD_CONSTRAINT")
      iw = cp_print_key_unit_nr(logger, print_section, "PROGRAM_RUN_INFO", &
                                extension=".cdftLog")
      unit_nr = cp_print_key_unit_nr(logger, print_section, "PROGRAM_RUN_INFO", &
                                     extension=".cube", log_filename=.FALSE.)
      CALL qs_rho_get(rho, rho_r=rho_r)

      IF (dft_control%qs_control%cdft) THEN
         dvol = cdft_control%weight%pw%pw_grid%dvol
         ! Build weight function
         IF (calc_pot) THEN
            cdft_control%weight%pw%cr3d = 0.0_dp

            ALLOCATE (catom(cdft_control%hirshfeld_control%natoms))
            ALLOCATE (coefficients(natom))
            ALLOCATE (is_constraint(natom))
            is_constraint = .FALSE.
            coefficients(:) = 0.0_dp

            DO i = 1, cdft_control%hirshfeld_control%natoms
               catom(i) = cdft_control%hirshfeld_control%atoms(i)
               is_constraint(catom(i)) = .TRUE.
               coefficients(catom(i)) = cdft_control%hirshfeld_control%coeff(i)
            END DO

            hirshfeld_env => cdft_control%hirshfeld_control%hirshfeld_env
            CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, pw_env=pw_env)
            nkind = SIZE(qs_kind_set)

            ! Setup Hirshfeld
            IF (.NOT. ASSOCIATED(hirshfeld_env%kind_shape_fn)) THEN
               CALL create_shape_function(hirshfeld_env, qs_kind_set, atomic_kind_set)
               CPASSERT(.NOT. ASSOCIATED(hirshfeld_env%charges))
               ALLOCATE (hirshfeld_env%charges(natom))
               SELECT CASE (hirshfeld_env%ref_charge)
               CASE (ref_charge_atomic)
                  DO ikind = 1, nkind
                     CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
                     atomic_kind => atomic_kind_set(ikind)
                     CALL get_atomic_kind(atomic_kind, atom_list=atom_list)
                     DO iat = 1, SIZE(atom_list)
                        i = atom_list(iat)
                        hirshfeld_env%charges(i) = zeff
                     END DO
                  END DO
               CASE (ref_charge_mulliken)
                  CALL cp_abort(__LOCATION__, &
                                "Mulliken reference charge NYI for Hirshfeld partitioning.")
               CASE DEFAULT
                  CALL cp_abort(__LOCATION__, &
                                "Unknown type of reference charge for Hirshfeld partitioning.")
               END SELECT
            END IF
            ! Calculate densities
            ! rs_rho_all: Promolecular density i.e. superposition of isolated Gaussian densities
            ! rs_rho_constr: Sum of isolated Gaussian densities over constraint atoms multiplied by coeff
            CALL pw_env_get(pw_env, auxbas_rs_desc=auxbas_rs_desc, auxbas_rs_grid=rs_rho_all, &
                            auxbas_pw_pool=auxbas_pw_pool)
            cube_info = pw_env%cube_info(1)
            CALL rs_grid_retain(rs_rho_all)
            CALL rs_grid_zero(rs_rho_all)
            CALL rs_grid_create(rs_rho_constr, auxbas_rs_desc)
            CALL rs_grid_zero(rs_rho_constr)
            eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
            ALLOCATE (pab(1, 1))
            nthread = 1
            ithread = 0

            DO ikind = 1, SIZE(atomic_kind_set)
               numexp = hirshfeld_env%kind_shape_fn(ikind)%numexp
               IF (numexp <= 0) CYCLE
               CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
               ALLOCATE (cores(natom))

               DO iex = 1, numexp
                  alpha = hirshfeld_env%kind_shape_fn(ikind)%zet(iex)
                  coef = hirshfeld_env%kind_shape_fn(ikind)%coef(iex)
                  npme = 0
                  cores = 0
                  DO iatom = 1, natom
                     atom_a = atom_list(iatom)
                     ra(:) = pbc(particle_set(atom_a)%r, cell)
                     IF (rs_rho_all%desc%parallel .AND. .NOT. rs_rho_all%desc%distributed) THEN
                        ! replicated realspace grid, split the atoms up between procs
                        IF (MODULO(iatom, rs_rho_all%desc%group_size) == rs_rho_all%desc%my_pos) THEN
                           npme = npme+1
                           cores(npme) = iatom
                        ENDIF
                     ELSE
                        npme = npme+1
                        cores(npme) = iatom
                     ENDIF
                  END DO
                  DO j = 1, npme
                     iatom = cores(j)
                     atom_a = atom_list(iatom)
                     pab(1, 1) = hirshfeld_env%charges(atom_a)*coef
                     ra(:) = pbc(particle_set(atom_a)%r, cell)
                     subpatch_pattern = 0
                     CALL collocate_pgf_product_rspace(0, alpha, 0, 0, 0.0_dp, 0, ra, &
                                                       (/0.0_dp, 0.0_dp, 0.0_dp/), 0.0_dp, 1.0_dp, pab, 0, 0, &
                                                       rs_rho_all, cell, cube_info, eps_rho_rspace, &
                                                       ga_gb_function=FUNC_AB, ithread=ithread, use_subpatch=.TRUE., &
                                                       subpatch_pattern=subpatch_pattern, lmax_global=0)
                     IF (is_constraint(atom_a)) &
                        CALL collocate_pgf_product_rspace(0, alpha, 0, 0, 0.0_dp, 0, ra, &
                                                          (/0.0_dp, 0.0_dp, 0.0_dp/), 0.0_dp, coefficients(atom_a), &
                                                          pab, 0, 0, rs_rho_constr, cell, cube_info, eps_rho_rspace, &
                                                          ga_gb_function=FUNC_AB, ithread=ithread, use_subpatch=.TRUE., &
                                                          subpatch_pattern=subpatch_pattern, lmax_global=0)
                  END DO
               END DO
               DEALLOCATE (cores)
            END DO
            DEALLOCATE (pab)
            DEALLOCATE (catom)
            DEALLOCATE (is_constraint)
            DEALLOCATE (coefficients)
            CALL get_hirshfeld_info(hirshfeld_env, fnorm=fnorm)
            IF (ASSOCIATED(fnorm)) THEN
               CALL pw_release(fnorm%pw)
               DEALLOCATE (fnorm)
            ENDIF
            ALLOCATE (fnorm)

            CALL pw_pool_create_pw(auxbas_pw_pool, fnorm%pw, use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL pw_pool_create_pw(auxbas_pw_pool, tmp%pw, use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL set_hirshfeld_info(hirshfeld_env, fnorm=fnorm)
            CALL rs_pw_transfer(rs_rho_all, fnorm%pw, rs2pw)
            CALL rs_pw_transfer(rs_rho_constr, tmp%pw, rs2pw)
            IF (debug_this_module) THEN
               ALLOCATE (stride(3))
               stride = (/2, 2, 2/)
               CALL qs_subsys_get(subsys, particles=particles)
               CALL cp_pw_to_cube(tmp%pw, unit_nr, "RHO_CONSTR", &
                                  particles=particles, &
                                  stride=stride)
               DEALLOCATE (stride)
            END IF
            tnfun2 = pw_integrate_function(tmp%pw)
            CALL rs_grid_release(rs_rho_constr)
            CALL hfun_scale(cdft_control%weight%pw%cr3d, tmp%pw%cr3d, &
                            hirshfeld_env%fnorm%pw%cr3d, divide)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, tmp%pw)
            CALL rs_grid_release(rs_rho_all)
            tnfun = pw_integrate_function(hirshfeld_env%fnorm%pw)
         END IF

         dE = 0.0_dp
         DO i = 1, dft_control%nspins
            dE = dE+accurate_sum(cdft_control%weight%pw%cr3d*rho_r(i)%pw%cr3d)*dvol
         END DO

         CALL get_qs_env(qs_env, energy=energy)
         CALL mp_sum(dE, para_env%group)

         cdft_control%value = dE
         energy%cdft = (cdft_control%value(1)-cdft_control%target(1))*cdft_control%strength(1)

         IF (iw > 0) THEN
            WRITE (iw, '(/,T3,A,T60)') &
               '------------------- Hirshfeld constraint information -------------------'
            WRITE (iw, '(T3,A,T60,(3X,F18.12))') &
               'Target value of constraint        :', cdft_control%target
            WRITE (iw, '(T3,A,T60,(3X,F18.12))') &
               'Current value of constraint       :', cdft_control%value
            WRITE (iw, '(T3,A,T60,(3X,ES12.2))') &
               'Deviation from target             :', ABS(cdft_control%value-cdft_control%target)
            WRITE (iw, '(T3,A,T54,(3X,F18.12))') &
               'Strength of constraint            :', cdft_control%strength
            WRITE (iw, '(T3,A)') &
               '------------------------------------------------------------------------'
         END IF
         CALL cp_print_key_finished_output(iw, logger, print_section, "PROGRAM_RUN_INFO")
         CALL cp_print_key_finished_output(unit_nr, logger, print_section, "PROGRAM_RUN_INFO")

         IF (calculate_forces) &
            CPABORT("Forces with Hirshfeld constraint NYI.")
         !CALL hirshfeld_force(qs_env,cdft_control)
      END IF
      CALL timestop(handle)

   END SUBROUTINE hirshfeld_constraint

! **************************************************************************************************
!> \brief calculates forces due to a Gaussian Hirshfeld constraint
!> \param qs_env the qs_env where to build the constraint
!> \param cdft_control the container for constraint related structures
!> \author  Nico Holmberg (01.2016)
! **************************************************************************************************
   SUBROUTINE hirshfeld_force(qs_env, cdft_control)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cdft_control_type)                            :: cdft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'hirshfeld_force', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, handle, i, iatom, iex, ikind, &
                                                            is, ithread, j, katom, natom, npme, &
                                                            nthread, numexp
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, catom, kind_of
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores
      LOGICAL                                            :: rho_r_valid
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: is_acceptor, is_donor
      REAL(kind=dp)                                      :: alpha, coef, dvol, eps_rho_rspace, &
                                                            force_a(3), force_b(3), ra(3)
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: coefficients
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fderiv, fval, hab, pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hirshfeld_type), POINTER                      :: hirshfeld_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: acceptor_weight, donor_weight, rhonorm, &
                                                            tmp, weight
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(realspace_grid_desc_type), POINTER            :: auxbas_rs_desc
      TYPE(realspace_grid_type), POINTER                 :: acceptor_weight_rs, donor_weight_rs, &
                                                            weight_rs

      NULLIFY (atom_list, cores, pab, hab, fderiv, fval, atomic_kind_set, &
               cell, para_env, dft_control, hirshfeld_env, particle_set, &
               pw_env, rho_r, auxbas_pw_pool, force, rho, auxbas_rs_desc, &
               weight_rs, acceptor_weight_rs, donor_weight_rs)
      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      particle_set=particle_set, &
                      rho=rho, &
                      force=force, &
                      natom=natom, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      atomic_kind_set=atomic_kind_set, &
                      pw_env=pw_env)

      IF (dft_control%qs_control%cdft) THEN
         dvol = cdft_control%weight%pw%pw_grid%dvol
         ALLOCATE (catom(cdft_control%hirshfeld_control%natoms))
         ALLOCATE (coefficients(natom))
         ALLOCATE (is_acceptor(natom))
         ALLOCATE (is_donor(natom))
         is_acceptor = .FALSE.
         is_donor = .FALSE.
         coefficients(:) = 0.0_dp

         DO i = 1, cdft_control%hirshfeld_control%natoms
            catom(i) = cdft_control%hirshfeld_control%atoms(i)
            coefficients(catom(i)) = cdft_control%hirshfeld_control%coeff(i)
            IF (coefficients(catom(i)) < 0.0_dp) THEN
               is_donor(catom(i)) = .TRUE.
            ELSE
               is_acceptor(catom(i)) = .TRUE.
            END IF
         END DO

         hirshfeld_env => cdft_control%hirshfeld_control%hirshfeld_env
         ! Initialize working arrays
         CALL qs_rho_get(rho, rho_r=rho_r, rho_r_valid=rho_r_valid)
         CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                         auxbas_rs_desc=auxbas_rs_desc, auxbas_rs_grid=weight_rs)
         CALL rs_grid_retain(weight_rs)
         CALL pw_pool_create_pw(auxbas_pw_pool, weight%pw, use_data=REALDATA3D)
         CALL pw_pool_create_pw(auxbas_pw_pool, donor_weight%pw, use_data=REALDATA3D)
         CALL pw_pool_create_pw(auxbas_pw_pool, acceptor_weight%pw, use_data=REALDATA3D)
         CALL pw_pool_create_pw(auxbas_pw_pool, tmp%pw, use_data=REALDATA3D)
         ! Here, we assume that all donor atoms have coeff +1 and acceptors -1
         tmp%pw%cr3d = 1.0_dp
         CALL pw_copy(cdft_control%weight%pw, weight%pw)
         CALL pw_scale(weight%pw, -1.0_dp)
         CALL pw_copy(weight%pw, donor_weight%pw)
         CALL pw_copy(weight%pw, acceptor_weight%pw)
         CALL pw_axpy(tmp%pw, donor_weight%pw, 1.0_dp)
         CALL pw_axpy(tmp%pw, acceptor_weight%pw, -1.0_dp)

         CALL rs_grid_create(donor_weight_rs, auxbas_rs_desc)
         CALL rs_grid_create(acceptor_weight_rs, auxbas_rs_desc)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhonorm%pw, use_data=REALDATA3D)
         ALLOCATE (fderiv(3, natom))
         fderiv = 0.0_dp
         ALLOCATE (fval(3, natom))
         fval = 0.0_dp
         eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
         nthread = 1
         ithread = 0
         ALLOCATE (hab(1, 1), pab(1, 1))
         ! loop over spins
         DO is = 1, SIZE(rho_r)
            IF (rho_r_valid) THEN
               CALL hfun_scale(rhonorm%pw%cr3d, rho_r(is)%pw%cr3d, &
                               hirshfeld_env%fnorm%pw%cr3d, divide)
            ELSE
               CPABORT("We need rho in real space")
            END IF
            ! Weight for non-constraint atoms
            CALL hfun_scale(tmp%pw%cr3d, rhonorm%pw%cr3d, weight%pw%cr3d, multiply)
            CALL rs_pw_transfer(weight_rs, tmp%pw, pw2rs)
            ! Weight for donor atoms
            CALL hfun_scale(tmp%pw%cr3d, rhonorm%pw%cr3d, donor_weight%pw%cr3d, multiply)
            CALL rs_pw_transfer(donor_weight_rs, tmp%pw, pw2rs)
            ! Weight for acceptor atoms
            CALL hfun_scale(tmp%pw%cr3d, rhonorm%pw%cr3d, acceptor_weight%pw%cr3d, multiply)
            CALL rs_pw_transfer(acceptor_weight_rs, tmp%pw, pw2rs)
            ! Calculate gradients
            DO ikind = 1, SIZE(atomic_kind_set)
               numexp = hirshfeld_env%kind_shape_fn(ikind)%numexp
               IF (numexp <= 0) CYCLE
               CALL get_atomic_kind(atomic_kind_set(ikind), natom=katom, atom_list=atom_list)
               ALLOCATE (cores(katom))

               DO iex = 1, numexp
                  alpha = hirshfeld_env%kind_shape_fn(ikind)%zet(iex)
                  coef = hirshfeld_env%kind_shape_fn(ikind)%coef(iex)
                  npme = 0
                  cores = 0
                  DO iatom = 1, katom
                     atom_a = atom_list(iatom)
                     ra(:) = pbc(particle_set(atom_a)%r, cell)
                     IF (weight_rs%desc%parallel .AND. .NOT. weight_rs%desc%distributed) THEN
                        ! replicated realspace grid, split the atoms up between procs
                        IF (MODULO(iatom, weight_rs%desc%group_size) == weight_rs%desc%my_pos) THEN
                           npme = npme+1
                           cores(npme) = iatom
                        ENDIF
                     ELSE
                        npme = npme+1
                        cores(npme) = iatom
                     ENDIF
                  END DO

                  DO j = 1, npme
                     iatom = cores(j)
                     atom_a = atom_list(iatom)
                     ra(:) = pbc(particle_set(atom_a)%r, cell)
                     pab(1, 1) = coef*hirshfeld_env%charges(atom_a)
                     hab(1, 1) = 0.0_dp
                     force_a(:) = 0.0_dp
                     force_b(:) = 0.0_dp
                     IF (is_donor(atom_a)) THEN
                        CALL integrate_pgf_product_rspace(0, alpha, 0, &
                                                          0, 0.0_dp, 0, ra, (/0.0_dp, 0.0_dp, 0.0_dp/), 0.0_dp, &
                                                          donor_weight_rs, cell, pw_env%cube_info(1), hab, pab=pab, &
                                                          o1=0, o2=0, eps_gvg_rspace=eps_rho_rspace, &
                                                          calculate_forces=.TRUE., force_a=force_a, force_b=force_b, &
                                                          use_virial=.FALSE., use_subpatch=.TRUE., &
                                                          subpatch_pattern=0_int_8)
                     ELSE IF (is_acceptor(atom_a)) THEN
                        CALL integrate_pgf_product_rspace(0, alpha, 0, &
                                                          0, 0.0_dp, 0, ra, (/0.0_dp, 0.0_dp, 0.0_dp/), 0.0_dp, &
                                                          acceptor_weight_rs, cell, pw_env%cube_info(1), hab, pab=pab, &
                                                          o1=0, o2=0, eps_gvg_rspace=eps_rho_rspace, &
                                                          calculate_forces=.TRUE., force_a=force_a, force_b=force_b, &
                                                          use_virial=.FALSE., use_subpatch=.TRUE., &
                                                          subpatch_pattern=0_int_8)
                     ELSE
                        CALL integrate_pgf_product_rspace(0, alpha, 0, &
                                                          0, 0.0_dp, 0, ra, (/0.0_dp, 0.0_dp, 0.0_dp/), 0.0_dp, &
                                                          weight_rs, cell, pw_env%cube_info(1), hab, pab=pab, &
                                                          o1=0, o2=0, eps_gvg_rspace=eps_rho_rspace, &
                                                          calculate_forces=.TRUE., force_a=force_a, force_b=force_b, &
                                                          use_virial=.FALSE., use_subpatch=.TRUE., &
                                                          subpatch_pattern=0_int_8)
                     END IF
                     fderiv(:, atom_a) = fderiv(:, atom_a)+force_a(:)*dvol
                     fval(:, atom_a) = fval(:, atom_a)+hab(1, 1)*dvol*coef*hirshfeld_env%charges(atom_a)
                  END DO
               END DO
               DEALLOCATE (cores)
            END DO
         END DO
         ! Cleanup working arrays
         CALL pw_pool_give_back_pw(auxbas_pw_pool, weight%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhonorm%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, tmp%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, acceptor_weight%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, donor_weight%pw)
         CALL rs_grid_release(weight_rs)
         CALL rs_grid_release(donor_weight_rs)
         CALL rs_grid_release(acceptor_weight_rs)
         DEALLOCATE (hab, pab)
         DEALLOCATE (catom)
         DEALLOCATE (is_acceptor, is_donor)
         DEALLOCATE (coefficients)
         ! Update forces. This is actually the gradient, so the sign needs to be positive!
         CALL mp_sum(fderiv, para_env%group)
         CALL mp_sum(fval, para_env%group)
         ALLOCATE (atom_of_kind(natom))
         ALLOCATE (kind_of(natom))
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  atom_of_kind=atom_of_kind, &
                                  kind_of=kind_of)
         ! Update only on master proc, see becke_force
         IF (para_env%mepos == para_env%source) THEN
            DO iatom = 1, natom
               ikind = kind_of(iatom)
               i = atom_of_kind(iatom)
               force(ikind)%rho_elec(:, i) = &
                  force(ikind)%rho_elec(:, i)+fderiv(:, iatom)*cdft_control%strength
            END DO
         END IF
         DEALLOCATE (fderiv)
         DEALLOCATE (kind_of)
         DEALLOCATE (atom_of_kind)
      END IF
      CALL timestop(handle)

   END SUBROUTINE hirshfeld_force

! **************************************************************************************************
!> \brief Calculate fout = fun1/fun2 or fout = fun1*fun2
!> \param fout the output 3D potential
!> \param fun1 the first input 3D potential
!> \param fun2 the second input 3D potential
!> \param mode decides whether to divide or multiply the input potentials
! **************************************************************************************************
   SUBROUTINE hfun_scale(fout, fun1, fun2, mode)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: fout
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: fun1, fun2
      INTEGER, INTENT(IN)                                :: mode

      CHARACTER(len=*), PARAMETER :: routineN = 'hfun_scale', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: small = 1.0e-12_dp

      INTEGER                                            :: i1, i2, i3, n1, n2, n3

      n1 = SIZE(fout, 1)
      n2 = SIZE(fout, 2)
      n3 = SIZE(fout, 3)
      CPASSERT(n1 == SIZE(fun1, 1))
      CPASSERT(n2 == SIZE(fun1, 2))
      CPASSERT(n3 == SIZE(fun1, 3))
      CPASSERT(n1 == SIZE(fun2, 1))
      CPASSERT(n2 == SIZE(fun2, 2))
      CPASSERT(n3 == SIZE(fun2, 3))

      SELECT CASE (mode)
      CASE (divide)
         DO i3 = 1, n3
            DO i2 = 1, n2
               DO i1 = 1, n1
                  IF (fun2(i1, i2, i3) > small) THEN
                     fout(i1, i2, i3) = fun1(i1, i2, i3)/fun2(i1, i2, i3)
                  ELSE
                     fout(i1, i2, i3) = 0.0_dp
                  END IF
               END DO
            END DO
         END DO
      CASE (multiply)
         DO i3 = 1, n3
            DO i2 = 1, n2
               DO i1 = 1, n1
                  fout(i1, i2, i3) = fun1(i1, i2, i3)*fun2(i1, i2, i3)
               END DO
            END DO
         END DO
      END SELECT

   END SUBROUTINE hfun_scale

! **************************************************************************************************
!> \brief Determine confinement bounds along confinement dir (hardcoded to be z)
!>        and optionally zero entries below a given threshold
!> \param fun input 3D potential (real space)
!> \param th threshold for screening values
!> \param just_bounds if the bounds should be computed without zeroing values
!> \param bounds the confinement bounds: fun is nonzero only between these values along 3rd dimension
! **************************************************************************************************
   SUBROUTINE hfun_zero(fun, th, just_bounds, bounds)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: fun
      REAL(KIND=dp), INTENT(IN)                          :: th
      LOGICAL                                            :: just_bounds
      INTEGER, OPTIONAL                                  :: bounds(2)

      CHARACTER(len=*), PARAMETER :: routineN = 'hfun_zero', routineP = moduleN//':'//routineN

      INTEGER                                            :: i1, i2, i3, lb, n1, n2, n3, nzeroed, ub
      LOGICAL                                            :: lb_final, ub_final

      n1 = SIZE(fun, 1)
      n2 = SIZE(fun, 2)
      n3 = SIZE(fun, 3)
      IF (just_bounds) THEN
         CPASSERT(PRESENT(bounds))
         lb = 1
         lb_final = .FALSE.
         ub_final = .FALSE.
      END IF

      DO i3 = 1, n3
         IF (just_bounds) nzeroed = 0
         DO i2 = 1, n2
            DO i1 = 1, n1
               IF (fun(i1, i2, i3) < th) THEN
                  IF (just_bounds) THEN
                     nzeroed = nzeroed+1
                  ELSE
                     fun(i1, i2, i3) = 0.0_dp
                  END IF
               ELSE
                  IF (just_bounds) EXIT
               END IF
            END DO
            IF (just_bounds .AND. nzeroed < n1) EXIT
         END DO
         IF (just_bounds) THEN
            IF (nzeroed == (n2*n1)) THEN
               IF (.NOT. lb_final) THEN
                  lb = i3
               ELSE IF (.NOT. ub_final) THEN
                  ub = i3
                  ub_final = .TRUE.
               END IF
            ELSE
               IF (.NOT. lb_final) lb_final = .TRUE.
               IF (ub_final) ub_final = .FALSE. ! Safeguard against "holes"
            END IF
         END IF
      END DO
      IF (just_bounds) THEN
         IF (.NOT. ub_final) ub = n3
         bounds(1) = lb
         bounds(2) = ub
         bounds = bounds-(n3/2)-1
      END IF

   END SUBROUTINE hfun_zero

END MODULE et_coupling

