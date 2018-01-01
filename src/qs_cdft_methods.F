!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Subroutines for building CDFT constraints
!> \par   History
!>                 separated from et_coupling [03.2017]
!> \author Nico Holmberg [03.2017]
! **************************************************************************************************
MODULE qs_cdft_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE cube_utils,                      ONLY: cube_info_type
   USE hirshfeld_methods,               ONLY: create_shape_function
   USE hirshfeld_types,                 ONLY: get_hirshfeld_info,&
                                              hirshfeld_type,&
                                              set_hirshfeld_info
   USE input_constants,                 ONLY: cdft_alpha_constraint,&
                                              cdft_beta_constraint,&
                                              cdft_charge_constraint,&
                                              cdft_magnetization_constraint,&
                                              ref_charge_atomic,&
                                              ref_charge_mulliken
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kahan_sum,                       ONLY: accurate_dot_product,&
                                              accurate_sum
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_sum
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integrate_function,&
                                              pw_scale
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_p_type,&
                                              pw_release
   USE qs_cdft_types,                   ONLY: becke_constraint_type,&
                                              becke_group_type,&
                                              cdft_control_type
   USE qs_cdft_utils,                   ONLY: becke_constraint_init,&
                                              hfun_scale
   USE qs_collocate_density,            ONLY: collocate_pgf_product_rspace
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integrate_potential,          ONLY: integrate_pgf_product_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_cdft_methods'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

! *** Public subroutines ***

   PUBLIC :: becke_constraint, hirshfeld_constraint

CONTAINS

! **************************************************************************************************
!> \brief Driver routine for calculating a Becke constraint
!> \param qs_env the qs_env where to build the constraint
!> \param calc_pot if the potential needs to be recalculated or just integrated
!> \param calculate_forces logical if potential has to be calculated or only_energy
!> \par   History
!>        Created 01.2007 [fschiff]
!>        Extended functionality 12/15-12/16 [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE becke_constraint(qs_env, calc_pot, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: calc_pot, calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'becke_constraint', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      IF (dft_control%qs_control%becke_constraint) THEN
         IF (calc_pot) THEN
            ! Initialize the Becke constraint environment
            CALL becke_constraint_init(qs_env)
            ! Calculate the Becke weight function and possibly the gradients
            CALL becke_constraint_low(qs_env)
         END IF
         ! Integrate the Becke constraint
         CALL becke_constraint_integrate(qs_env)
         ! Calculate gradients
         IF (calculate_forces) CALL becke_constraint_force(qs_env)
      END IF
      CALL timestop(handle)

   END SUBROUTINE becke_constraint

! **************************************************************************************************
!> \brief Low level routine to build a Becke weight function and its gradients
!> \param qs_env the qs_env where to build the constraint
!> \param just_gradients optional logical which determines if only the gradients should be calculated
!> \par   History
!>        Created 03.2017 [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE becke_constraint_low(qs_env, just_gradients)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, OPTIONAL                                  :: just_gradients

      CHARACTER(len=*), PARAMETER :: routineN = 'becke_constraint_low', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iatom, igroup, ind(3), ip, j, &
                                                            jatom, jp, k, natom, np(3), nskipped
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: catom
      INTEGER, DIMENSION(2, 3)                           :: bo, bo_conf
      LOGICAL                                            :: in_memory, my_just_gradients
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: is_constraint, skip_me
      LOGICAL, ALLOCATABLE, DIMENSION(:, :)              :: atom_in_group
      REAL(kind=dp)                                      :: dist1, dist2, dmyexp, dvol, eps_cavity, &
                                                            my1, my1_homo, myexp, sum_cell_f_all, &
                                                            th, tmp_const
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: cell_functions, ds_dR_i, ds_dR_j, &
                                                            sum_cell_f_group
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: d_sum_Pm_dR, dP_i_dRi
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dP_i_dRj
      REAL(kind=dp), DIMENSION(3)                        :: cell_v, dist_vec, dmy_dR_i, dmy_dR_j, &
                                                            dr, dr1_r2, dr_i_dR, dr_ij_dR, &
                                                            dr_j_dR, grid_p, r, r1, shift
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cutoffs
      TYPE(becke_constraint_type), POINTER               :: becke_control
      TYPE(becke_group_type), DIMENSION(:), POINTER      :: group
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: charge

      NULLIFY (cutoffs, cell, dft_control, particle_set, group, charge)

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      particle_set=particle_set, &
                      natom=natom, &
                      dft_control=dft_control)
      becke_control => dft_control%qs_control%becke_control
      group => becke_control%group
      cutoffs => becke_control%cutoffs
      IF (becke_control%atomic_charges) &
         charge => becke_control%charge
      in_memory = .FALSE.
      IF (becke_control%save_pot) THEN
         in_memory = becke_control%in_memory
      END IF
      eps_cavity = becke_control%eps_cavity
      ! Decide if only gradients need to be calculated
      my_just_gradients = .FALSE.
      IF (PRESENT(just_gradients)) my_just_gradients = just_gradients
      IF (my_just_gradients) THEN
         in_memory = .TRUE.
         !  Pairwise distances need to be recalculated
         IF (becke_control%vector_buffer%store_vectors) THEN
            ALLOCATE (becke_control%vector_buffer%distances(natom))
            ALLOCATE (becke_control%vector_buffer%distance_vecs(3, natom))
            IF (in_memory) ALLOCATE (becke_control%vector_buffer%pair_dist_vecs(3, natom, natom))
            ALLOCATE (becke_control%vector_buffer%position_vecs(3, natom))
         END IF
         ALLOCATE (becke_control%vector_buffer%R12(natom, natom))
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
               IF (becke_control%vector_buffer%store_vectors) THEN
                  becke_control%vector_buffer%position_vecs(:, iatom) = r(:)
                  IF (iatom == 1 .AND. jatom == natom) becke_control%vector_buffer%position_vecs(:, jatom) = r1(:)
                  IF (in_memory) THEN
                     becke_control%vector_buffer%pair_dist_vecs(:, iatom, jatom) = dist_vec(:)
                     becke_control%vector_buffer%pair_dist_vecs(:, jatom, iatom) = -dist_vec(:)
                  END IF
               END IF
               becke_control%vector_buffer%R12(iatom, jatom) = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
               becke_control%vector_buffer%R12(jatom, iatom) = becke_control%vector_buffer%R12(iatom, jatom)
            END DO
         END DO
      END IF
      ALLOCATE (catom(becke_control%natoms))
      IF (becke_control%save_pot .OR. &
          becke_control%cavity_confine .OR. &
          becke_control%should_skip) THEN
         ALLOCATE (is_constraint(natom))
         is_constraint = .FALSE.
      END IF
      ! This boolean is needed to prevent calculation of atom pairs ji when the pair ij has
      ! already been calculated (data for pair ji is set using symmetry)
      ! With gradient precomputation, symmetry exploited for both weight function and gradients
      ALLOCATE (skip_me(natom))
      DO i = 1, becke_control%natoms
         catom(i) = becke_control%atoms(i)
         ! Notice that here is_constraint=.TRUE. also for dummy atoms to properly compute their Becke charges
         ! A subsequent check (atom_in_group) ensures that the gradients of these dummy atoms are correct
         IF (becke_control%save_pot .OR. &
             becke_control%cavity_confine .OR. &
             becke_control%should_skip) &
            is_constraint(catom(i)) = .TRUE.
      END DO
      bo = group(1)%weight%pw%pw_grid%bounds_local
      np = group(1)%weight%pw%pw_grid%npts
      dvol = group(1)%weight%pw%pw_grid%dvol
      dr = group(1)%weight%pw%pw_grid%dr
      np = group(1)%weight%pw%pw_grid%npts
      shift = -REAL(MODULO(np, 2), dp)*dr/2.0_dp
      DO i = 1, 3
         cell_v(i) = cell%hmat(i, i)
      END DO
      ! If requested, allocate storage for gradients
      IF (in_memory) THEN
         bo_conf = bo
         ! With confinement active, we dont need to store gradients outside
         ! the confinement bounds since they vanish for all particles
         IF (becke_control%cavity_confine) THEN
            bo_conf(1, 3) = becke_control%confine_bounds(1)
            bo_conf(2, 3) = becke_control%confine_bounds(2)
         END IF
         ALLOCATE (atom_in_group(SIZE(group), natom))
         atom_in_group = .FALSE.
         DO igroup = 1, SIZE(group)
            ALLOCATE (group(igroup)%gradients(3*natom, bo_conf(1, 1):bo_conf(2, 1), &
                                              bo_conf(1, 2):bo_conf(2, 2), &
                                              bo_conf(1, 3):bo_conf(2, 3)))
            group(igroup)%gradients = 0.0_dp
            ALLOCATE (group(igroup)%d_sum_const_dR(3, natom))
            group(igroup)%d_sum_const_dR = 0.0_dp
            DO ip = 1, SIZE(group(igroup)%atoms)
               atom_in_group(igroup, group(igroup)%atoms(ip)) = .TRUE.
            END DO
         END DO
      END IF
      ! Allocate remaining work
      ALLOCATE (sum_cell_f_group(SIZE(group)))
      ALLOCATE (cell_functions(natom))
      IF (in_memory) THEN
         ALLOCATE (ds_dR_j(3))
         ALLOCATE (ds_dR_i(3))
         ALLOCATE (d_sum_Pm_dR(3, natom))
         ALLOCATE (dP_i_dRj(3, natom, natom))
         ALLOCATE (dP_i_dRi(3, natom))
         th = 1.0e-8_dp
      END IF
      ! Build constraint
      DO k = bo(1, 1), bo(2, 1)
         DO j = bo(1, 2), bo(2, 2)
            DO i = bo(1, 3), bo(2, 3)
               ! If the grid point is too far from all constraint atoms and cavity confinement is active,
               ! we can skip this grid point as it does not contribute to the weight or gradients
               IF (becke_control%cavity_confine) THEN
                  IF (becke_control%cavity%pw%cr3d(k, j, i) < eps_cavity) CYCLE
               END IF
               ind = (/k, j, i/)
               grid_p(1) = k*dr(1)+shift(1)
               grid_p(2) = j*dr(2)+shift(2)
               grid_p(3) = i*dr(3)+shift(3)
               nskipped = 0
               cell_functions = 1.0_dp
               skip_me = .FALSE.
               IF (becke_control%vector_buffer%store_vectors) becke_control%vector_buffer%distances = 0.0_dp
               IF (in_memory) THEN
                  d_sum_Pm_dR = 0.0_dp
                  DO igroup = 1, SIZE(group)
                     group(igroup)%d_sum_const_dR = 0.0_dp
                  END DO
                  dP_i_dRi = 0.0_dp
               END IF
               ! Iterate over all atoms in the system
               DO iatom = 1, natom
                  IF (skip_me(iatom)) THEN
                     cell_functions(iatom) = 0.0_dp
                     IF (becke_control%should_skip) THEN
                        IF (is_constraint(iatom)) nskipped = nskipped+1
                        IF (nskipped == becke_control%natoms) THEN
                           IF (in_memory) THEN
                              IF (becke_control%cavity_confine) THEN
                                 becke_control%cavity%pw%cr3d(k, j, i) = 0.0_dp
                              END IF
                           END IF
                           EXIT
                        END IF
                     END IF
                     CYCLE
                  END IF
                  IF (becke_control%vector_buffer%store_vectors) THEN
                     IF (becke_control%vector_buffer%distances(iatom) .EQ. 0.0_dp) THEN
                        r = becke_control%vector_buffer%position_vecs(:, iatom)
                        dist_vec = (r-grid_p)-ANINT((r-grid_p)/cell_v)*cell_v
                        dist1 = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
                        becke_control%vector_buffer%distance_vecs(:, iatom) = dist_vec
                        becke_control%vector_buffer%distances(iatom) = dist1
                     ELSE
                        dist_vec = becke_control%vector_buffer%distance_vecs(:, iatom)
                        dist1 = becke_control%vector_buffer%distances(iatom)
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
                           IF (becke_control%vector_buffer%store_vectors) THEN
                              IF (becke_control%vector_buffer%distances(jatom) .EQ. 0.0_dp) THEN
                                 r1 = becke_control%vector_buffer%position_vecs(:, jatom)
                                 dist_vec = (r1-grid_p)-ANINT((r1-grid_p)/cell_v)*cell_v
                                 dist2 = SQRT(DOT_PRODUCT(dist_vec, dist_vec))
                                 becke_control%vector_buffer%distance_vecs(:, jatom) = dist_vec
                                 becke_control%vector_buffer%distances(jatom) = dist2
                              ELSE
                                 dist_vec = becke_control%vector_buffer%distance_vecs(:, jatom)
                                 dist2 = becke_control%vector_buffer%distances(jatom)
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
                              IF (becke_control%vector_buffer%store_vectors) THEN
                                 dr1_r2 = becke_control%vector_buffer%pair_dist_vecs(:, iatom, jatom)
                              ELSE
                                 dr1_r2 = (r-r1)-ANINT((r-r1)/cell_v)*cell_v
                              END IF
                              IF (dist2 .LE. th) dist2 = th
                              tmp_const = (becke_control%vector_buffer%R12(iatom, jatom)**3)
                              dr_ij_dR(:) = dr1_r2(:)/tmp_const
                              !derivative w.r.t. Rj
                              dr_j_dR = dist_vec(:)/dist2
                              dmy_dR_j(:) = -(dr_j_dR(:)/becke_control%vector_buffer%R12(iatom, jatom)-(dist1-dist2)*dr_ij_dR(:))
                              !derivative w.r.t. Ri
                              dmy_dR_i(:) = dr_i_dR(:)/becke_control%vector_buffer%R12(iatom, jatom)-(dist1-dist2)*dr_ij_dR(:)
                           END IF
                           ! myij
                           my1 = (dist1-dist2)/becke_control%vector_buffer%R12(iatom, jatom)
                           IF (becke_control%adjust) THEN
                              my1_homo = my1 ! Homonuclear quantity needed for gradient
                              my1 = my1+becke_control%aij(iatom, jatom)*(1.0_dp-my1**2)
                           END IF
                           ! f(myij)
                           myexp = 1.5_dp*my1-0.5_dp*my1**3
                           IF (in_memory) THEN
                              dmyexp = 1.5_dp-1.5_dp*my1**2
                              tmp_const = (1.5_dp**2)*dmyexp*(1-myexp**2)* &
                                          (1.0_dp-((1.5_dp*myexp-0.5_dp*(myexp**3))**2))
                              ! d s(myij)/d R_i
                              ds_dR_i(:) = -0.5_dp*tmp_const*dmy_dR_i(:)
                              ! d s(myij)/d R_j
                              ds_dR_j(:) = -0.5_dp*tmp_const*dmy_dR_j(:)
                              IF (becke_control%adjust) THEN
                                 tmp_const = 1.0_dp-2.0_dp*my1_homo* &
                                             becke_control%aij(iatom, jatom)
                                 ds_dR_i(:) = ds_dR_i(:)*tmp_const
                                 ! tmp_const is same for both since aij=-aji and myij=-myji
                                 ds_dR_j(:) = ds_dR_j(:)*tmp_const
                              END IF
                           END IF
                           ! s(myij) = f[f(f{myij})]
                           myexp = 1.5_dp*myexp-0.5_dp*myexp**3
                           myexp = 1.5_dp*myexp-0.5_dp*myexp**3
                           tmp_const = 0.5_dp*(1.0_dp-myexp)
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
                        ! Final value of dP_i_dRi
                        dP_i_dRi(:, iatom) = cell_functions(iatom)*dP_i_dRi(:, iatom)
                        ! Update relevant sums with value
                        d_sum_Pm_dR(:, iatom) = d_sum_Pm_dR(:, iatom)+dP_i_dRi(:, iatom)
                        IF (is_constraint(iatom)) THEN
                           DO igroup = 1, SIZE(group)
                              IF (.NOT. atom_in_group(igroup, iatom)) CYCLE
                              DO jp = 1, SIZE(group(igroup)%atoms)
                                 IF (iatom == group(igroup)%atoms(jp)) THEN
                                    ip = jp
                                    EXIT
                                 END IF
                              END DO
                              group(igroup)%d_sum_const_dR(1:3, iatom) = group(igroup)%d_sum_const_dR(1:3, iatom)+ &
                                                                         group(igroup)%coeff(ip)*dP_i_dRi(:, iatom)
                           END DO
                        END IF
                        DO jatom = 1, natom
                           IF (jatom .NE. iatom) THEN
                              ! Final value of dP_i_dRj
                              dP_i_dRj(:, iatom, jatom) = cell_functions(iatom)*dP_i_dRj(:, iatom, jatom)
                              ! Update where needed
                              d_sum_Pm_dR(:, jatom) = d_sum_Pm_dR(:, jatom)+dP_i_dRj(:, iatom, jatom)
                              IF (is_constraint(iatom)) THEN
                                 DO igroup = 1, SIZE(group)
                                    IF (.NOT. atom_in_group(igroup, iatom)) CYCLE
                                    ip = -1
                                    DO jp = 1, SIZE(group(igroup)%atoms)
                                       IF (iatom == group(igroup)%atoms(jp)) THEN
                                          ip = jp
                                          EXIT
                                       END IF
                                    END DO
                                    group(igroup)%d_sum_const_dR(1:3, jatom) = group(igroup)%d_sum_const_dR(1:3, jatom)+ &
                                                                               group(igroup)%coeff(ip)* &
                                                                               dP_i_dRj(:, iatom, jatom)
                                 END DO
                              END IF
                           END IF
                        END DO
                     END IF
                  ELSE
                     cell_functions(iatom) = 0.0_dp
                     skip_me(iatom) = .TRUE.
                     IF (becke_control%should_skip) THEN
                        IF (is_constraint(iatom)) nskipped = nskipped+1
                        IF (nskipped == becke_control%natoms) THEN
                           IF (in_memory) THEN
                              IF (becke_control%cavity_confine) THEN
                                 becke_control%cavity%pw%cr3d(k, j, i) = 0.0_dp
                              END IF
                           END IF
                           EXIT
                        END IF
                     END IF
                  END IF
               END DO !iatom
               IF (nskipped == becke_control%natoms) CYCLE
               ! Sum up cell functions
               sum_cell_f_group = 0.0_dp
               DO igroup = 1, SIZE(group)
                  DO ip = 1, SIZE(group(igroup)%atoms)
                     sum_cell_f_group(igroup) = sum_cell_f_group(igroup)+group(igroup)%coeff(ip)* &
                                                cell_functions(group(igroup)%atoms(ip))
                  END DO
               END DO
               sum_cell_f_all = 0.0_dp
               DO ip = 1, natom
                  sum_cell_f_all = sum_cell_f_all+cell_functions(ip)
               END DO
               ! Gradients at (k,j,i)
               IF (in_memory .AND. ABS(sum_cell_f_all) .GT. 0.0_dp) THEN
                  DO igroup = 1, SIZE(group)
                     DO iatom = 1, natom
                        group(igroup)%gradients(3*(iatom-1)+1:3*(iatom-1)+3, k, j, i) = &
                           group(igroup)%d_sum_const_dR(1:3, iatom)/sum_cell_f_all-sum_cell_f_group(igroup)* &
                           d_sum_Pm_dR(1:3, iatom)/(sum_cell_f_all**2)
                     END DO
                  END DO
               END IF
               ! Weight function(s) at (k,j,i)
               IF (.NOT. my_just_gradients .AND. ABS(sum_cell_f_all) .GT. 0.000001) THEN
                  DO igroup = 1, SIZE(group)
                     group(igroup)%weight%pw%cr3d(k, j, i) = sum_cell_f_group(igroup)/sum_cell_f_all
                  END DO
                  IF (becke_control%atomic_charges) THEN
                     DO iatom = 1, becke_control%natoms
                        charge(iatom)%pw%cr3d(k, j, i) = cell_functions(catom(iatom))/sum_cell_f_all
                     END DO
                  END IF
               END IF
            END DO
         END DO
      END DO
      ! Release storage
      IF (in_memory) THEN
         DEALLOCATE (ds_dR_j)
         DEALLOCATE (ds_dR_i)
         DEALLOCATE (d_sum_Pm_dR)
         DEALLOCATE (dP_i_dRj)
         DEALLOCATE (dP_i_dRi)
         DO igroup = 1, SIZE(group)
            DEALLOCATE (group(igroup)%d_sum_const_dR)
         END DO
         DEALLOCATE (atom_in_group)
         IF (becke_control%vector_buffer%store_vectors) THEN
            DEALLOCATE (becke_control%vector_buffer%pair_dist_vecs)
         END IF
      END IF
      NULLIFY (cutoffs)
      IF (ALLOCATED(is_constraint)) &
         DEALLOCATE (is_constraint)
      DEALLOCATE (catom)
      DEALLOCATE (cell_functions)
      DEALLOCATE (skip_me)
      DEALLOCATE (sum_cell_f_group)
      DEALLOCATE (becke_control%vector_buffer%R12)
      IF (becke_control%vector_buffer%store_vectors) THEN
         DEALLOCATE (becke_control%vector_buffer%distances)
         DEALLOCATE (becke_control%vector_buffer%distance_vecs)
         DEALLOCATE (becke_control%vector_buffer%position_vecs)
      END IF
      CALL timestop(handle)

   END SUBROUTINE becke_constraint_low

! **************************************************************************************************
!> \brief Calculates atomic forces due to a Becke constraint
!> \param qs_env the qs_env where to build the gradients
!> \par   History
!>        Created 01.2007 [fschiff]
!>        Extended functionality 12/15-12/16 [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE becke_constraint_force(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'becke_constraint_force', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iatom, igroup, ikind, ind(3), &
                                                            ispin, j, k, natom, np(3), nvar
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(kind=dp)                                      :: dvol, eps_cavity, sign, th
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: strength
      REAL(kind=dp), DIMENSION(3)                        :: dr, shift
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cutoffs
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(becke_constraint_type), POINTER               :: becke_control
      TYPE(becke_group_type), DIMENSION(:), POINTER      :: group
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)
      NULLIFY (atomic_kind_set, cell, para_env, dft_control, particle_set, &
               rho, rho_r, force, cutoffs, becke_control, group)

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
      becke_control => dft_control%qs_control%becke_control
      group => becke_control%group
      nvar = SIZE(becke_control%target)
      ALLOCATE (strength(nvar))
      strength(:) = becke_control%strength(:)
      cutoffs => becke_control%cutoffs
      eps_cavity = becke_control%eps_cavity
      ALLOCATE (atom_of_kind(natom))
      ALLOCATE (kind_of(natom))

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, &
                               kind_of=kind_of)
      DO igroup = 1, SIZE(group)
         ALLOCATE (group(igroup)%integrated(3, natom))
         group(igroup)%integrated = 0.0_dp
      END DO
      bo = group(1)%weight%pw%pw_grid%bounds_local
      np = group(1)%weight%pw%pw_grid%npts
      dvol = group(1)%weight%pw%pw_grid%dvol
      dr = group(1)%weight%pw%pw_grid%dr
      np = group(1)%weight%pw%pw_grid%npts
      shift = -REAL(MODULO(np, 2), dp)*dr/2.0_dp

      IF (.NOT. becke_control%in_memory) &
         ! Gradients were not precomputed => calculate them now
         CALL becke_constraint_low(qs_env, just_gradients=.TRUE.)
      ! Integrate gradients
      IF (.NOT. ASSOCIATED(becke_control%cavity_mat)) THEN
         ! No external control
         DO k = bo(1, 1), bo(2, 1)
            DO j = bo(1, 2), bo(2, 2)
               DO i = bo(1, 3), bo(2, 3)
                  ! First check if this grid point should be skipped
                  IF (becke_control%cavity_confine) THEN
                     IF (becke_control%cavity%pw%cr3d(k, j, i) < eps_cavity) CYCLE
                  END IF
                  ind = (/k, j, i/)
                  DO igroup = 1, SIZE(group)
                     DO iatom = 1, natom
                        DO ispin = 1, dft_control%nspins
                           SELECT CASE (group (igroup)%constraint_type)
                           CASE (cdft_charge_constraint)
                              sign = 1.0_dp
                           CASE (cdft_magnetization_constraint)
                              IF (ispin == 1) THEN
                                 sign = 1.0_dp
                              ELSE
                                 sign = -1.0_dp
                              END IF
                           CASE (cdft_alpha_constraint)
                              sign = 1.0_dp
                              IF (ispin == 2) CYCLE
                           CASE (cdft_beta_constraint)
                              sign = 1.0_dp
                              IF (ispin == 1) CYCLE
                           CASE DEFAULT
                              CPABORT("Unknown constraint type.")
                           END SELECT
                           group(igroup)%integrated(:, iatom) = group(igroup)%integrated(:, iatom)+sign* &
                                                                group(igroup)%gradients(3*(iatom-1)+1:3*(iatom-1)+3, k, j, i)* &
                                                                rho_r(ispin)%pw%cr3d(k, j, i)*dvol
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      ELSE
         DO k = LBOUND(becke_control%cavity_mat, 1), UBOUND(becke_control%cavity_mat, 1)
            DO j = LBOUND(becke_control%cavity_mat, 2), UBOUND(becke_control%cavity_mat, 2)
               DO i = LBOUND(becke_control%cavity_mat, 3), UBOUND(becke_control%cavity_mat, 3)
                  ! First check if this grid point should be skipped
                  IF (becke_control%cavity_mat(k, j, i) < eps_cavity) CYCLE
                  DO igroup = 1, SIZE(group)
                     DO iatom = 1, natom
                        DO ispin = 1, dft_control%nspins
                           SELECT CASE (group (igroup)%constraint_type)
                           CASE (cdft_charge_constraint)
                              sign = 1.0_dp
                           CASE (cdft_magnetization_constraint)
                              IF (ispin == 1) THEN
                                 sign = 1.0_dp
                              ELSE
                                 sign = -1.0_dp
                              END IF
                           CASE (cdft_alpha_constraint)
                              sign = 1.0_dp
                              IF (ispin == 2) CYCLE
                           CASE (cdft_beta_constraint)
                              sign = 1.0_dp
                              IF (ispin == 1) CYCLE
                           CASE DEFAULT
                              CPABORT("Unknown constraint type.")
                           END SELECT
                           group(igroup)%integrated(:, iatom) = group(igroup)%integrated(:, iatom)+sign* &
                                                                group(igroup)%gradients(3*(iatom-1)+1:3*(iatom-1)+3, k, j, i)* &
                                                                rho_r(ispin)%pw%cr3d(k, j, i)*dvol
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END IF
      IF (.NOT. becke_control%transfer_pot) THEN
         DO igroup = 1, SIZE(group)
            DEALLOCATE (group(igroup)%gradients)
         END DO
      END IF
      DO igroup = 1, SIZE(group)
         CALL mp_sum(group(igroup)%integrated, para_env%group)
      END DO
      ! Update force only on master process. Otherwise force due to constraint becomes multiplied
      ! by the number of processes when the final force%rho_elec is constructed in qs_force
      ! by mp_summing [the final integrated(:,:) is distributed on all processors]
      IF (para_env%mepos == para_env%source) THEN
         DO igroup = 1, SIZE(group)
            DO iatom = 1, natom
               ikind = kind_of(iatom)
               i = atom_of_kind(iatom)
               force(ikind)%rho_elec(:, i) = force(ikind)%rho_elec(:, i)+group(igroup)%integrated(:, iatom)*strength(igroup)
            END DO
         END DO
      END IF
      DEALLOCATE (strength)
      DEALLOCATE (atom_of_kind)
      DEALLOCATE (kind_of)
      DO igroup = 1, SIZE(group)
         DEALLOCATE (group(igroup)%integrated)
      END DO
      NULLIFY (group)

      CALL timestop(handle)

   END SUBROUTINE becke_constraint_force

! **************************************************************************************************
!> \brief Calculates the value of a Becke constraint by integrating the product of the Becke
!>        weight function and the realspace electron density
!> \param qs_env the qs_env where to build the constraint
!> \par   History
!>        Created 3.2017 [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE becke_constraint_integrate(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'becke_constraint_integrate', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=2)                                   :: element_symbol
      INTEGER                                            :: handle, i, iatom, igroup, ikind, ivar, &
                                                            iw, jatom, natom, nelectron_total, &
                                                            nfrag_spins, nvar, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: stride
      LOGICAL                                            :: mpi_io, needs_spin_density, paw_atom
      REAL(kind=dp)                                      :: dvol, eps_cavity, multiplier(2), &
                                                            nelectron_frag, sign, tc(2), zeff
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: dE, strength, target_val
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: electronic_charge, gapw_offset
      TYPE(becke_constraint_type), POINTER               :: becke_control
      TYPE(becke_group_type), DIMENSION(:), POINTER      :: group
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mpole_rho_atom), DIMENSION(:), POINTER        :: mp_rho
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: charge, rho_frag, rho_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(rho0_mpole_type), POINTER                     :: rho0_mpole
      TYPE(section_vals_type), POINTER                   :: becke_constraint_section

      NULLIFY (para_env, dft_control, particle_set, rho_r, energy, rho, logger, becke_constraint_section, qs_kind_set, &
               particles, subsys, pw_env, auxbas_pw_pool, mp_rho, rho0_mpole, group, charge, rho_frag)
      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      CALL get_qs_env(qs_env, &
                      particle_set=particle_set, &
                      rho=rho, &
                      natom=natom, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      qs_kind_set=qs_kind_set)
      CALL qs_rho_get(rho, rho_r=rho_r)
      CPASSERT(ASSOCIATED(qs_kind_set))
      becke_constraint_section => section_vals_get_subs_vals(qs_env%input, "DFT%QS%BECKE_CONSTRAINT")
      iw = cp_print_key_unit_nr(logger, becke_constraint_section, "PROGRAM_RUN_INFO", extension=".cdftLog")
      becke_control => dft_control%qs_control%becke_control
      group => becke_control%group
      !
      nvar = SIZE(becke_control%target)
      ALLOCATE (strength(nvar))
      ALLOCATE (target_val(nvar))
      ALLOCATE (dE(nvar))
      strength(:) = becke_control%strength(:)
      target_val(:) = becke_control%target(:)
      dE = 0.0_dp
      dvol = group(1)%weight%pw%pw_grid%dvol
      eps_cavity = becke_control%eps_cavity
      IF (becke_control%atomic_charges) THEN
         charge => becke_control%charge
         ALLOCATE (electronic_charge(becke_control%natoms, dft_control%nspins))
         electronic_charge = 0.0_dp
      END IF
      ! Calculate value of constraint i.e. int ( rho(r) w(r) dr)
      DO i = 1, dft_control%nspins
         DO igroup = 1, SIZE(group)
            SELECT CASE (group (igroup)%constraint_type)
            CASE (cdft_charge_constraint)
               sign = 1.0_dp
            CASE (cdft_magnetization_constraint)
               IF (i == 1) THEN
                  sign = 1.0_dp
               ELSE
                  sign = -1.0_dp
               END IF
            CASE (cdft_alpha_constraint)
               sign = 1.0_dp
               IF (i == 2) CYCLE
            CASE (cdft_beta_constraint)
               sign = 1.0_dp
               IF (i == 1) CYCLE
            CASE DEFAULT
               CPABORT("Unknown constraint type.")
            END SELECT
            IF (becke_control%external_control .AND. becke_control%cavity_confine) THEN
               ! With external control, we can use cavity_mat as a mask to kahan sum
               IF (igroup /= 1) &
                  CALL cp_abort(__LOCATION__, &
                                "Multiple constraints not yet supported by parallel mixed calculations.")
               dE(igroup) = dE(igroup)+sign*accurate_dot_product(group(igroup)%weight%pw%cr3d, rho_r(i)%pw%cr3d, &
                                                                 becke_control%cavity_mat, eps_cavity)*dvol
            ELSE
               dE(igroup) = dE(igroup)+sign*accurate_sum(group(igroup)%weight%pw%cr3d*rho_r(i)%pw%cr3d)*dvol
            END IF
         END DO
         IF (becke_control%atomic_charges) THEN
            DO iatom = 1, becke_control%natoms
               electronic_charge(iatom, i) = accurate_sum(charge(iatom)%pw%cr3d*rho_r(i)%pw%cr3d)*dvol
            END DO
         END IF
      END DO
      CALL get_qs_env(qs_env, energy=energy)
      CALL mp_sum(dE, para_env%group)
      IF (becke_control%atomic_charges) THEN
         CALL mp_sum(electronic_charge, para_env%group)
      END IF
      ! Use fragment densities as reference value (= Becke deformation density)
      IF (becke_control%fragment_density) THEN
         ! Fragment densities are meaningful only for some calculation types
         IF (.NOT. qs_env%single_point_run) &
            CALL cp_abort(__LOCATION__, &
                          "Becke deformation density only compatible with single "// &
                          "point calculations (run_type ENERGY or ENERGY_FORCE).")
         IF (dft_control%qs_control%gapw) &
            CALL cp_abort(__LOCATION__, &
                          "Becke deformation density not compatible with GAPW.")
         needs_spin_density = .FALSE.
         multiplier = 1.0_dp
         nfrag_spins = 1
         DO igroup = 1, SIZE(group)
            SELECT CASE (group (igroup)%constraint_type)
            CASE (cdft_charge_constraint)
               ! Do nothing
            CASE (cdft_magnetization_constraint)
               needs_spin_density = .TRUE.
            CASE (cdft_alpha_constraint, cdft_beta_constraint)
               CALL cp_abort(__LOCATION__, &
                             "Becke deformation density not yet compatible with "// &
                             "spin specific constraints.")
            CASE DEFAULT
               CPABORT("Unknown constraint type.")
            END SELECT
         END DO
         IF (needs_spin_density) THEN
            nfrag_spins = 2
            DO i = 1, 2
               IF (becke_control%flip_fragment(i)) multiplier(i) = -1.0_dp
            END DO
         END IF
         IF (.NOT. becke_control%fragments_integrated) THEN
            ALLOCATE (becke_control%fragments(nfrag_spins, 2))
            ALLOCATE (rho_frag(nfrag_spins))
            CALL get_qs_env(qs_env, pw_env=pw_env)
            CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
            ! Read fragment reference densities
            ! Total density (rho_alpha + rho_beta)
            CALL pw_pool_create_pw(auxbas_pw_pool, becke_control%fragments(1, 1)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL cube_to_pw(becke_control%fragments(1, 1)%pw, &
                            becke_control%fragment_a_fname, 1.0_dp)
            CALL pw_pool_create_pw(auxbas_pw_pool, becke_control%fragments(1, 2)%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL cube_to_pw(becke_control%fragments(1, 2)%pw, &
                            becke_control%fragment_b_fname, 1.0_dp)
            ! Spin difference density (rho_alpha - rho_beta) if needed
            IF (needs_spin_density) THEN
               CALL pw_pool_create_pw(auxbas_pw_pool, becke_control%fragments(2, 1)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL cube_to_pw(becke_control%fragments(2, 1)%pw, &
                               becke_control%fragment_a_spin_fname, multiplier(1))
               CALL pw_pool_create_pw(auxbas_pw_pool, becke_control%fragments(2, 2)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL cube_to_pw(becke_control%fragments(2, 2)%pw, &
                               becke_control%fragment_b_spin_fname, multiplier(2))
            END IF
            ! Sum up fragments
            DO i = 1, nfrag_spins
               CALL pw_pool_create_pw(auxbas_pw_pool, rho_frag(i)%pw, use_data=REALDATA3D, &
                                      in_space=REALSPACE)
               CALL pw_copy(becke_control%fragments(i, 1)%pw, rho_frag(i)%pw)
               CALL pw_axpy(becke_control%fragments(i, 2)%pw, rho_frag(i)%pw, 1.0_dp)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, becke_control%fragments(i, 1)%pw)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, becke_control%fragments(i, 2)%pw)
            END DO
            DEALLOCATE (becke_control%fragments)
            ! Check that the number of electrons is consistent
            CALL get_qs_env(qs_env, subsys=subsys)
            CALL qs_subsys_get(subsys, nelectron_total=nelectron_total)
            nelectron_frag = pw_integrate_function(rho_frag(1)%pw)
            IF (NINT(nelectron_frag) /= nelectron_total) &
               CALL cp_abort(__LOCATION__, &
                             "The number of electrons in the reference and interacting "// &
                             " configurations does not match. Check your fragment cube files.")
            ! Update constraint target value i.e. perform integration w_i*rho_frag_{tot/spin}*dr
            becke_control%target = 0.0_dp
            DO igroup = 1, SIZE(group)
               IF (group(igroup)%constraint_type == cdft_charge_constraint) THEN
                  i = 1
               ELSE
                  i = 2
               END IF
               IF (becke_control%external_control .AND. becke_control%cavity_confine) THEN
                  becke_control%target(igroup) = becke_control%target(igroup)+ &
                                                 accurate_dot_product(group(igroup)%weight%pw%cr3d, rho_frag(i)%pw%cr3d, &
                                                                      becke_control%cavity_mat, eps_cavity)*dvol
               ELSE
                  becke_control%target(igroup) = becke_control%target(igroup)+ &
                                                 accurate_sum(group(igroup)%weight%pw%cr3d*rho_frag(i)%pw%cr3d)*dvol
               END IF
            END DO
            CALL mp_sum(becke_control%target, para_env%group)
            IF (debug_this_module) THEN
               ALLOCATE (stride(3))
               stride = (/1, 1, 1/)
               mpi_io = .TRUE.
               unit_nr = cp_print_key_unit_nr(logger, becke_constraint_section, "PROGRAM_RUN_INFO", &
                                              middle_name="BECKE_WEIGHT", &
                                              extension=".cube", file_position="REWIND", &
                                              log_filename=.FALSE., mpi_io=mpi_io)
               CALL get_qs_env(qs_env, subsys=subsys)
               CALL qs_subsys_get(subsys, particles=particles)
               CALL cp_pw_to_cube(becke_control%group(1)%weight%pw, &
                                  unit_nr, "BECKE_WEIGHT", particles=particles, &
                                  stride=stride, mpi_io=mpi_io)
               CALL cp_print_key_finished_output(unit_nr, logger, becke_constraint_section, "PROGRAM_RUN_INFO", mpi_io=mpi_io)
               unit_nr = cp_print_key_unit_nr(logger, becke_constraint_section, "PROGRAM_RUN_INFO", &
                                              middle_name="FRAG_DENS", &
                                              extension=".cube", file_position="REWIND", &
                                              log_filename=.FALSE., mpi_io=mpi_io)
               CALL cp_pw_to_cube(rho_frag(1)%pw, &
                                  unit_nr, "FRAG_DENS", particles=particles, &
                                  stride=stride, mpi_io=mpi_io)
               CALL cp_print_key_finished_output(unit_nr, logger, becke_constraint_section, "PROGRAM_RUN_INFO", mpi_io=mpi_io)
               DEALLOCATE (stride)
            END IF
            ! Calculate reference atomic charges int( w_i * rho_frag * dr )
            IF (becke_control%atomic_charges) THEN
               ALLOCATE (becke_control%charges_fragment(becke_control%natoms, nfrag_spins))
               DO i = 1, nfrag_spins
                  DO iatom = 1, becke_control%natoms
                     becke_control%charges_fragment(iatom, i) = accurate_sum(charge(iatom)%pw%cr3d*rho_frag(i)%pw%cr3d)*dvol
                  END DO
               END DO
               CALL mp_sum(becke_control%charges_fragment, para_env%group)
            END IF
            DO i = 1, nfrag_spins
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_frag(i)%pw)
            END DO
            DEALLOCATE (rho_frag)
            target_val(:) = becke_control%target(:)
            IF (dft_control%qs_control%cdft) dft_control%qs_control%cdft_control%target = target_val
            becke_control%fragments_integrated = .TRUE.
         END IF
      END IF
      IF (dft_control%qs_control%gapw) THEN
         ! GAPW: add core charges (rho_hard - rho_soft)
         IF (becke_control%fragment_density) &
            CALL cp_abort(__LOCATION__, &
                          "Becke deformation density not yet compatible with GAPW.")
         ALLOCATE (gapw_offset(nvar, dft_control%nspins))
         gapw_offset = 0.0_dp
         CALL get_qs_env(qs_env, rho0_mpole=rho0_mpole)
         CALL get_rho0_mpole(rho0_mpole, mp_rho=mp_rho)
         DO i = 1, dft_control%nspins
            DO igroup = 1, SIZE(group)
               DO iatom = 1, SIZE(group(igroup)%atoms)
                  SELECT CASE (group (igroup)%constraint_type)
                  CASE (cdft_charge_constraint)
                     sign = 1.0_dp
                  CASE (cdft_magnetization_constraint)
                     IF (i == 1) THEN
                        sign = 1.0_dp
                     ELSE
                        sign = -1.0_dp
                     END IF
                  CASE (cdft_alpha_constraint)
                     sign = 1.0_dp
                     IF (i == 2) CYCLE
                  CASE (cdft_beta_constraint)
                     sign = 1.0_dp
                     IF (i == 1) CYCLE
                  CASE DEFAULT
                     CPABORT("Unknown constraint type.")
                  END SELECT
                  jatom = group(igroup)%atoms(iatom)
                  CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=ikind)
                  CALL get_qs_kind(qs_kind_set(ikind), paw_atom=paw_atom)
                  IF (paw_atom) THEN
                     gapw_offset(igroup, i) = gapw_offset(igroup, i)+sign*group(igroup)%coeff(iatom)*mp_rho(jatom)%q0(i)
                  END IF
               END DO
            END DO
         END DO
         IF (becke_control%atomic_charges) THEN
            DO iatom = 1, becke_control%natoms
               jatom = becke_control%atoms(iatom)
               CALL get_atomic_kind(particle_set(jatom)%atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), paw_atom=paw_atom)
               IF (paw_atom) THEN
                  DO i = 1, dft_control%nspins
                     electronic_charge(iatom, i) = electronic_charge(iatom, i)+mp_rho(jatom)%q0(i)
                  END DO
               END IF
            END DO
         END IF
         DO i = 1, dft_control%nspins
            DO ivar = 1, nvar
               dE(ivar) = dE(ivar)+gapw_offset(ivar, i)
            END DO
         END DO
         DEALLOCATE (gapw_offset)
      END IF
      ! Update constraint value and energy
      becke_control%becke_order_p(:) = dE(:)
      energy%becke = 0.0_dp
      DO ivar = 1, nvar
         energy%becke = energy%becke+(dE(ivar)-target_val(ivar))*strength(ivar)
      END DO
      IF (dft_control%qs_control%cdft) dft_control%qs_control%cdft_control%value(:) = dE(:)
      ! Print constraint info and atomic Becke charges
      IF (iw > 0) THEN
         WRITE (iw, '(/,T3,A)') &
            '--------------------- Becke constraint information ---------------------'
         DO igroup = 1, SIZE(group)
            IF (igroup > 1) WRITE (iw, '(T3,A)') ' '
            WRITE (iw, '(T3,A,T54,(3X,I18))') &
               'Atomic group                :', igroup
            SELECT CASE (group (igroup)%constraint_type)
            CASE (cdft_charge_constraint)
               IF (group(igroup)%is_fragment_constraint) THEN
                  WRITE (iw, '(T3,A,T42,A)') &
                     'Type of constraint          :', ADJUSTR('Charge density constraint (frag.)')
               ELSE
                  WRITE (iw, '(T3,A,T50,A)') &
                     'Type of constraint          :', ADJUSTR('Charge density constraint')
               END IF
            CASE (cdft_magnetization_constraint)
               IF (group(igroup)%is_fragment_constraint) THEN
                  WRITE (iw, '(T3,A,T35,A)') &
                     'Type of constraint          :', ADJUSTR('Magnetization density constraint (frag.)')
               ELSE
                  WRITE (iw, '(T3,A,T43,A)') &
                     'Type of constraint          :', ADJUSTR('Magnetization density constraint')
               END IF
            CASE (cdft_alpha_constraint)
               IF (group(igroup)%is_fragment_constraint) THEN
                  WRITE (iw, '(T3,A,T38,A)') &
                     'Type of constraint          :', ADJUSTR('Alpha spin density constraint (frag.)')
               ELSE
                  WRITE (iw, '(T3,A,T46,A)') &
                     'Type of constraint          :', ADJUSTR('Alpha spin density constraint')
               END IF
            CASE (cdft_beta_constraint)
               IF (group(igroup)%is_fragment_constraint) THEN
                  WRITE (iw, '(T3,A,T39,A)') &
                     'Type of constraint          :', ADJUSTR('Beta spin density constraint (frag.)')
               ELSE
                  WRITE (iw, '(T3,A,T47,A)') &
                     'Type of constraint          :', ADJUSTR('Beta spin density constraint')
               END IF
            CASE DEFAULT
               CPABORT("Unknown constraint type.")
            END SELECT
            WRITE (iw, '(T3,A,T54,(3X,F18.12))') &
               'Target value of constraint  :', target_val(igroup)
            WRITE (iw, '(T3,A,T54,(3X,F18.12))') &
               'Current value of constraint :', dE(igroup)
            WRITE (iw, '(T3,A,T59,(3X,ES13.3))') &
               'Deviation from target       :', dE(igroup)-target_val(igroup)
            WRITE (iw, '(T3,A,T54,(3X,F18.12))') &
               'Strength of constraint      :', strength(igroup)
         END DO
         WRITE (iw, '(T3,A)') &
            '------------------------------------------------------------------------'
         IF (becke_control%atomic_charges) THEN
            IF (.NOT. becke_control%fragment_density) THEN
               IF (dft_control%nspins == 1) THEN
                  WRITE (iw, '(/,T3,A)') &
                     '-------------------------------- Becke atomic charges --------------------------------'
                  WRITE (iw, '(T3,A,A)') &
                     '#Atom  Element   Is_constraint', '   Core charge    Population (total)'// &
                     '          Net charge'
                  tc = 0.0_dp
                  DO iatom = 1, becke_control%natoms
                     jatom = becke_control%atoms(iatom)
                     CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, &
                                          element_symbol=element_symbol, &
                                          kind_number=ikind)
                     CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
                     WRITE (iw, "(i7,T15,A2,T23,L10,T39,F8.3,T61,F8.3,T81,F8.3)") &
                        jatom, ADJUSTR(element_symbol), becke_control%is_constraint(iatom), zeff, electronic_charge(iatom, 1), &
                        (zeff-electronic_charge(iatom, 1))
                     tc(1) = tc(1)+(zeff-electronic_charge(iatom, 1))
                  END DO
                  WRITE (iw, '(/,T3,A,T81,F8.3,/)') "Total Charge: ", tc(1)
               ELSE
                  WRITE (iw, '(/,T3,A)') &
                     '------------------------------------------ Becke atomic charges -------------------------------------------'
                  WRITE (iw, '(T3,A,A)') &
                     '#Atom  Element   Is_constraint', '   Core charge    Population (alpha, beta)'// &
                     '    Net charge      Spin population'
                  tc = 0.0_dp
                  DO iatom = 1, becke_control%natoms
                     jatom = becke_control%atoms(iatom)
                     CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, &
                                          element_symbol=element_symbol, &
                                          kind_number=ikind)
                     CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
                     WRITE (iw, "(i7,T15,A2,T23,L10,T39,F8.3,T53,F8.3,T67,F8.3,T81,F8.3,T102,F8.3)") &
                        jatom, ADJUSTR(element_symbol), &
                        becke_control%is_constraint(iatom), &
                        zeff, electronic_charge(iatom, 1), electronic_charge(iatom, 2), &
                        (zeff-electronic_charge(iatom, 1)-electronic_charge(iatom, 2)), &
                        electronic_charge(iatom, 1)-electronic_charge(iatom, 2)
                     tc(1) = tc(1)+(zeff-electronic_charge(iatom, 1)-electronic_charge(iatom, 2))
                     tc(2) = tc(2)+(electronic_charge(iatom, 1)-electronic_charge(iatom, 2))
                  END DO
                  WRITE (iw, '(/,T3,A,T81,F8.3,T102,F8.3/)') "Total Charge and Spin Moment: ", tc(1), tc(2)
               END IF
            ELSE
               IF (.NOT. needs_spin_density) THEN
                  WRITE (iw, '(/,T3,A)') &
                     '-------------------------------- Becke atomic charges --------------------------------'
                  IF (dft_control%nspins == 1) THEN
                     WRITE (iw, '(T3,A,A)') &
                        '#Atom  Element   Is_constraint', '   Fragment charge    Population (total)'// &
                        '      Net charge'
                  ELSE
                     WRITE (iw, '(T3,A,A)') &
                        '#Atom  Element   Is_constraint', '   Fragment charge  Population (alpha, beta)'// &
                        '  Net charge'
                  END IF
                  tc = 0.0_dp
                  DO iatom = 1, becke_control%natoms
                     jatom = becke_control%atoms(iatom)
                     CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, &
                                          element_symbol=element_symbol, &
                                          kind_number=ikind)
                     IF (dft_control%nspins == 1) THEN
                        WRITE (iw, "(i7,T15,A2,T23,L10,T43,F8.3,T65,F8.3,T81,F8.3)") &
                           jatom, ADJUSTR(element_symbol), &
                           becke_control%is_constraint(iatom), &
                           becke_control%charges_fragment(iatom, 1), &
                           electronic_charge(iatom, 1), &
                           (electronic_charge(iatom, 1)- &
                            becke_control%charges_fragment(iatom, 1))
                        tc(1) = tc(1)+(electronic_charge(iatom, 1)- &
                                       becke_control%charges_fragment(iatom, 1))
                     ELSE
                        WRITE (iw, "(i7,T15,A2,T23,L10,T43,F8.3,T57,F8.3,T69,F8.3,T81,F8.3)") &
                           jatom, ADJUSTR(element_symbol), &
                           becke_control%is_constraint(iatom), &
                           becke_control%charges_fragment(iatom, 1), &
                           electronic_charge(iatom, 1), electronic_charge(iatom, 2), &
                           (electronic_charge(iatom, 1)+electronic_charge(iatom, 2)- &
                            becke_control%charges_fragment(iatom, 1))
                        tc(1) = tc(1)+(electronic_charge(iatom, 1)+electronic_charge(iatom, 2)- &
                                       becke_control%charges_fragment(iatom, 1))
                     END IF
                  END DO
                  WRITE (iw, '(/,T3,A,T81,F8.3,/)') "Total Charge: ", tc(1)
               ELSE
                  WRITE (iw, '(/,T3,A)') &
                     '------------------------------------------ Becke atomic charges -------------------------------------------'
                  WRITE (iw, '(T3,A,A)') &
                     '#Atom  Element  Is_constraint', ' Fragment charge/spin moment'// &
                     '  Population (alpha, beta)  Net charge/spin moment'
                  tc = 0.0_dp
                  DO iatom = 1, becke_control%natoms
                     jatom = becke_control%atoms(iatom)
                     CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, &
                                          element_symbol=element_symbol, &
                                          kind_number=ikind)
                     WRITE (iw, "(i7,T15,A2,T22,L10,T40,F8.3,T52,F8.3,T66,F8.3,T78,F8.3,T90,F8.3,T102,F8.3)") &
                        jatom, ADJUSTR(element_symbol), &
                        becke_control%is_constraint(iatom), &
                        becke_control%charges_fragment(iatom, 1), &
                        becke_control%charges_fragment(iatom, 2), &
                        electronic_charge(iatom, 1), electronic_charge(iatom, 2), &
                        (electronic_charge(iatom, 1)+electronic_charge(iatom, 2)- &
                         becke_control%charges_fragment(iatom, 1)), &
                        (electronic_charge(iatom, 1)-electronic_charge(iatom, 2)- &
                         becke_control%charges_fragment(iatom, 2))
                     tc(1) = tc(1)+(electronic_charge(iatom, 1)+electronic_charge(iatom, 2)- &
                                    becke_control%charges_fragment(iatom, 1))
                     tc(2) = tc(2)+(electronic_charge(iatom, 1)-electronic_charge(iatom, 2)- &
                                    becke_control%charges_fragment(iatom, 2))
                  END DO
                  WRITE (iw, '(/,T3,A,T90,F8.3,T102,F8.3/)') "Total Charge and Spin Moment: ", tc(1), tc(2)
               END IF
            END IF
         END IF
      END IF
      DEALLOCATE (dE, strength, target_val)
      IF (becke_control%atomic_charges) &
         DEALLOCATE (electronic_charge)
      CALL cp_print_key_finished_output(iw, logger, becke_constraint_section, "PROGRAM_RUN_INFO")
      CALL timestop(handle)

   END SUBROUTINE becke_constraint_integrate

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
                            hirshfeld_env%fnorm%pw%cr3d, divide=.TRUE.)
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
                               hirshfeld_env%fnorm%pw%cr3d, divide=.TRUE.)
            ELSE
               CPABORT("We need rho in real space")
            END IF
            ! Weight for non-constraint atoms
            CALL hfun_scale(tmp%pw%cr3d, rhonorm%pw%cr3d, weight%pw%cr3d, divide=.FALSE.)
            CALL rs_pw_transfer(weight_rs, tmp%pw, pw2rs)
            ! Weight for donor atoms
            CALL hfun_scale(tmp%pw%cr3d, rhonorm%pw%cr3d, donor_weight%pw%cr3d, divide=.FALSE.)
            CALL rs_pw_transfer(donor_weight_rs, tmp%pw, pw2rs)
            ! Weight for acceptor atoms
            CALL hfun_scale(tmp%pw%cr3d, rhonorm%pw%cr3d, acceptor_weight%pw%cr3d, divide=.FALSE.)
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

END MODULE qs_cdft_methods
