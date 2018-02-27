!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility subroutines for CDFT calculations
!> \par   History
!>                 separated from et_coupling [03.2017]
!> \author Nico Holmberg [03.2017]
! **************************************************************************************************
MODULE qs_cdft_utils
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE bibliography,                    ONLY: Becke1988b,&
                                              Holmberg2017,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type,&
                                              qs_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE hirshfeld_methods,               ONLY: create_shape_function
   USE hirshfeld_types,                 ONLY: create_hirshfeld_type,&
                                              hirshfeld_type,&
                                              set_hirshfeld_info
   USE input_constants,                 ONLY: &
        becke_cutoff_element, becke_cutoff_global, cdft_charge_constraint, &
        outer_scf_becke_constraint, outer_scf_cdft_constraint, outer_scf_hirshfeld_constraint, &
        outer_scf_none, outer_scf_optimizer_broyden, outer_scf_optimizer_newton, &
        outer_scf_optimizer_newton_ls, radius_covalent, radius_user, shape_function_gaussian
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE memory_utilities,                ONLY: reallocate
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE
   USE qs_cdft_types,                   ONLY: becke_constraint_type,&
                                              becke_group_type
   USE qs_collocate_density,            ONLY: collocate_pgf_product_rspace
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_modify_pab_block,             ONLY: FUNC_AB
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                              rs2pw,&
                                              rs_grid_release,&
                                              rs_grid_retain,&
                                              rs_grid_zero,&
                                              rs_pw_transfer
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_cdft_utils'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

! *** Public subroutines ***
   PUBLIC :: becke_constraint_init, read_becke_section, read_cdft_control_section
   PUBLIC :: hfun_scale

CONTAINS

! **************************************************************************************************
!> \brief Initializes the Becke constraint environment
!> \param qs_env the qs_env where to build the constraint
!> \par   History
!>        Created 01.2007 [fschiff]
!>        Extended functionality 12/15-12/16 [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE becke_constraint_init(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'becke_constraint_init', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=2)                                   :: element_symbol
      INTEGER :: atom_a, bounds(2), handle, i, iatom, iex, igroup, ikind, ip, ithread, iw, j, &
         jatom, katom, natom, nkind, npme, nthread, numexp, unit_nr
      INTEGER, DIMENSION(2, 3)                           :: bo
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores, stride
      LOGICAL                                            :: build, in_memory, mpi_io
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: is_constraint
      REAL(KIND=dp)                                      :: alpha, chi, coef, eps_cavity, ircov, &
                                                            jrcov, uij
      REAL(KIND=dp), DIMENSION(3)                        :: cell_v, dist_vec, r, r1, ra
      REAL(KIND=dp), DIMENSION(:), POINTER               :: radii_list
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(becke_constraint_type), POINTER               :: becke_control
      TYPE(becke_group_type), DIMENSION(:), POINTER      :: group
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hirshfeld_type), POINTER                      :: cavity_env
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(realspace_grid_type), POINTER                 :: rs_cavity
      TYPE(section_vals_type), POINTER                   :: becke_constraint_section

      NULLIFY (cores, stride, atom_list, cell, para_env, dft_control, &
               particle_set, logger, becke_constraint_section, qs_kind_set, &
               particles, subsys, pab, pw_env, rs_cavity, cavity_env, &
               auxbas_pw_pool, atomic_kind_set, group, radii_list)
      logger => cp_get_default_logger()
      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      particle_set=particle_set, &
                      natom=natom, &
                      dft_control=dft_control, &
                      para_env=para_env)
      becke_constraint_section => section_vals_get_subs_vals(qs_env%input, "DFT%QS%BECKE_CONSTRAINT")
      iw = cp_print_key_unit_nr(logger, becke_constraint_section, "PROGRAM_RUN_INFO", extension=".cdftLog")
      becke_control => dft_control%qs_control%becke_control
      group => becke_control%group
      ! Sanity checks
      DO igroup = 1, SIZE(group)
         IF (group(igroup)%constraint_type /= cdft_charge_constraint .AND. dft_control%nspins == 1) &
            CALL cp_abort(__LOCATION__, &
                          "Becke magnetization density constraint requires UKS calculation.")
      END DO
      IF (natom < becke_control%natoms) &
         CALL cp_abort(__LOCATION__, &
                       "The number of constraint atoms exceeds the total number of atoms.")
      in_memory = .FALSE.
      IF (becke_control%save_pot) THEN
         in_memory = becke_control%in_memory
      END IF
      IF (becke_control%save_pot .OR. &
          becke_control%cavity_confine .OR. &
          becke_control%should_skip) THEN
         ALLOCATE (is_constraint(natom))
         is_constraint = .FALSE.
      END IF
      DO i = 1, becke_control%natoms
         ! Notice that here is_constraint=.TRUE. also for dummy atoms to properly compute their Becke charges
         ! A subsequent check (atom_in_group) ensures that the gradients of these dummy atoms are correct
         IF (becke_control%save_pot .OR. &
             becke_control%cavity_confine .OR. &
             becke_control%should_skip) &
            is_constraint(becke_control%atoms(i)) = .TRUE.
      END DO
      eps_cavity = becke_control%eps_cavity
      ! Setup atomic radii for adjusting cell boundaries
      IF (becke_control%adjust) THEN
         IF (.NOT. ASSOCIATED(becke_control%radii)) THEN
            CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
            IF (.NOT. SIZE(atomic_kind_set) == SIZE(becke_control%radii_tmp)) &
               CALL cp_abort(__LOCATION__, &
                             "Length of keyword BECKE_CONSTRAINT\ATOMIC_RADII does not "// &
                             "match number of atomic kinds in the input coordinate file.")
            ALLOCATE (becke_control%radii(SIZE(atomic_kind_set)))
            becke_control%radii(:) = becke_control%radii_tmp(:)
            DEALLOCATE (becke_control%radii_tmp)
         END IF
      END IF
      ! Setup cutoff scheme
      IF (.NOT. ASSOCIATED(becke_control%cutoffs)) THEN
         CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
         ALLOCATE (becke_control%cutoffs(natom))
         SELECT CASE (becke_control%cutoff_type)
         CASE (becke_cutoff_global)
            becke_control%cutoffs(:) = becke_control%rglobal
         CASE (becke_cutoff_element)
            IF (.NOT. SIZE(atomic_kind_set) == SIZE(becke_control%cutoffs_tmp)) &
               CALL cp_abort(__LOCATION__, &
                             "Length of keyword BECKE_CONSTRAINT\ELEMENT_CUTOFFS does not "// &
                             "match number of atomic kinds in the input coordinate file.")
            DO ikind = 1, SIZE(atomic_kind_set)
               CALL get_atomic_kind(atomic_kind_set(ikind), natom=katom, atom_list=atom_list)
               DO iatom = 1, katom
                  atom_a = atom_list(iatom)
                  becke_control%cutoffs(atom_a) = becke_control%cutoffs_tmp(ikind)
               END DO
            END DO
            DEALLOCATE (becke_control%cutoffs_tmp)
         END SELECT
      END IF
      ! Zero weight functions
      DO igroup = 1, SIZE(group)
         group(igroup)%weight%pw%cr3d = 0.0_dp
      END DO
      IF (becke_control%atomic_charges) THEN
         DO iatom = 1, becke_control%natoms
            becke_control%charge(iatom)%pw%cr3d = 0.0_dp
         END DO
      END IF
      ! Allocate storage for cell adjustment coefficients and needed distance vectors
      build = .FALSE.
      IF (becke_control%adjust .AND. .NOT. ASSOCIATED(becke_control%aij)) THEN
         ALLOCATE (becke_control%aij(natom, natom))
         build = .TRUE.
      END IF
      IF (becke_control%vector_buffer%store_vectors) THEN
         ALLOCATE (becke_control%vector_buffer%distances(natom))
         ALLOCATE (becke_control%vector_buffer%distance_vecs(3, natom))
         IF (in_memory) ALLOCATE (becke_control%vector_buffer%pair_dist_vecs(3, natom, natom))
         ALLOCATE (becke_control%vector_buffer%position_vecs(3, natom))
      END IF
      ALLOCATE (becke_control%vector_buffer%R12(natom, natom))
      ! Calculate pairwise distances between each atom pair
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
            ! Set up heteronuclear cell partitioning using user defined radii
            IF (build) THEN
               CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, kind_number=ikind)
               ircov = becke_control%radii(ikind)
               CALL get_atomic_kind(atomic_kind=particle_set(jatom)%atomic_kind, kind_number=ikind)
               jrcov = becke_control%radii(ikind)
               IF (ircov .NE. jrcov) THEN
                  chi = ircov/jrcov
                  uij = (chi-1.0_dp)/(chi+1.0_dp)
                  becke_control%aij(iatom, jatom) = uij/(uij**2-1.0_dp)
                  IF (becke_control%aij(iatom, jatom) .GT. 0.5_dp) THEN
                     becke_control%aij(iatom, jatom) = 0.5_dp
                  ELSE IF (becke_control%aij(iatom, jatom) .LT. -0.5_dp) THEN
                     becke_control%aij(iatom, jatom) = -0.5_dp
                  END IF
               ELSE
                  becke_control%aij(iatom, jatom) = 0.0_dp
               END IF
               ! Note change of sign
               becke_control%aij(jatom, iatom) = -becke_control%aij(iatom, jatom)
            END IF
         END DO
      END DO
      ! Dump some additional information about the calculation
      IF (becke_control%first_iteration) THEN
         IF (iw > 0) THEN
            WRITE (iw, '(/,T3,A)') &
               '----------------------- Becke atomic parameters ------------------------'
            IF (becke_control%adjust) THEN
               WRITE (iw, '(T3,A)') &
                  'Atom   Element           Cutoff (angstrom)        CDFT Radius (angstrom)'
               DO iatom = 1, natom
                  CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, element_symbol=element_symbol, &
                                       kind_number=ikind)
                  ircov = cp_unit_from_cp2k(becke_control%radii(ikind), "angstrom")
                  WRITE (iw, "(i6,T15,A2,T37,F8.3,T67,F8.3)") &
                     iatom, ADJUSTR(element_symbol), cp_unit_from_cp2k(becke_control%cutoffs(iatom), "angstrom"), &
                     ircov
               END DO
            ELSE
               WRITE (iw, '(T3,A)') &
                  'Atom   Element           Cutoff (angstrom)'
               DO iatom = 1, natom
                  CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, element_symbol=element_symbol)
                  WRITE (iw, "(i7,T15,A2,T37,F8.3)") &
                     iatom, ADJUSTR(element_symbol), cp_unit_from_cp2k(becke_control%cutoffs(iatom), "angstrom")
               END DO
            END IF
            WRITE (iw, '(T3,A)') &
               '------------------------------------------------------------------------'
            WRITE (iw, '(/,T3,A,T60)') &
               '----------------------- Becke group definitions ------------------------'
            DO igroup = 1, SIZE(group)
               IF (igroup > 1) WRITE (iw, '(T3,A)') ' '
               WRITE (iw, '(T5,A,I5,A,I5)') &
                  'Atomic group', igroup, ' of ', SIZE(group)
               WRITE (iw, '(T5,A)') 'Atom  Element  Coefficient'
               DO ip = 1, SIZE(group(igroup)%atoms)
                  iatom = group(igroup)%atoms(ip)
                  CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, element_symbol=element_symbol)
                  WRITE (iw, '(i8,T16,A2,T23,F8.3)') iatom, ADJUSTR(element_symbol), group(igroup)%coeff(ip)
               END DO
            END DO
            WRITE (iw, '(T3,A)') &
               '------------------------------------------------------------------------'
         END IF
         becke_control%first_iteration = .FALSE.
      END IF
      ! Setup cavity confinement using spherical Gaussians
      IF (becke_control%cavity_confine) THEN
         cavity_env => becke_control%cavity_env
         CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, pw_env=pw_env, qs_kind_set=qs_kind_set)
         CPASSERT(ASSOCIATED(qs_kind_set))
         nkind = SIZE(qs_kind_set)
         ! Setup the Gaussian shape function
         IF (.NOT. ASSOCIATED(cavity_env%kind_shape_fn)) THEN
            IF (ASSOCIATED(becke_control%radii)) THEN
               ALLOCATE (radii_list(SIZE(becke_control%radii)))
               DO ikind = 1, SIZE(becke_control%radii)
                  IF (cavity_env%use_bohr) THEN
                     radii_list(ikind) = becke_control%radii(ikind)
                  ELSE
                     radii_list(ikind) = cp_unit_from_cp2k(becke_control%radii(ikind), "angstrom")
                  END IF
               END DO
            END IF
            CALL create_shape_function(cavity_env, qs_kind_set, atomic_kind_set, &
                                       radius=becke_control%rcavity, &
                                       radii_list=radii_list)
            IF (ASSOCIATED(radii_list)) &
               DEALLOCATE (radii_list)
         END IF
         ! Form cavity by summing isolated Gaussian densities over constraint atoms
         NULLIFY (rs_cavity)
         CALL pw_env_get(pw_env, auxbas_rs_grid=rs_cavity, auxbas_pw_pool=auxbas_pw_pool)
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
                  IF (becke_control%vector_buffer%store_vectors) THEN
                     ra(:) = becke_control%vector_buffer%position_vecs(:, atom_a)+cell_v(:)/2._dp
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
         CALL pw_pool_create_pw(auxbas_pw_pool, becke_control%cavity%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL rs_pw_transfer(rs_cavity, becke_control%cavity%pw, rs2pw)
         CALL rs_grid_release(rs_cavity)
         ! Grid points where the Gaussian density falls below eps_cavity are ignored
         ! We can calculate the smallest/largest values along z-direction outside
         ! which the cavity is zero at every point (x, y)
         ! If gradients are needed storage needs to be allocated only for grid points within
         ! these bounds
         IF (in_memory .OR. becke_control%save_pot) THEN
            CALL hfun_zero(becke_control%cavity%pw%cr3d, eps_cavity, just_bounds=.TRUE., bounds=bounds)
            ! Save bounds (first nonzero grid point indices)
            bo = group(1)%weight%pw%pw_grid%bounds_local
            IF (bounds(2) .LT. bo(2, 3)) THEN
               bounds(2) = bounds(2)-1
            ELSE
               bounds(2) = bo(2, 3)
            END IF
            IF (bounds(1) .GT. bo(1, 3)) THEN
               ! In the special case bounds(1) == bounds(2) == bo(2, 3), after this check
               ! bounds(1) > bounds(2) and the subsequent gradient allocation (:, :, :, bounds(1):bounds(2))
               ! will correctly allocate a 0-sized array
               bounds(1) = bounds(1)+1
            ELSE
               bounds(1) = bo(1, 3)
            END IF
            becke_control%confine_bounds = bounds
         END IF
         ! Optional printing of cavity (meant for testing, so options currently hardcoded...)
         IF (becke_control%print_cavity) THEN
            CALL hfun_zero(becke_control%cavity%pw%cr3d, eps_cavity, just_bounds=.FALSE.)
            ALLOCATE (stride(3))
            stride = (/2, 2, 2/)
            mpi_io = .TRUE.
            ! Note PROGRAM_RUN_INFO section neeeds to be active!
            unit_nr = cp_print_key_unit_nr(logger, becke_constraint_section, "PROGRAM_RUN_INFO", &
                                           middle_name="BECKE_CAVITY", &
                                           extension=".cube", file_position="REWIND", &
                                           log_filename=.FALSE., mpi_io=mpi_io)
            IF (para_env%mepos == para_env%source .AND. unit_nr .LT. 1) &
               CALL cp_abort(__LOCATION__, &
                             "Please turn on PROGRAM_RUN_INFO to print cavity")
            CALL get_qs_env(qs_env, subsys=subsys)
            CALL qs_subsys_get(subsys, particles=particles)
            CALL cp_pw_to_cube(becke_control%cavity%pw, unit_nr, "CAVITY", particles=particles, stride=stride, mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr, logger, becke_constraint_section, "PROGRAM_RUN_INFO", mpi_io=mpi_io)
            DEALLOCATE (stride)
         END IF
      END IF
      IF (ALLOCATED(is_constraint)) &
         DEALLOCATE (is_constraint)
      CALL timestop(handle)

   END SUBROUTINE becke_constraint_init

! **************************************************************************************************
!> \brief reads the input parameters needed for evaluating a becke weight population constraint
!> \param qs_control the qs_control which holds the Becke control type
!> \param becke_section the input section containing Becke constraint information
!> \author fschiff
! **************************************************************************************************
   SUBROUTINE read_becke_section(qs_control, becke_section)

      TYPE(qs_control_type), INTENT(INOUT)               :: qs_control
      TYPE(section_vals_type), POINTER                   :: becke_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_becke_section', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, jj, k, n_rep, natoms, nvar, &
                                                            tot_natoms
      INTEGER, DIMENSION(:), POINTER                     :: atomlist, dummylist, tmplist
      LOGICAL                                            :: exists, is_duplicate, is_et_coupling
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rtmplist
      TYPE(section_vals_type), POINTER                   :: group_section

      NULLIFY (tmplist, rtmplist, atomlist, dummylist, group_section)

      group_section => section_vals_get_subs_vals(becke_section, "ATOM_GROUP")
      CALL section_vals_get(group_section, n_repetition=nvar, explicit=exists)
      IF (.NOT. exists) CPABORT("Section ATOM_GROUP is missing.")
      ! Special considerations in case ET_COUPLING calculation (in qs_energy_utils)
      is_et_coupling = .FALSE.
      IF (ASSOCIATED(qs_control%becke_control%group)) &
         is_et_coupling = .TRUE.
      IF (is_et_coupling) THEN
         DO i = 1, SIZE(qs_control%becke_control%group)
            DEALLOCATE (qs_control%becke_control%group(i)%atoms)
            DEALLOCATE (qs_control%becke_control%group(i)%coeff)
         END DO
      END IF
      IF (.NOT. is_et_coupling) &
         ALLOCATE (qs_control%becke_control%group(nvar))
      tot_natoms = 0
      ! Parse all ATOM_GROUP sections
      DO k = 1, nvar
         ! First determine how much storage is needed
         natoms = 0
         CALL section_vals_val_get(group_section, "ATOMS", i_rep_section=k, n_rep_val=n_rep)
         DO j = 1, n_rep
            CALL section_vals_val_get(group_section, "ATOMS", i_rep_section=k, i_rep_val=j, i_vals=tmplist)
            IF (SIZE(tmplist) < 1) &
               CPABORT("Each ATOM_GROUP must contain at least 1 atom.")
            natoms = natoms+SIZE(tmplist)
         END DO
         ALLOCATE (qs_control%becke_control%group(k)%atoms(natoms))
         ALLOCATE (qs_control%becke_control%group(k)%coeff(natoms))
         IF (.NOT. is_et_coupling) THEN
            NULLIFY (qs_control%becke_control%group(k)%weight%pw)
            NULLIFY (qs_control%becke_control%group(k)%gradients)
            NULLIFY (qs_control%becke_control%group(k)%integrated)
         END IF
         tot_natoms = tot_natoms+natoms
         ! Now parse
         jj = 0
         DO j = 1, n_rep
            CALL section_vals_val_get(group_section, "ATOMS", i_rep_section=k, i_rep_val=j, i_vals=tmplist)
            DO i = 1, SIZE(tmplist)
               jj = jj+1
               qs_control%becke_control%group(k)%atoms(jj) = tmplist(i)
            END DO
         END DO
         CALL section_vals_val_get(group_section, "COEFF", i_rep_section=k, n_rep_val=n_rep)
         jj = 0
         DO j = 1, n_rep
            CALL section_vals_val_get(group_section, "COEFF", i_rep_section=k, i_rep_val=j, r_vals=rtmplist)
            DO i = 1, SIZE(rtmplist)
               jj = jj+1
               IF (jj > natoms) &
                  CPABORT("Length of keywords ATOMS and COEFF must match.")
               IF (ABS(rtmplist(i)) /= 1.0_dp) &
                  CPABORT("Keyword COEFF accepts only values +/-1.0")
               qs_control%becke_control%group(k)%coeff(jj) = rtmplist(i)
            END DO
         END DO
         IF (jj < natoms) &
            CPABORT("Length of keywords ATOMS and COEFF must match.")
         CALL section_vals_val_get(group_section, "CONSTRAINT_TYPE", i_rep_section=k, &
                                   i_val=qs_control%becke_control%group(k)%constraint_type)
         CALL section_vals_val_get(group_section, "FRAGMENT_CONSTRAINT", i_rep_section=k, &
                                   l_val=qs_control%becke_control%group(k)%is_fragment_constraint)
         IF (qs_control%becke_control%group(k)%is_fragment_constraint) &
            qs_control%becke_control%fragment_density = .TRUE.
      END DO
      ! Create a list containing all constraint atoms
      ALLOCATE (atomlist(tot_natoms))
      atomlist = -1
      jj = 0
      DO k = 1, nvar
         DO j = 1, SIZE(qs_control%becke_control%group(k)%atoms)
            is_duplicate = .FALSE.
            DO i = 1, jj+1
               IF (qs_control%becke_control%group(k)%atoms(j) == atomlist(i)) THEN
                  is_duplicate = .TRUE.
                  EXIT
               END IF
            END DO
            IF (.NOT. is_duplicate) THEN
               jj = jj+1
               atomlist(jj) = qs_control%becke_control%group(k)%atoms(j)
            END IF
         END DO
      END DO
      CALL reallocate(atomlist, 1, jj)
      CALL section_vals_val_get(becke_section, "ATOMIC_CHARGES", &
                                l_val=qs_control%becke_control%atomic_charges)
      ! Parse any dummy atoms (no constraint, just charges)
      IF (qs_control%becke_control%atomic_charges) THEN
         group_section => section_vals_get_subs_vals(becke_section, "DUMMY_ATOMS")
         CALL section_vals_get(group_section, explicit=exists)
         IF (exists) THEN
            ! First determine how many atoms there are
            natoms = 0
            CALL section_vals_val_get(group_section, "ATOMS", n_rep_val=n_rep)
            DO j = 1, n_rep
               CALL section_vals_val_get(group_section, "ATOMS", i_rep_val=j, i_vals=tmplist)
               IF (SIZE(tmplist) < 1) &
                  CPABORT("DUMMY_ATOMS must contain at least 1 atom.")
               natoms = natoms+SIZE(tmplist)
            END DO
            ALLOCATE (dummylist(natoms))
            ! Now parse
            jj = 0
            DO j = 1, n_rep
               CALL section_vals_val_get(group_section, "ATOMS", i_rep_val=j, i_vals=tmplist)
               DO i = 1, SIZE(tmplist)
                  jj = jj+1
                  dummylist(jj) = tmplist(i)
               END DO
            END DO
            ! Check for duplicates
            DO j = 1, natoms
               DO i = j+1, natoms
                  IF (dummylist(i) == dummylist(j)) &
                     CPABORT("Duplicate atoms defined in section DUMMY_ATOMS.")
               END DO
            END DO
            ! Check that a dummy atom is not included in any ATOM_GROUP
            DO j = 1, SIZE(atomlist)
               DO i = 1, SIZE(dummylist)
                  IF (dummylist(i) == atomlist(j)) &
                     CALL cp_abort(__LOCATION__, &
                                   "Duplicate atoms defined in sections ATOM_GROUP and DUMMY_ATOMS.")
               END DO
            END DO
         END IF
      END IF
      ! Join dummy atoms and constraint atoms into one list
      IF (ASSOCIATED(qs_control%becke_control%atoms)) &
         DEALLOCATE (qs_control%becke_control%atoms)
      IF (ASSOCIATED(qs_control%becke_control%is_constraint)) &
         DEALLOCATE (qs_control%becke_control%is_constraint)
      IF (ASSOCIATED(dummylist)) THEN
         qs_control%becke_control%natoms = SIZE(atomlist)+SIZE(dummylist)
      ELSE
         qs_control%becke_control%natoms = SIZE(atomlist)
      END IF
      ALLOCATE (qs_control%becke_control%atoms(qs_control%becke_control%natoms))
      ALLOCATE (qs_control%becke_control%is_constraint(qs_control%becke_control%natoms))
      qs_control%becke_control%atoms(1:SIZE(atomlist)) = atomlist
      IF (ASSOCIATED(dummylist)) THEN
         qs_control%becke_control%atoms(1+SIZE(atomlist):) = dummylist
         DEALLOCATE (dummylist)
      END IF
      qs_control%becke_control%is_constraint = .FALSE.
      qs_control%becke_control%is_constraint(1:SIZE(atomlist)) = .TRUE.
      DEALLOCATE (atomlist)
      ! Need to deallocate first in case ET_COUPLING calculation (in qs_energy_utils)
      IF (ASSOCIATED(qs_control%becke_control%strength)) &
         DEALLOCATE (qs_control%becke_control%strength)
      IF (ASSOCIATED(qs_control%becke_control%becke_order_p)) &
         DEALLOCATE (qs_control%becke_control%becke_order_p)
      IF (ASSOCIATED(qs_control%becke_control%target)) &
         DEALLOCATE (qs_control%becke_control%target)
      !
      ALLOCATE (qs_control%becke_control%strength(nvar))
      ALLOCATE (qs_control%becke_control%becke_order_p(nvar))
      ALLOCATE (qs_control%becke_control%target(nvar))
      CALL section_vals_val_get(becke_section, "STRENGTH", r_vals=rtmplist)
      IF (SIZE(rtmplist) /= nvar) &
         CALL cp_abort(__LOCATION__, &
                       "The length of keyword STRENGTH is incorrect. "// &
                       "Expected "//TRIM(ADJUSTL(cp_to_string(nvar)))// &
                       " value(s), got "// &
                       TRIM(ADJUSTL(cp_to_string(SIZE(rtmplist))))//" value(s).")
      DO j = 1, nvar
         qs_control%becke_control%strength(j) = rtmplist(j)
      END DO
      CALL section_vals_val_get(becke_section, "TARGET", r_vals=rtmplist)
      IF (SIZE(rtmplist) /= nvar) &
         CALL cp_abort(__LOCATION__, &
                       "The length of keyword TARGET is incorrect. "// &
                       "Expected "//TRIM(ADJUSTL(cp_to_string(nvar)))// &
                       " value(s), got "// &
                       TRIM(ADJUSTL(cp_to_string(SIZE(rtmplist))))//" value(s).")
      DO j = 1, nvar
         qs_control%becke_control%target(j) = rtmplist(j)
      END DO
      IF (qs_control%becke_control%atomic_charges) THEN
         IF (ASSOCIATED(qs_control%becke_control%charge)) &
            DEALLOCATE (qs_control%becke_control%charge)
         ALLOCATE (qs_control%becke_control%charge(qs_control%becke_control%natoms))
      END IF

      CALL section_vals_val_get(becke_section, "ADJUST_SIZE", &
                                l_val=qs_control%becke_control%adjust)
      IF (qs_control%becke_control%adjust) THEN
         CALL section_vals_val_get(becke_section, "ATOMIC_RADII", explicit=exists)
         IF (.NOT. exists) CPABORT("Keyword ATOMIC_RADII is missing.")
         CALL section_vals_val_get(becke_section, "ATOMIC_RADII", r_vals=rtmplist)
         CPASSERT(SIZE(rtmplist) > 0)
         ALLOCATE (qs_control%becke_control%radii_tmp(SIZE(rtmplist)))
         DO j = 1, SIZE(rtmplist)
            qs_control%becke_control%radii_tmp(j) = rtmplist(j)
         END DO
      END IF

      CALL section_vals_val_get(becke_section, "CUTOFF_TYPE", &
                                i_val=qs_control%becke_control%cutoff_type)
      SELECT CASE (qs_control%becke_control%cutoff_type)
      CASE (becke_cutoff_global)
         CALL section_vals_val_get(becke_section, "GLOBAL_CUTOFF", &
                                   r_val=qs_control%becke_control%rglobal)
      CASE (becke_cutoff_element)
         CALL section_vals_val_get(becke_section, "ELEMENT_CUTOFF", r_vals=rtmplist)
         CPASSERT(SIZE(rtmplist) > 0)
         ALLOCATE (qs_control%becke_control%cutoffs_tmp(SIZE(rtmplist)))
         DO j = 1, SIZE(rtmplist)
            qs_control%becke_control%cutoffs_tmp(j) = rtmplist(j)
         END DO
      END SELECT
      CALL section_vals_val_get(becke_section, "CAVITY_CONFINE", &
                                l_val=qs_control%becke_control%cavity_confine)
      CALL section_vals_val_get(becke_section, "SHOULD_SKIP", &
                                l_val=qs_control%becke_control%should_skip)
      CALL section_vals_val_get(becke_section, "IN_MEMORY", &
                                l_val=qs_control%becke_control%in_memory)
      IF (qs_control%becke_control%cavity_confine) THEN
         CALL section_vals_val_get(becke_section, "CAVITY_SHAPE", &
                                   i_val=qs_control%becke_control%cavity_shape)
         IF (qs_control%becke_control%cavity_shape == radius_user .AND. &
             .NOT. qs_control%becke_control%adjust) &
            CALL cp_abort(__LOCATION__, &
                          "Activate keyword ADJUST_SIZE to use cavity shape USER.")
         CALL section_vals_val_get(becke_section, "CAVITY_RADIUS", &
                                   r_val=qs_control%becke_control%rcavity)
         CALL section_vals_val_get(becke_section, "EPS_CAVITY", &
                                   r_val=qs_control%becke_control%eps_cavity)
         CALL section_vals_val_get(becke_section, "CAVITY_PRINT", &
                                   l_val=qs_control%becke_control%print_cavity)
         CALL section_vals_val_get(becke_section, "CAVITY_USE_BOHR", &
                                   l_val=qs_control%becke_control%use_bohr)
         IF (.NOT. qs_control%becke_control%use_bohr) THEN
            qs_control%becke_control%rcavity = cp_unit_from_cp2k(qs_control%becke_control%rcavity, "angstrom")
         END IF
         CALL create_hirshfeld_type(qs_control%becke_control%cavity_env)
         CALL set_hirshfeld_info(qs_control%becke_control%cavity_env, &
                                 shape_function_type=shape_function_gaussian, iterative=.FALSE., &
                                 radius_type=qs_control%becke_control%cavity_shape, &
                                 use_bohr=qs_control%becke_control%use_bohr)
      END IF

      IF (qs_control%becke_control%fragment_density) THEN
         CALL section_vals_val_get(becke_section, "FRAGMENT_A_FILE_NAME", &
                                   c_val=qs_control%becke_control%fragment_a_fname)
         CALL section_vals_val_get(becke_section, "FRAGMENT_B_FILE_NAME", &
                                   c_val=qs_control%becke_control%fragment_b_fname)
         CALL section_vals_val_get(becke_section, "FRAGMENT_A_SPIN_FILE", &
                                   c_val=qs_control%becke_control%fragment_a_spin_fname)
         CALL section_vals_val_get(becke_section, "FRAGMENT_B_SPIN_FILE", &
                                   c_val=qs_control%becke_control%fragment_b_spin_fname)
         CALL section_vals_val_get(becke_section, "FLIP_FRAGMENT_A", &
                                   l_val=qs_control%becke_control%flip_fragment(1))
         CALL section_vals_val_get(becke_section, "FLIP_FRAGMENT_B", &
                                   l_val=qs_control%becke_control%flip_fragment(2))
      END IF

      CALL cite_reference(Becke1988b)

   END SUBROUTINE read_becke_section

! **************************************************************************************************
!> \brief reads the input parameters needed for CDFT with OT
!> \param qs_control the qs_control which holds the CDFT control type
!> \param cdft_control_section the input section for CDFT
!> \author Nico Holmberg [12.2015]
! **************************************************************************************************
   SUBROUTINE read_cdft_control_section(qs_control, cdft_control_section)
      TYPE(qs_control_type), INTENT(INOUT)               :: qs_control
      TYPE(section_vals_type), POINTER                   :: cdft_control_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_cdft_control_section', &
         routineP = moduleN//':'//routineN

      INTEGER, DIMENSION(:), POINTER                     :: tmplist
      LOGICAL                                            :: exists
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rtmplist
      TYPE(section_vals_type), POINTER                   :: hirshfeld_constraint_section, &
                                                            outer_scf_section

      NULLIFY (outer_scf_section, hirshfeld_constraint_section)
      CALL section_vals_val_get(cdft_control_section, "TYPE_OF_CONSTRAINT", &
                                i_val=qs_control%cdft_control%type)
      IF (qs_control%cdft_control%type /= outer_scf_none) THEN
         CALL section_vals_val_get(cdft_control_section, "REUSE_PRECOND", &
                                   l_val=qs_control%cdft_control%reuse_precond)
         CALL section_vals_val_get(cdft_control_section, "PRECOND_FREQ", &
                                   i_val=qs_control%cdft_control%precond_freq)
         CALL section_vals_val_get(cdft_control_section, "MAX_REUSE", &
                                   i_val=qs_control%cdft_control%max_reuse)
         CALL section_vals_val_get(cdft_control_section, "PURGE_HISTORY", &
                                   l_val=qs_control%cdft_control%purge_history)
         CALL section_vals_val_get(cdft_control_section, "PURGE_FREQ", &
                                   i_val=qs_control%cdft_control%purge_freq)
         CALL section_vals_val_get(cdft_control_section, "PURGE_OFFSET", &
                                   i_val=qs_control%cdft_control%purge_offset)
         CALL section_vals_val_get(cdft_control_section, "COUNTER", &
                                   i_val=qs_control%cdft_control%ienergy)
         outer_scf_section => section_vals_get_subs_vals(cdft_control_section, "OUTER_SCF")
         CALL section_vals_val_get(outer_scf_section, "_SECTION_PARAMETERS_", &
                                   l_val=qs_control%cdft_control%constraint_control%have_scf)
         IF (qs_control%cdft_control%constraint_control%have_scf) THEN
            CALL section_vals_val_get(outer_scf_section, "TYPE", &
                                      i_val=qs_control%cdft_control%constraint_control%type)

            IF (qs_control%cdft_control%constraint_control%type /= outer_scf_becke_constraint .AND. &
                qs_control%cdft_control%constraint_control%type /= outer_scf_cdft_constraint) &
               CPABORT("Unsupported CDFT constraint.")

            CALL section_vals_val_get(outer_scf_section, "EPS_SCF", &
                                      r_val=qs_control%cdft_control%constraint_control%eps_scf)
            CALL section_vals_val_get(outer_scf_section, "STEP_SIZE", &
                                      r_val=qs_control%cdft_control%constraint_control%step_size, explicit=exists)
            CALL section_vals_val_get(outer_scf_section, "DIIS_BUFFER_LENGTH", &
                                      i_val=qs_control%cdft_control%constraint_control%diis_buffer_length)
            CALL section_vals_val_get(outer_scf_section, "BISECT_TRUST_COUNT", &
                                      i_val=qs_control%cdft_control%constraint_control%bisect_trust_count)
            CALL section_vals_val_get(outer_scf_section, "OPTIMIZER", &
                                      i_val=qs_control%cdft_control%constraint_control%optimizer)
            CALL section_vals_val_get(outer_scf_section, "MAX_SCF", &
                                      i_val=qs_control%cdft_control%constraint_control%max_scf)
            CALL section_vals_val_get(outer_scf_section, "MAX_LS", &
                                      i_val=qs_control%cdft_control%constraint_control%max_ls)
            CALL section_vals_val_get(outer_scf_section, "CONTINUE_LS", &
                                      l_val=qs_control%cdft_control%constraint_control%continue_ls)
            CALL section_vals_val_get(outer_scf_section, "FACTOR_LS", &
                                      r_val=qs_control%cdft_control%constraint_control%factor_ls)
            IF (qs_control%cdft_control%constraint_control%factor_ls .LE. 0.0_dp .OR. &
                qs_control%cdft_control%constraint_control%factor_ls .GE. 1.0_dp) &
               CALL cp_abort(__LOCATION__, &
                             "Illegal value for keyword FACTOR_LS.")
            CALL section_vals_val_get(outer_scf_section, "EXTRAPOLATION_ORDER", &
                                      i_val=qs_control%cdft_control%constraint_control%extrapolation_order)
            CALL section_vals_val_get(outer_scf_section, "JACOBIAN_TYPE", &
                                      i_val=qs_control%cdft_control%constraint_control%jacobian_type)
            CALL section_vals_val_get(outer_scf_section, "JACOBIAN_STEP", &
                                      r_val=qs_control%cdft_control%constraint_control%jacobian_step)
            CALL section_vals_val_get(outer_scf_section, "BROYDEN_TYPE", &
                                      i_val=qs_control%cdft_control%constraint_control%broyden_type)
            IF (exists) THEN
               qs_control%cdft_control%constraint_control%newton_step = &
                  ABS(qs_control%cdft_control%constraint_control%step_size)
               ! Permanent copy needed in case line search is performed
               qs_control%cdft_control%constraint_control%newton_step_save = &
                  qs_control%cdft_control%constraint_control%newton_step
            END IF
            CALL section_vals_val_get(outer_scf_section, "JACOBIAN_FREQ", explicit=exists)
            IF (exists) THEN
               CALL section_vals_val_get(outer_scf_section, "JACOBIAN_FREQ", &
                                         i_vals=tmplist)
               IF (SIZE(tmplist) /= 2) &
                  CALL cp_abort(__LOCATION__, &
                                "Keyword JACOBIAN_FREQ takes exactly two input values.")
               IF (ANY(tmplist .LT. 0)) &
                  CALL cp_abort(__LOCATION__, &
                                "Keyword JACOBIAN_FREQ takes only positive values.")
               IF (ALL(tmplist .EQ. 0)) &
                  CALL cp_abort(__LOCATION__, &
                                "Both values to keyword JACOBIAN_FREQ cannot be zero.")
               qs_control%cdft_control%constraint_control%jacobian_freq(:) = tmplist(1:2)
            END IF
            CALL section_vals_val_get(outer_scf_section, "JACOBIAN_RESTART", &
                                      l_val=qs_control%cdft_control%constraint_control%jacobian_restart)
            SELECT CASE (qs_control%cdft_control%constraint_control%optimizer)
            CASE DEFAULT
               qs_control%cdft_control%constraint_control%jacobian_restart = .FALSE.
            CASE (outer_scf_optimizer_newton, outer_scf_optimizer_newton_ls, outer_scf_optimizer_broyden)
               ! Use parsed value
            END SELECT
            IF (qs_control%cdft_control%constraint_control%jacobian_restart) THEN
               CALL section_vals_val_get(outer_scf_section, "JACOBIAN_VECTOR", &
                                         r_vals=rtmplist)
               ALLOCATE (qs_control%cdft_control%constraint_control%jacobian_vector(SIZE(rtmplist)))
               qs_control%cdft_control%constraint_control%jacobian_vector = rtmplist
            END IF
            SELECT CASE (qs_control%cdft_control%type)
            CASE (outer_scf_hirshfeld_constraint)
               IF (qs_control%cdft_control%constraint_control%type == outer_scf_cdft_constraint) THEN
                  hirshfeld_constraint_section => section_vals_get_subs_vals(cdft_control_section, "HIRSHFELD_CONSTRAINT")
                  CALL section_vals_get(hirshfeld_constraint_section, explicit=exists)
                  IF (exists) THEN
                     CALL read_hirshfeld_constraint_section(qs_control, hirshfeld_constraint_section)
                     ALLOCATE (qs_control%cdft_control%constraint_type(1))
                     qs_control%cdft_control%constraint_type(1) = cdft_charge_constraint
                  ELSE
                     CPABORT("HIRSHFELD_CONSTRAINT section is missing.")
                  END IF
               ELSE
                  CPABORT("Mismatch in defining constraint.")
               END IF
            CASE (outer_scf_becke_constraint)
               ! Do nothing, yet. For now, constraint is defined QS%BECKE_CONSTRAINT
            END SELECT

            CALL cite_reference(Holmberg2017)
         ELSE
            qs_control%cdft = .FALSE.
         END IF
      ELSE
         qs_control%cdft = .FALSE.
      END IF

   END SUBROUTINE read_cdft_control_section

! **************************************************************************************************
!> \brief reads the input parameters needed for Hirshfeld constraint
!> \param qs_control the qs_control which holds the Hirshfeld constraint
!> \param hirshfeld_constraint_section the input section for a Hirshfeld constraint
!> \author Nico Holmberg [12.2015]
! **************************************************************************************************
   SUBROUTINE read_hirshfeld_constraint_section(qs_control, hirshfeld_constraint_section)
      TYPE(qs_control_type), INTENT(INOUT)               :: qs_control
      TYPE(section_vals_type), POINTER                   :: hirshfeld_constraint_section

      CHARACTER(len=*), PARAMETER :: routineN = 'read_hirshfeld_constraint_section', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: j, jj, k, n_rep, nvar, radius_type, &
                                                            refc, shapef
      INTEGER, DIMENSION(:), POINTER                     :: tmplist
      LOGICAL                                            :: do_radius, do_sc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rtmplist

      NULLIFY (tmplist, rtmplist)
      CPASSERT(.NOT. ASSOCIATED(qs_control%cdft_control%hirshfeld_control))
      ALLOCATE (qs_control%cdft_control%hirshfeld_control)

      CALL section_vals_val_get(hirshfeld_constraint_section, "ATOMS", n_rep_val=n_rep)
      jj = 0
      DO k = 1, n_rep
         CALL section_vals_val_get(hirshfeld_constraint_section, "ATOMS", i_rep_val=k, i_vals=tmplist)
         DO j = 1, SIZE(tmplist)
            jj = jj+1
         END DO
      END DO
      IF (jj < 1) CPABORT("Need at least 1 atom to use Hirshfeld constraints.")
      qs_control%cdft_control%hirshfeld_control%natoms = jj
      ALLOCATE (qs_control%cdft_control%hirshfeld_control%atoms( &
                qs_control%cdft_control%hirshfeld_control%natoms))
      jj = 0
      DO k = 1, n_rep
         CALL section_vals_val_get(hirshfeld_constraint_section, "ATOMS", i_rep_val=k, i_vals=tmplist)
         DO j = 1, SIZE(tmplist)
            jj = jj+1
            qs_control%cdft_control%hirshfeld_control%atoms(jj) = tmplist(j)
         END DO
      END DO

      ALLOCATE (qs_control%cdft_control%hirshfeld_control%coeff(qs_control%cdft_control%hirshfeld_control%natoms))
      qs_control%cdft_control%hirshfeld_control%coeff = 1.0_dp

      CALL section_vals_val_get(hirshfeld_constraint_section, "COEFF", n_rep_val=n_rep)
      jj = 0
      DO k = 1, n_rep
         CALL section_vals_val_get(hirshfeld_constraint_section, "COEFF", i_rep_val=k, r_vals=rtmplist)
         DO j = 1, SIZE(rtmplist)
            jj = jj+1
            IF (jj > qs_control%cdft_control%hirshfeld_control%natoms) &
               CPABORT("Need the same number of coeff as there are atoms.")
            qs_control%cdft_control%hirshfeld_control%coeff(jj) = rtmplist(j)
            IF (ABS(rtmplist(j)) /= 1.0_dp) &
               CPABORT("Illegal coefficient. Only -1.0 or 1.0 allowed.")
         END DO
      END DO
      IF (jj < qs_control%cdft_control%hirshfeld_control%natoms .AND. jj .NE. 0) &
         CPABORT("Need no or the same number of coeff as there are atoms.")

      nvar = 1
      ALLOCATE (qs_control%cdft_control%target(nvar))
      ALLOCATE (qs_control%cdft_control%strength(nvar))
      ALLOCATE (qs_control%cdft_control%value(nvar))
      CALL section_vals_val_get(hirshfeld_constraint_section, "STRENGTH", r_vals=rtmplist)
      DO j = 1, nvar
         qs_control%cdft_control%strength(j) = rtmplist(j)
      END DO
      CALL section_vals_val_get(hirshfeld_constraint_section, "TARGET", r_vals=rtmplist)
      DO j = 1, nvar
         qs_control%cdft_control%target(j) = rtmplist(j)
      END DO

      NULLIFY (qs_control%cdft_control%hirshfeld_control%hirshfeld_env)
      CALL create_hirshfeld_type(qs_control%cdft_control%hirshfeld_control%hirshfeld_env)
      CALL section_vals_val_get(hirshfeld_constraint_section, "SELF_CONSISTENT", l_val=do_sc)
      CALL section_vals_val_get(hirshfeld_constraint_section, "USER_RADIUS", l_val=do_radius)
      CALL section_vals_val_get(hirshfeld_constraint_section, "SHAPE_FUNCTION", i_val=shapef)
      CALL section_vals_val_get(hirshfeld_constraint_section, "REFERENCE_CHARGE", i_val=refc)
      IF (do_radius) THEN
         radius_type = radius_user
      ELSE
         radius_type = radius_covalent
      END IF
      CALL set_hirshfeld_info(qs_control%cdft_control%hirshfeld_control%hirshfeld_env, &
                              shape_function_type=shapef, iterative=do_sc, ref_charge=refc, radius_type=radius_type)

   END SUBROUTINE read_hirshfeld_constraint_section

! **************************************************************************************************
!> \brief Calculate fout = fun1/fun2 or fout = fun1*fun2
!> \param fout the output 3D potential
!> \param fun1 the first input 3D potential
!> \param fun2 the second input 3D potential
!> \param divide logical that decides whether to divide or multiply the input potentials
! **************************************************************************************************
   SUBROUTINE hfun_scale(fout, fun1, fun2, divide)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: fout
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: fun1, fun2
      LOGICAL, INTENT(IN)                                :: divide

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

      IF (divide) THEN
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
      ELSE
         DO i3 = 1, n3
            DO i2 = 1, n2
               DO i1 = 1, n1
                  fout(i1, i2, i3) = fun1(i1, i2, i3)*fun2(i1, i2, i3)
               END DO
            END DO
         END DO
      END IF

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

      INTEGER                                            :: i1, i2, i3, lb, n1, n2, n3, nzeroed, &
                                                            nzeroed_inner, ub
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
            IF (just_bounds) nzeroed_inner = 0
            DO i1 = 1, n1
               IF (fun(i1, i2, i3) < th) THEN
                  IF (just_bounds) THEN
                     nzeroed_inner = nzeroed_inner+1
                  ELSE
                     fun(i1, i2, i3) = 0.0_dp
                  END IF
               ELSE
                  IF (just_bounds) EXIT
               END IF
            END DO
            IF (just_bounds) THEN
               IF (nzeroed_inner < n1) EXIT
               nzeroed = nzeroed+nzeroed_inner
            END IF
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

END MODULE qs_cdft_utils
