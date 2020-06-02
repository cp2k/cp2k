!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief This module deals with all the integrals done on local atomic grids in xas_tdp. This is
!>        mostly used to compute the xc kernel matrix elements wrt two RI basis elements (centered
!>        on the same excited atom) <P|fxc(r)|Q>, where the kernel fxc is purely a function of the
!>        ground state density and r. This is also used to compute the SOC matrix elements in the
!>        orbital basis
! **************************************************************************************************
MODULE xas_tdp_atom
   USE ai_contraction_sphi,             ONLY: ab_contract
   USE atom_operators,                  ONLY: calculate_model_potential
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_array_utils,                  ONLY: cp_1d_i_p_type,&
                                              cp_1d_r_p_type,&
                                              cp_2d_r_p_type,&
                                              cp_3d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              qs_control_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_dist2d_to_dist,&
                                              dbcsr_deallocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cp_para_env,                     ONLY: cp_para_env_create
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_complete_redistribute, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, &
        dbcsr_get_block_p, dbcsr_get_stored_coordinates, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_replicate_all, dbcsr_type
   USE dbcsr_tensor_api,                ONLY: dbcsr_t_destroy,&
                                              dbcsr_t_get_block,&
                                              dbcsr_t_iterator_blocks_left,&
                                              dbcsr_t_iterator_next_block,&
                                              dbcsr_t_iterator_start,&
                                              dbcsr_t_iterator_stop,&
                                              dbcsr_t_iterator_type,&
                                              dbcsr_t_type
   USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                              distribution_2d_type
   USE input_constants,                 ONLY: do_potential_id
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                              get_number_of_lebedev_grid,&
                                              init_lebedev_grids,&
                                              lebedev_grid
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE mathconstants,                   ONLY: dfac,&
                                              pi
   USE mathlib,                         ONLY: get_diag
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_comm_split_direct,&
                                              mp_irecv,&
                                              mp_isend,&
                                              mp_sum,&
                                              mp_sync,&
                                              mp_waitall
   USE orbital_pointers,                ONLY: indco,&
                                              indso,&
                                              nco,&
                                              ncoset,&
                                              nso,&
                                              nsoset
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: c_light_au
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                              create_grid_atom,&
                                              grid_atom_type
   USE qs_harmonics_atom,               ONLY: allocate_harmonics_atom,&
                                              create_harmonics_atom,&
                                              get_maxl_CG,&
                                              get_none0_cg_list,&
                                              harmonics_atom_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_overlap,                      ONLY: build_overlap_matrix_simple
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE spherical_harmonics,             ONLY: clebsch_gordon,&
                                              clebsch_gordon_deallocate,&
                                              clebsch_gordon_init
   USE util,                            ONLY: get_limit,&
                                              sort_unique
   USE xas_tdp_integrals,               ONLY: build_xas_tdp_3c_nl,&
                                              build_xas_tdp_ovlp_nl,&
                                              create_pqX_tensor,&
                                              fill_pqx_tensor,&
                                              get_opt_3c_dist2d
   USE xas_tdp_types,                   ONLY: batch_info_type,&
                                              get_proc_batch_sizes,&
                                              release_batch_info,&
                                              xas_atom_env_type,&
                                              xas_tdp_control_type,&
                                              xas_tdp_env_type
   USE xc,                              ONLY: divide_by_norm_drho
   USE xc_atom,                         ONLY: xc_rho_set_atom_update
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_create,&
                                              xc_dset_get_derivative,&
                                              xc_dset_release
   USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                              xc_derivative_type
   USE xc_derivatives,                  ONLY: xc_functionals_eval,&
                                              xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_atom'

   PUBLIC :: init_xas_atom_env, integrate_fxc_atoms, integrate_soc_atoms

CONTAINS

! **************************************************************************************************
!> \brief Initializes a xas_atom_env type given the qs_enxas_atom_env, qs_envv
!> \param xas_atom_env the xas_atom_env to initialize
!> \param xas_tdp_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE init_xas_atom_env(xas_atom_env, xas_tdp_env, xas_tdp_control, qs_env)

      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'init_xas_atom_env', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikind, natom, nex_atoms, &
                                                            nex_kinds, nkind, nspins
      LOGICAL                                            :: do_soc, do_xc
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, particle_set)

      CALL timeset(routineN, handle)

!  Initializing the type
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, natom=natom, particle_set=particle_set)

      nkind = SIZE(qs_kind_set)
      nex_kinds = xas_tdp_env%nex_kinds
      nex_atoms = xas_tdp_env%nex_atoms
      do_xc = xas_tdp_control%do_xc
      do_soc = xas_tdp_control%do_soc
      nspins = 1; IF (xas_tdp_control%do_uks .OR. xas_tdp_control%do_roks) nspins = 2

      xas_atom_env%nspins = nspins
      xas_atom_env%ri_radius = xas_tdp_control%ri_radius
      ALLOCATE (xas_atom_env%grid_atom_set(nkind))
      ALLOCATE (xas_atom_env%harmonics_atom_set(nkind))
      ALLOCATE (xas_atom_env%ri_sphi_so(nkind))
      ALLOCATE (xas_atom_env%orb_sphi_so(nkind))
      IF (do_xc) THEN
         ALLOCATE (xas_atom_env%gr(nkind))
         ALLOCATE (xas_atom_env%ga(nkind))
         ALLOCATE (xas_atom_env%dgr1(nkind))
         ALLOCATE (xas_atom_env%dgr2(nkind))
         ALLOCATE (xas_atom_env%dga1(nkind))
         ALLOCATE (xas_atom_env%dga2(nkind))
      END IF

      xas_atom_env%excited_atoms => xas_tdp_env%ex_atom_indices
      xas_atom_env%excited_kinds => xas_tdp_env%ex_kind_indices

!  Allocate and initialize the atomic grids and harmonics
      CALL init_xas_atom_grid_harmo(xas_atom_env, xas_tdp_control%grid_info, do_xc, qs_env)

!  Compute the contraction coefficients for spherical orbitals
      DO ikind = 1, nkind
         NULLIFY (xas_atom_env%orb_sphi_so(ikind)%array, xas_atom_env%ri_sphi_so(ikind)%array)
         IF (do_soc) THEN
            CALL compute_sphi_so(ikind, "ORB", xas_atom_env%orb_sphi_so(ikind)%array, qs_env)
         END IF
         IF (ANY(xas_atom_env%excited_kinds == ikind) .AND. do_xc) THEN
            CALL compute_sphi_so(ikind, "RI_XAS", xas_atom_env%ri_sphi_so(ikind)%array, qs_env)
         END IF
      END DO !ikind

!  Compute the coefficients to expand the density in the RI_XAS basis, if requested
      IF (do_xc) CALL calculate_density_coeffs(xas_atom_env, qs_env)

      CALL timestop(handle)

   END SUBROUTINE init_xas_atom_env

! **************************************************************************************************
!> \brief Initializes the atomic grids and harmonics for the RI atomic calculations
!> \param xas_atom_env ...
!> \param grid_info ...
!> \param do_xc Whether the xc kernel will ne computed on the atomic grids. If not, the harmonics
!>        are built for the orbital basis for all kinds.
!> \param qs_env ...
!> \note Largely inspired by init_rho_atom subroutine
! **************************************************************************************************
   SUBROUTINE init_xas_atom_grid_harmo(xas_atom_env, grid_info, do_xc, qs_env)

      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      CHARACTER(len=default_string_length), &
         DIMENSION(:, :), POINTER                        :: grid_info
      LOGICAL, INTENT(IN)                                :: do_xc
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'init_xas_atom_grid_harmo', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: kind_name
      INTEGER :: igrid, ikind, il, iso, iso1, iso2, l1, l1l2, l2, la, lc1, lc2, lcleb, ll, llmax, &
         lp, m1, m2, max_s_harm, max_s_set, maxl, maxlgto, maxs, mm, mp, na, nr, quadrature, stat
      REAL(dp)                                           :: kind_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rga
      REAL(dp), DIMENSION(:, :, :), POINTER              :: my_CG
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: tmp_basis
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(qs_control_type), POINTER                     :: qs_control
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (my_CG, qs_kind_set, tmp_basis, grid_atom, harmonics, qs_control, dft_control)

!  Initialization of some integer for the CG coeff generation
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, dft_control=dft_control)
      IF (do_xc) THEN
         CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto, basis_type="RI_XAS")
      ELSE
         CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto, basis_type="ORB")
      END IF
      qs_control => dft_control%qs_control

      !maximum expansion
      llmax = 2*maxlgto
      max_s_harm = nsoset(llmax)
      max_s_set = nsoset(maxlgto)
      lcleb = llmax

!  Allocate and compute the CG coeffs (copied from init_rho_atom)
      CALL clebsch_gordon_init(lcleb)
      CALL reallocate(my_CG, 1, max_s_set, 1, max_s_set, 1, max_s_harm)

      ALLOCATE (rga(lcleb, 2))
      DO lc1 = 0, maxlgto
         DO iso1 = nsoset(lc1 - 1) + 1, nsoset(lc1)
            l1 = indso(1, iso1)
            m1 = indso(2, iso1)
            DO lc2 = 0, maxlgto
               DO iso2 = nsoset(lc2 - 1) + 1, nsoset(lc2)
                  l2 = indso(1, iso2)
                  m2 = indso(2, iso2)
                  CALL clebsch_gordon(l1, m1, l2, m2, rga)
                  IF (l1 + l2 > llmax) THEN
                     l1l2 = llmax
                  ELSE
                     l1l2 = l1 + l2
                  END IF
                  mp = m1 + m2
                  mm = m1 - m2
                  IF (m1*m2 < 0 .OR. (m1*m2 == 0 .AND. (m1 < 0 .OR. m2 < 0))) THEN
                     mp = -ABS(mp)
                     mm = -ABS(mm)
                  ELSE
                     mp = ABS(mp)
                     mm = ABS(mm)
                  END IF
                  DO lp = MOD(l1 + l2, 2), l1l2, 2
                     il = lp/2 + 1
                     IF (ABS(mp) <= lp) THEN
                     IF (mp >= 0) THEN
                        iso = nsoset(lp - 1) + lp + 1 + mp
                     ELSE
                        iso = nsoset(lp - 1) + lp + 1 - ABS(mp)
                     END IF
                     my_CG(iso1, iso2, iso) = rga(il, 1)
                     ENDIF
                     IF (mp /= mm .AND. ABS(mm) <= lp) THEN
                     IF (mm >= 0) THEN
                        iso = nsoset(lp - 1) + lp + 1 + mm
                     ELSE
                        iso = nsoset(lp - 1) + lp + 1 - ABS(mm)
                     END IF
                     my_CG(iso1, iso2, iso) = rga(il, 2)
                     ENDIF
                  END DO
               ENDDO ! iso2
            ENDDO ! lc2
         ENDDO ! iso1
      ENDDO ! lc1
      DEALLOCATE (rga)
      CALL clebsch_gordon_deallocate()

!  Create the Lebedev grids and compute the spherical harmonics
      CALL init_lebedev_grids()
      quadrature = qs_control%gapw_control%quadrature

      DO ikind = 1, SIZE(xas_atom_env%grid_atom_set)

!        Allocate the grid and the harmonics for this kind
         NULLIFY (xas_atom_env%grid_atom_set(ikind)%grid_atom)
         NULLIFY (xas_atom_env%harmonics_atom_set(ikind)%harmonics_atom)
         CALL allocate_grid_atom(xas_atom_env%grid_atom_set(ikind)%grid_atom)
         CALL allocate_harmonics_atom(xas_atom_env%harmonics_atom_set(ikind)%harmonics_atom)

         NULLIFY (grid_atom, harmonics)
         grid_atom => xas_atom_env%grid_atom_set(ikind)%grid_atom
         harmonics => xas_atom_env%harmonics_atom_set(ikind)%harmonics_atom

!        Initialize some integers
         CALL get_qs_kind(qs_kind_set(ikind), ngrid_rad=nr, ngrid_ang=na, name=kind_name)

         !take the grid dimension given as input, if none, take the GAPW ones above
         IF (ANY(grid_info == kind_name)) THEN
            DO igrid = 1, SIZE(grid_info, 1)
               IF (grid_info(igrid, 1) == kind_name) THEN

                  !hack to convert string into integer
                  READ (grid_info(igrid, 2), *, iostat=stat) na
                  IF (stat .NE. 0) CPABORT("The 'na' value for the GRID keyword must be an integer")
                  READ (grid_info(igrid, 3), *, iostat=stat) nr
                  IF (stat .NE. 0) CPABORT("The 'nr' value for the GRID keyword must be an integer")
                  EXIT
               END IF
            END DO
         ELSE IF (do_xc .AND. ANY(xas_atom_env%excited_kinds == ikind)) THEN
            !need good integration grids for the xc kernel, but taking the default GAPW grid
            CPWARN("The default (and possibly small) GAPW grid is being used for one excited KIND")
         END IF

         ll = get_number_of_lebedev_grid(n=na)
         na = lebedev_grid(ll)%n
         la = lebedev_grid(ll)%l
         grid_atom%ng_sphere = na
         grid_atom%nr = nr

!        If this is an excited kind, create the harmonics with the RI_XAS basis, otherwise the ORB
         IF (ANY(xas_atom_env%excited_kinds == ikind) .AND. do_xc) THEN
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=tmp_basis, basis_type="RI_XAS")
         ELSE
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=tmp_basis, basis_type="ORB")
         END IF
         CALL get_gto_basis_set(gto_basis_set=tmp_basis, maxl=maxl, kind_radius=kind_radius)

         CALL create_grid_atom(grid_atom, nr, na, llmax, ll, quadrature)
         CALL truncate_radial_grid(grid_atom, kind_radius)

         maxs = nsoset(maxl)
         CALL create_harmonics_atom(harmonics, &
                                    my_CG, na, llmax, maxs, max_s_harm, ll, grid_atom%wa, &
                                    grid_atom%azi, grid_atom%pol)
         CALL get_maxl_CG(harmonics, tmp_basis, llmax, max_s_harm)
      END DO

      CALL deallocate_lebedev_grids()
      DEALLOCATE (my_CG)

   END SUBROUTINE init_xas_atom_grid_harmo

! **************************************************************************************************
!> \brief Reduces the radial extension of an atomic grid such that it only covers a given radius
!> \param grid_atom the atomic grid from which we truncate the radial part
!> \param max_radius the maximal radial extension of the resulting grid
!> \note Since the RI density used for <P|fxc|Q> is only of quality close to the atom, and the
!>       integrand only non-zero within the radius of the gaussian P,Q. One can reduce the grid
!>       extansion to the largest radius of the RI basis set
! **************************************************************************************************
   SUBROUTINE truncate_radial_grid(grid_atom, max_radius)

      TYPE(grid_atom_type), POINTER                      :: grid_atom
      REAL(dp), INTENT(IN)                               :: max_radius

      CHARACTER(len=*), PARAMETER :: routineN = 'truncate_radial_grid', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: first_ir, ir, llmax_p1, na, new_nr, nr

      nr = grid_atom%nr
      na = grid_atom%ng_sphere
      llmax_p1 = SIZE(grid_atom%rad2l, 2) - 1

!     Find the index corresponding to the limiting radius (small ir => large radius)
      DO ir = 1, nr
         IF (grid_atom%rad(ir) < max_radius) THEN
            first_ir = ir
            EXIT
         END IF
      END DO
      new_nr = nr - first_ir + 1

!     Reallcoate everything that depends on r
      grid_atom%nr = new_nr

      grid_atom%rad(1:new_nr) = grid_atom%rad(first_ir:nr)
      grid_atom%rad2(1:new_nr) = grid_atom%rad2(first_ir:nr)
      grid_atom%wr(1:new_nr) = grid_atom%wr(first_ir:nr)
      grid_atom%weight(:, 1:new_nr) = grid_atom%weight(:, first_ir:nr)
      grid_atom%rad2l(1:new_nr, :) = grid_atom%rad2l(first_ir:nr, :)
      grid_atom%oorad2l(1:new_nr, :) = grid_atom%oorad2l(first_ir:nr, :)

      CALL reallocate(grid_atom%rad, 1, new_nr)
      CALL reallocate(grid_atom%rad2, 1, new_nr)
      CALL reallocate(grid_atom%wr, 1, new_nr)
      CALL reallocate(grid_atom%weight, 1, na, 1, new_nr)
      CALL reallocate(grid_atom%rad2l, 1, new_nr, 0, llmax_p1)
      CALL reallocate(grid_atom%oorad2l, 1, new_nr, 0, llmax_p1)

   END SUBROUTINE truncate_radial_grid

! **************************************************************************************************
!> \brief Computes the contraction coefficients to go from spherical orbitals to sgf for a given
!>        atomic kind and a given basis type.
!> \param ikind the kind for which we compute the coefficients
!> \param basis_type the type of basis for which we compute
!> \param sphi_so where the new contraction coefficients are stored (not yet allocated)
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_sphi_so(ikind, basis_type, sphi_so, qs_env)

      INTEGER, INTENT(IN)                                :: ikind
      CHARACTER(len=*), INTENT(IN)                       :: basis_type
      REAL(dp), DIMENSION(:, :), POINTER                 :: sphi_so
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_sphi_so', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ico, ipgf, iset, iso, l, lx, ly, lz, &
                                                            maxso, nset, sgfi, start_c, start_s
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf, nsgf_set
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf
      REAL(dp)                                           :: factor
      REAL(dp), DIMENSION(:, :), POINTER                 :: sphi
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (basis, lmax, lmin, npgf, nsgf_set, qs_kind_set, first_sgf, sphi)

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis, basis_type=basis_type)
      CALL get_gto_basis_set(basis, lmax=lmax, nset=nset, npgf=npgf, maxso=maxso, lmin=lmin, &
                             nsgf_set=nsgf_set, sphi=sphi, first_sgf=first_sgf)

      ALLOCATE (sphi_so(maxso, SUM(nsgf_set)))
      sphi_so = 0.0_dp

      DO iset = 1, nset
         sgfi = first_sgf(1, iset)
         DO ipgf = 1, npgf(iset)
            start_s = (ipgf - 1)*nsoset(lmax(iset))
            start_c = (ipgf - 1)*ncoset(lmax(iset))

            DO l = lmin(iset), lmax(iset)
               DO iso = 1, nso(l)
                  DO ico = 1, nco(l)

                     lx = indco(1, ico + ncoset(l - 1))
                     ly = indco(2, ico + ncoset(l - 1))
                     lz = indco(3, ico + ncoset(l - 1))

                     factor = orbtramat(l)%s2c(iso, ico) &
                              *SQRT(4.0_dp*pi/dfac(2*l + 1)*dfac(2*lx - 1)*dfac(2*ly - 1)*dfac(2*lz - 1))

                     CALL daxpy(nsgf_set(iset), factor, &
                                sphi(start_c + ncoset(l - 1) + ico, sgfi:sgfi + nsgf_set(iset) - 1), 1, &
                                sphi_so(start_s + nsoset(l - 1) + iso, sgfi:sgfi + nsgf_set(iset) - 1), 1)

                  END DO !ico
               END DO !iso
            END DO !l
         END DO !ipgf
      END DO !iset

   END SUBROUTINE compute_sphi_so

! **************************************************************************************************
!> \brief Find the neighbors of a given set of atoms based on the non-zero blocks of a provided
!>        overlap matrix. Optionally returns an array containing the indices of all involved atoms
!>        (the given subset plus all their neighbors, without repetition) AND/OR an array of arrays
!>        providing the indices of the neighbors of each input atom.
!> \param base_atoms the set of atoms for which we search neighbors
!> \param mat_s the overlap matrix used to find neighbors
!> \param radius the cutoff radius after which atoms are not considered neighbors
!> \param qs_env ...
!> \param all_neighbors the array uniquely contatining all indices of all atoms involved
!> \param neighbor_set array of arrays containing the neighbors of all given atoms
! **************************************************************************************************
   SUBROUTINE find_neighbors(base_atoms, mat_s, radius, qs_env, all_neighbors, neighbor_set)

      INTEGER, DIMENSION(:), INTENT(INOUT)               :: base_atoms
      TYPE(dbcsr_type), INTENT(IN)                       :: mat_s
      REAL(dp)                                           :: radius
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: all_neighbors
      TYPE(cp_1d_i_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: neighbor_set

      CHARACTER(len=*), PARAMETER :: routineN = 'find_neighbors', routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, i, iat, ibase, iblk, jblk, mepos, &
                                                            natom, nb, nbase
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_to_base, inb, who_is_there
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: n_neighbors
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: is_base_atom
      REAL(dp)                                           :: dist2, rad2, ri(3), rij(3), rj(3)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER        :: my_neighbor_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (particle_set, para_env, my_neighbor_set, cell)

      ! Initialization
      CALL get_qs_env(qs_env, para_env=para_env, natom=natom, particle_set=particle_set, cell=cell)
      mepos = para_env%mepos
      nbase = SIZE(base_atoms)
      !work with the neighbor_set structure, see at the end if we keep it
      ALLOCATE (my_neighbor_set(nbase))
      rad2 = radius**2

      ALLOCATE (blk_to_base(natom), is_base_atom(natom))
      blk_to_base = 0; is_base_atom = .FALSE.
      DO ibase = 1, nbase
         blk_to_base(base_atoms(ibase)) = ibase
         is_base_atom(base_atoms(ibase)) = .TRUE.
      END DO

      ! First loop over S => count the number of neighbors
      ALLOCATE (n_neighbors(nbase, 0:para_env%num_pe - 1))
      n_neighbors = 0

      CALL dbcsr_iterator_start(iter, mat_s)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

         !avoid self-neighbors
         IF (iblk == jblk) CYCLE

         !test distance
         ri = pbc(particle_set(iblk)%r, cell)
         rj = pbc(particle_set(jblk)%r, cell)
         rij = pbc(ri, rj, cell)
         dist2 = SUM(rij**2)
         IF (dist2 > rad2) CYCLE

         IF (is_base_atom(iblk)) THEN
            ibase = blk_to_base(iblk)
            n_neighbors(ibase, mepos) = n_neighbors(ibase, mepos) + 1
         END IF
         IF (is_base_atom(jblk)) THEN
            ibase = blk_to_base(jblk)
            n_neighbors(ibase, mepos) = n_neighbors(ibase, mepos) + 1
         END IF

      END DO !iter
      CALL dbcsr_iterator_stop(iter)
      CALL mp_sum(n_neighbors, para_env%group)

      ! Allocate the neighbor_set arrays at the correct length
      DO ibase = 1, nbase
         ALLOCATE (my_neighbor_set(ibase)%array(SUM(n_neighbors(ibase, :))))
         my_neighbor_set(ibase)%array = 0
      END DO

      ! Loop a second time over S, this time fill the neighbors details
      CALL dbcsr_iterator_start(iter, mat_s)
      ALLOCATE (inb(nbase))
      inb = 1
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
         IF (iblk == jblk) CYCLE

         !test distance
         ri = pbc(particle_set(iblk)%r, cell)
         rj = pbc(particle_set(jblk)%r, cell)
         rij = pbc(ri, rj, cell)
         dist2 = SUM(rij**2)
         IF (dist2 > rad2) CYCLE

         IF (is_base_atom(iblk)) THEN
            ibase = blk_to_base(iblk)
            my_neighbor_set(ibase)%array(SUM(n_neighbors(ibase, 0:mepos - 1)) + inb(ibase)) = jblk
            inb(ibase) = inb(ibase) + 1
         END IF
         IF (is_base_atom(jblk)) THEN
            ibase = blk_to_base(jblk)
            my_neighbor_set(ibase)%array(SUM(n_neighbors(ibase, 0:mepos - 1)) + inb(ibase)) = iblk
            inb(ibase) = inb(ibase) + 1
         END IF

      END DO !iter
      CALL dbcsr_iterator_stop(iter)

      ! Make sure the info is shared among the procs
      DO ibase = 1, nbase
         CALL mp_sum(my_neighbor_set(ibase)%array, para_env%group)
      END DO

      ! Gather all indices if asked for it
      IF (PRESENT(all_neighbors)) THEN
         ALLOCATE (who_is_there(natom))
         who_is_there = 0

         DO ibase = 1, nbase
            who_is_there(base_atoms(ibase)) = 1
            DO nb = 1, SIZE(my_neighbor_set(ibase)%array)
               who_is_there(my_neighbor_set(ibase)%array(nb)) = 1
            END DO
         END DO

         ALLOCATE (all_neighbors(natom))
         i = 0
         DO iat = 1, natom
            IF (who_is_there(iat) == 1) THEN
               i = i + 1
               all_neighbors(i) = iat
            END IF
         END DO
         CALL reallocate(all_neighbors, 1, i)
      END IF

      ! If not asked for the neighbor set, deallocate it
      IF (PRESENT(neighbor_set)) THEN
         neighbor_set => my_neighbor_set
      ELSE
         DO ibase = 1, nbase
            DEALLOCATE (my_neighbor_set(ibase)%array)
         END DO
         DEALLOCATE (my_neighbor_set)
      END IF

   END SUBROUTINE find_neighbors

! **************************************************************************************************
!> \brief Returns the RI inverse overlap for a subset of the RI_XAS matrix spaning a given
!>        excited atom and its neighbors.
!> \param ri_sinv the inverse overlap as a dbcsr matrix
!> \param whole_s the whole RI overlap matrix
!> \param neighbors the indeces of the excited atom and their neighbors
!> \param idx_to_nb array telling where any atom can be found in neighbors (if there at all)
!> \param basis_set_ri the RI basis set list for all kinds
!> \param qs_env ...
!> \note It is assumed that the neighbors are sorted, the output matrix is assumed to be small and
!>       is replicated on all processors
! **************************************************************************************************
   SUBROUTINE get_exat_ri_sinv(ri_sinv, whole_s, neighbors, idx_to_nb, basis_set_ri, qs_env)

      TYPE(dbcsr_type)                                   :: ri_sinv, whole_s
      INTEGER, DIMENSION(:), INTENT(IN)                  :: neighbors, idx_to_nb
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_ri
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'get_exat_ri_sinv', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, dest, group, handle, iat, ikind, &
                                                            inb, ir, is, jat, jnb, natom, nnb, &
                                                            npcols, nprows, source, tag
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: recv_req, send_req
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, nsgf, row_dist
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      LOGICAL                                            :: found_risinv, found_whole
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: is_neighbor
      REAL(dp)                                           :: ri(3), rij(3), rj(3)
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: radius
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_risinv, block_whole
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_2d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: recv_buff, send_buff
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type)                      :: sinv_dist
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (pgrid, dbcsr_dist, row_dist, col_dist, nsgf, particle_set, block_whole, block_risinv)
      NULLIFY (cell, para_env, blacs_env)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist, particle_set=particle_set, natom=natom, &
                      para_env=para_env, blacs_env=blacs_env, cell=cell)
      nnb = SIZE(neighbors)
      ALLOCATE (nsgf(nnb), is_neighbor(natom), radius(nnb))
      is_neighbor = .FALSE.
      DO inb = 1, nnb
         iat = neighbors(inb)
         ikind = particle_set(iat)%atomic_kind%kind_number
         CALL get_gto_basis_set(basis_set_ri(ikind)%gto_basis_set, nsgf=nsgf(inb), kind_radius=radius(inb))
         is_neighbor(iat) = .TRUE.
      END DO

      !Create the ri_sinv matrix based on some arbitrary dbcsr_dist
      CALL dbcsr_distribution_get(dbcsr_dist, group=group, pgrid=pgrid, nprows=nprows, npcols=npcols)

      ALLOCATE (row_dist(nnb), col_dist(nnb))
      DO inb = 1, nnb
         row_dist(inb) = MODULO(nprows - inb, nprows)
         col_dist(inb) = MODULO(npcols - inb, npcols)
      END DO

      CALL dbcsr_distribution_new(sinv_dist, group=group, pgrid=pgrid, row_dist=row_dist, &
                                  col_dist=col_dist)

      CALL dbcsr_create(matrix=ri_sinv, name="RI_SINV", matrix_type="S", dist=sinv_dist, &
                        row_blk_size=nsgf, col_blk_size=nsgf)
      !reserving the blocks in the correct pattern
      DO inb = 1, nnb
         ri = pbc(particle_set(neighbors(inb))%r, cell)
         DO jnb = inb, nnb

            !do the atom overlap ?
            rj = pbc(particle_set(neighbors(jnb))%r, cell)
            rij = pbc(ri, rj, cell)
            IF (SUM(rij**2) > (radius(inb) + radius(jnb))**2) CYCLE

            CALL dbcsr_get_stored_coordinates(ri_sinv, inb, jnb, blk)
            IF (para_env%mepos == blk) THEN
               ALLOCATE (block_risinv(nsgf(inb), nsgf(jnb)))
               block_risinv = 0.0_dp
               CALL dbcsr_put_block(ri_sinv, inb, jnb, block_risinv)
               DEALLOCATE (block_risinv)
            END IF
         END DO
      END DO
      CALL dbcsr_finalize(ri_sinv)

      CALL dbcsr_distribution_release(sinv_dist)
      DEALLOCATE (row_dist, col_dist)

      !prepare the send and recv buffers we will need for change of dist between the two matrices
      !worst case scenario: all neighbors are on same procs => need to send nnb**2 messages
      ALLOCATE (send_buff(nnb**2), recv_buff(nnb**2))
      ALLOCATE (send_req(nnb**2), recv_req(nnb**2))
      is = 0; ir = 0

      !Loop over the whole RI overlap matrix and pick the blocks we need
      CALL dbcsr_iterator_start(iter, whole_s)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iat, column=jat, blk=blk)
         CALL dbcsr_get_block_p(whole_s, iat, jat, block_whole, found_whole)

         !only interested in neighbors
         IF (.NOT. found_whole) CYCLE
         IF (.NOT. is_neighbor(iat)) CYCLE
         IF (.NOT. is_neighbor(jat)) CYCLE

         inb = idx_to_nb(iat)
         jnb = idx_to_nb(jat)

         !If blocks are on the same proc for both matrices, simply copy
         CALL dbcsr_get_block_p(ri_sinv, inb, jnb, block_risinv, found_risinv)
         IF (found_risinv) THEN
            CALL dcopy(nsgf(inb)*nsgf(jnb), block_whole, 1, block_risinv, 1)
         ELSE

            !send the block with unique tag to the proc where inb,jnb is in ri_sinv
            CALL dbcsr_get_stored_coordinates(ri_sinv, inb, jnb, dest)
            is = is + 1
            send_buff(is)%array => block_whole
            tag = natom*iat + jat
            CALL mp_isend(msgin=send_buff(is)%array, dest=dest, comm=group, request=send_req(is), tag=tag)

         END IF

      END DO !dbcsr iter
      CALL dbcsr_iterator_stop(iter)

      !Loop over ri_sinv and receive all those blocks
      CALL dbcsr_iterator_start(iter, ri_sinv)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=inb, column=jnb, blk=blk)
         CALL dbcsr_get_block_p(ri_sinv, inb, jnb, block_risinv, found_risinv)

         IF (.NOT. found_risinv) CYCLE

         iat = neighbors(inb)
         jat = neighbors(jnb)

         !If blocks are on the same proc on both matrices do nothing
         CALL dbcsr_get_stored_coordinates(whole_s, iat, jat, source)
         IF (para_env%mepos == source) CYCLE

         tag = natom*iat + jat
         ir = ir + 1
         recv_buff(ir)%array => block_risinv
         CALL mp_irecv(msgout=recv_buff(ir)%array, source=source, request=recv_req(ir), &
                       tag=tag, comm=group)

      END DO
      CALL dbcsr_iterator_stop(iter)

      !make sure that all comm is over before proceeding
      CALL mp_waitall(send_req(1:is))
      CALL mp_waitall(recv_req(1:ir))

      !Invert
      CALL cp_dbcsr_cholesky_decompose(ri_sinv, para_env=para_env, blacs_env=blacs_env)
      CALL cp_dbcsr_cholesky_invert(ri_sinv, para_env=para_env, blacs_env=blacs_env, upper_to_full=.TRUE.)
      CALL dbcsr_filter(ri_sinv, 1.E-10_dp) !make sure ri_sinv is sparse coming out of fm routines
      CALL dbcsr_replicate_all(ri_sinv)

      !clean-up
      DEALLOCATE (nsgf)

      CALL timestop(handle)

   END SUBROUTINE get_exat_ri_sinv

! **************************************************************************************************
!> \brief Compute the coefficients to project the density on a partial RI_XAS basis
!> \param xas_atom_env ...
!> \param qs_env ...
!> \note The density is n = sum_ab P_ab*phi_a*phi_b, the RI basis covers the products of orbital sgfs
!>       => n = sum_ab sum_cd P_ab (phi_a phi_b xi_c) S_cd^-1 xi_d
!>            = sum_d coeff_d xi_d , where xi are the RI basis func.
!>       In this case, with the partial RI projection, the RI basis is restricted to an excited atom
!>       and its neighbors at a time. Leads to smaller overlap matrix to invert and less 3-center
!>       overlap to compute. The procedure is repeated for each excited atom
!>       This is a test implementation for tensors
! **************************************************************************************************
   SUBROUTINE calculate_density_coeffs(xas_atom_env, qs_env)

      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_density_coeffs', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, exat, handle, i, iat, iatom, iex, &
                                                            inb, ind(3), ispin, jatom, jnb, katom, &
                                                            natom, nex, nkind, nnb, nspins, &
                                                            output_unit, ri_at
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_size_orb, blk_size_ri, idx_to_nb, &
                                                            neighbors
      INTEGER, DIMENSION(:), POINTER                     :: all_ri_atoms
      LOGICAL                                            :: pmat_found, pmat_foundt, sinv_found, &
                                                            sinv_foundt, tensor_found, unique
      REAL(dp)                                           :: factor, prefac
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: work2
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: work1
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: t_block
      REAL(dp), DIMENSION(:, :), POINTER                 :: pmat_block, pmat_blockt, sinv_block, &
                                                            sinv_blockt
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type)                      :: opt_dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: overlap, rho_ao, rho_redist
      TYPE(dbcsr_t_iterator_type)                        :: iter
      TYPE(dbcsr_t_type)                                 :: pqX
      TYPE(dbcsr_type)                                   :: ri_sinv
      TYPE(dbcsr_type), POINTER                          :: ri_mats
      TYPE(distribution_2d_type), POINTER                :: opt_dist2d
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_orb, basis_set_ri
      TYPE(libint_potential_type)                        :: pot
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ab_list, ac_list, sab_ri
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho

      NULLIFY (qs_kind_set, basis_set_ri, basis_set_orb, ac_list, rho, rho_ao, sab_ri, ri_mats)
      NULLIFY (particle_set, para_env, all_ri_atoms, overlap, pmat_blockt)
      NULLIFY (blacs_env, pmat_block, ab_list, opt_dist2d, rho_redist, sinv_block, sinv_blockt)

      !Idea: We don't do a full RI here as it would be too expensive in many ways (inversion of a
      !      large matrix, many 3-center overlaps, density coefficients for all atoms, etc...)
      !      Instead, we go excited atom by excited atom and only do a RI expansion on its basis
      !      and that of its closest neighbors (defined through RI_RADIUS), such that we only have
      !      very small matrices to invert and only a few (abc) overlp integrals with c on the
      !      excited atom its neighbors. This is physically sound since we only need the density
      !      well defined on the excited atom grid as we do (aI|P)*(P|Q)^-1*(Q|fxc|R)*(R|S)^-1*(S|Jb)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, nkind=nkind, qs_kind_set=qs_kind_set, rho=rho, &
                      natom=natom, particle_set=particle_set, &
                      para_env=para_env, blacs_env=blacs_env)
      nspins = xas_atom_env%nspins
      nex = SIZE(xas_atom_env%excited_atoms)
      CALL qs_rho_get(rho, rho_ao=rho_ao)

!  Create the needed neighbor list and basis set lists.
      ALLOCATE (basis_set_ri(nkind))
      ALLOCATE (basis_set_orb(nkind))
      CALL basis_set_list_setup(basis_set_ri, "RI_XAS", qs_kind_set)
      CALL basis_set_list_setup(basis_set_orb, "ORB", qs_kind_set)

!  Compute the RI overlap matrix on the whole system
      CALL build_xas_tdp_ovlp_nl(sab_ri, basis_set_ri, basis_set_ri, qs_env)
      CALL build_overlap_matrix_simple(qs_env%ks_env, overlap, basis_set_ri, basis_set_ri, sab_ri)
      ri_mats => overlap(1)%matrix
      CALL release_neighbor_list_sets(sab_ri)

!  Get the neighbors of the excited atoms (= all the atoms where density coeffs are needed)
      CALL find_neighbors(xas_atom_env%excited_atoms, ri_mats, xas_atom_env%ri_radius, &
                          qs_env, all_neighbors=all_ri_atoms, neighbor_set=xas_atom_env%exat_neighbors)

      !keep in mind that double occupation is included in rho_ao in case of closed-shell
      factor = 0.5_dp; IF (nspins == 2) factor = 1.0_dp

!  Allocate space for the projected density coefficients. On all ri_atoms.
!  Note: the sub-region where we project the density changes from excited atom to excited atom
!        => need different sets of RI coeffs
      ALLOCATE (blk_size_ri(natom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=blk_size_ri, basis=basis_set_ri)
      ALLOCATE (xas_atom_env%ri_dcoeff(natom, nspins, nex))
      DO iex = 1, nex
         DO ispin = 1, nspins
            DO iat = 1, natom
               NULLIFY (xas_atom_env%ri_dcoeff(iat, ispin, iex)%array)
               IF ((.NOT. ANY(xas_atom_env%exat_neighbors(iex)%array == iat)) &
                   .AND. (.NOT. xas_atom_env%excited_atoms(iex) == iat)) CYCLE
               ALLOCATE (xas_atom_env%ri_dcoeff(iat, ispin, iex)%array(blk_size_ri(iat)))
               xas_atom_env%ri_dcoeff(iat, ispin, iex)%array = 0.0_dp
            END DO
         END DO
      END DO

      output_unit = cp_logger_get_default_io_unit()
      IF (output_unit > 0) THEN
         WRITE (output_unit, FMT="(/,T7,A,/,T7,A)") &
            "Excited atom, natoms in RI_REGION:", &
            "---------------------------------"
      END IF

      !We go atom by atom, first compute the optimal distribution for the integrals, then
      !the integrals themselves that we put into a tensor, then the contraction with the
      !density

      ALLOCATE (blk_size_orb(natom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=blk_size_orb, basis=basis_set_orb)

      DO iex = 1, nex

         !get neighbors of current atom
         exat = xas_atom_env%excited_atoms(iex)
         nnb = 1 + SIZE(xas_atom_env%exat_neighbors(iex)%array)
         ALLOCATE (neighbors(nnb))
         neighbors(1) = exat
         neighbors(2:nnb) = xas_atom_env%exat_neighbors(iex)%array(:)
         CALL sort_unique(neighbors, unique)

         !link the atoms to their position in neighbors
         ALLOCATE (idx_to_nb(natom))
         idx_to_nb = 0
         DO inb = 1, nnb
            idx_to_nb(neighbors(inb)) = inb
         END DO

         IF (output_unit > 0) THEN
            WRITE (output_unit, FMT="(T7,I12,I21)") &
               exat, nnb
         END IF

         !Get the optimal distribution_2d for the overlap integrals (abc), centers c on the current
         !excited atom and its neighbors defined by RI_RADIUS
         CALL build_xas_tdp_ovlp_nl(ab_list, basis_set_orb, basis_set_orb, qs_env)
         CALL build_xas_tdp_3c_nl(ac_list, basis_set_orb, basis_set_ri, do_potential_id, &
                                  qs_env, excited_atoms=neighbors)
         CALL get_opt_3c_dist2d(opt_dist2d, ab_list, ac_list, basis_set_orb, basis_set_orb, &
                                basis_set_ri, qs_env)
         CALL release_neighbor_list_sets(ab_list)
         CALL release_neighbor_list_sets(ac_list)

         !get the ab and ac neighbor lists based on the optimized distribution
         CALL build_xas_tdp_ovlp_nl(ab_list, basis_set_orb, basis_set_orb, qs_env, ext_dist2d=opt_dist2d)
         CALL build_xas_tdp_3c_nl(ac_list, basis_set_orb, basis_set_ri, do_potential_id, &
                                  qs_env, excited_atoms=neighbors, ext_dist2d=opt_dist2d)

         !get the AO density matrix in the optimized distribution
         CALL cp_dbcsr_dist2d_to_dist(opt_dist2d, opt_dbcsr_dist)
         ALLOCATE (rho_redist(nspins))
         DO ispin = 1, nspins
            ALLOCATE (rho_redist(ispin)%matrix)
            CALL dbcsr_create(rho_redist(ispin)%matrix, template=rho_ao(ispin)%matrix, dist=opt_dbcsr_dist)
            CALL dbcsr_complete_redistribute(rho_ao(ispin)%matrix, rho_redist(ispin)%matrix)
         END DO

         !Compute the 3-center overlap integrals
         pot%potential_type = do_potential_id

         CALL create_pqX_tensor(pqX, ab_list, ac_list, opt_dbcsr_dist, blk_size_orb, blk_size_orb, &
                                blk_size_ri)
         CALL fill_pqX_tensor(pqX, ab_list, ac_list, basis_set_orb, basis_set_orb, basis_set_ri, &
                              pot, qs_env)

         !Compute the RI inverse overlap matrix on the reduced RI basis that spans the excited
         !atom and its neighbors, ri_sinv is replicated over all procs
         CALL get_exat_ri_sinv(ri_sinv, ri_mats, neighbors, idx_to_nb, basis_set_ri, qs_env)

         !Do the actual contraction: coeff_y = sum_pq sum_x P_pq (phi_p phi_q xi_x) S_xy^-1
         !TODO: re-enable OMP parallelization once/if the tensor iterator becomes thread safe

!!$OMP PARALLEL DEFAULT(NONE) &
!!$OMP SHARED (pqX,nspins,blk_size_ri,nnb,neighbors,ri_sinv,factor,xas_atom_env,iex,rho_redist,idx_to_nb)&
!!$OMP PRIVATE(iter,ind,blk,t_block,tensor_found,iatom,jatom,inb,prefac,ispin,work1,pmat_block, &
!!$OMP         pmat_blockt,pmat_found,pmat_foundt,i,jnb,ri_at,sinv_block,sinv_blockt,sinv_found,&
!!$OMP         sinv_foundt,work2,katom)

         CALL dbcsr_t_iterator_start(iter, pqX)
         DO WHILE (dbcsr_t_iterator_blocks_left(iter))
            CALL dbcsr_t_iterator_next_block(iter, ind, blk)
!!$OMP CRITICAL (get_t_block)
            CALL dbcsr_t_get_block(pqX, ind, t_block, tensor_found)
!!$OMP END CRITICAL (get_t_block)

            iatom = ind(1)
            jatom = ind(2)
            katom = ind(3)
            inb = idx_to_nb(katom)

            !non-diagonal elements need to be counted twice
            prefac = 2.0_dp
            IF (iatom == jatom) prefac = 1.0_dp

            DO ispin = 1, nspins

               !rho_redist is symmetric, block can be in either location
               CALL dbcsr_get_block_p(rho_redist(ispin)%matrix, iatom, jatom, pmat_block, pmat_found)
               CALL dbcsr_get_block_p(rho_redist(ispin)%matrix, jatom, iatom, pmat_blockt, pmat_foundt)
               IF ((.NOT. pmat_found) .AND. (.NOT. pmat_foundt)) CYCLE

               ALLOCATE (work1(blk_size_ri(katom), 1))
               work1 = 0.0_dp

               !first contraction with the density matrix
               IF (pmat_found) THEN
                  DO i = 1, blk_size_ri(katom)
                     work1(i, 1) = prefac*SUM(pmat_block(:, :)*t_block(:, :, i))
                  END DO
               ELSE
                  DO i = 1, blk_size_ri(katom)
                     work1(i, 1) = prefac*SUM(TRANSPOSE(pmat_blockt(:, :))*t_block(:, :, i))
                  END DO
               END IF

               !loop over neighbors
               DO jnb = 1, nnb

                  ri_at = neighbors(jnb)

                  !ri_sinv is a symmetric matrix => actual block is one of the two
                  CALL dbcsr_get_block_p(ri_sinv, inb, jnb, sinv_block, sinv_found)
                  CALL dbcsr_get_block_p(ri_sinv, jnb, inb, sinv_blockt, sinv_foundt)
                  IF ((.NOT. sinv_found) .AND. (.NOT. sinv_foundt)) CYCLE

                  !second contraction with the inverse RI overlap
                  ALLOCATE (work2(SIZE(xas_atom_env%ri_dcoeff(ri_at, ispin, iex)%array)))
                  work2 = 0.0_dp

                  IF (sinv_found) THEN
                     DO i = 1, blk_size_ri(katom)
                        CALL daxpy(SIZE(work2), factor*work1(i, 1), sinv_block(i, :), 1, work2(:), 1)
                     END DO
                  ELSE
                     DO i = 1, blk_size_ri(katom)
                        CALL daxpy(SIZE(work2), factor*work1(i, 1), sinv_blockt(:, i), 1, work2(:), 1)
                     END DO
                  END IF

!!$OMP CRITICAL (ri_dcoeff_update)
                  xas_atom_env%ri_dcoeff(ri_at, ispin, iex)%array(:) = &
                     xas_atom_env%ri_dcoeff(ri_at, ispin, iex)%array(:) + work2(:)
!!$OMP END CRITICAL (ri_dcoeff_update)

                  DEALLOCATE (work2)
               END DO !jnb

               DEALLOCATE (work1)
            END DO

            DEALLOCATE (t_block)
         END DO !iter
         CALL dbcsr_t_iterator_stop(iter)
!!$OMP END PARALLEL

         !clean-up
         CALL dbcsr_release(ri_sinv)
         CALL dbcsr_t_destroy(pqX)
         CALL distribution_2d_release(opt_dist2d)
         CALL dbcsr_deallocate_matrix_set(rho_redist)
         CALL dbcsr_distribution_release(opt_dbcsr_dist)
         CALL release_neighbor_list_sets(ab_list)
         CALL release_neighbor_list_sets(ac_list)
         DEALLOCATE (neighbors, idx_to_nb)

      END DO !iex

      !making sure all procs have the same info
      DO iex = 1, nex
         DO ispin = 1, nspins
            DO iat = 1, natom
               IF ((.NOT. ANY(xas_atom_env%exat_neighbors(iex)%array == iat)) &
                   .AND. (.NOT. xas_atom_env%excited_atoms(iex) == iat)) CYCLE
               CALL mp_sum(xas_atom_env%ri_dcoeff(iat, ispin, iex)%array, para_env%group)
            END DO !iat
         END DO !ispin
      END DO !iex

!  clean-up
      CALL dbcsr_deallocate_matrix_set(overlap)
      DEALLOCATE (basis_set_ri, basis_set_orb, all_ri_atoms)

      CALL timestop(handle)

   END SUBROUTINE calculate_density_coeffs

! **************************************************************************************************
!> \brief Evaluates the density on a given atomic grid
!> \param rho_set where the densities are stored
!> \param ri_dcoeff the arrays containing the RI density coefficients of this atom, for each spin
!> \param atom_kind the kind of the atom in question
!> \param do_gga whether the gradient of the density should also be put on the grid
!> \param batch_info how the so are distributed
!> \param xas_atom_env ...
!> \param qs_env ...
!> \note The density is expressed as n = sum_d coeff_d*xi_d. Knowing the coordinate of each grid
!>       grid point, one can simply evaluate xi_d(r)
! **************************************************************************************************
   SUBROUTINE put_density_on_atomic_grid(rho_set, ri_dcoeff, atom_kind, do_gga, batch_info, &
                                         xas_atom_env, qs_env)

      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      TYPE(cp_1d_r_p_type), DIMENSION(:)                 :: ri_dcoeff
      INTEGER, INTENT(IN)                                :: atom_kind
      LOGICAL, INTENT(IN)                                :: do_gga
      TYPE(batch_info_type)                              :: batch_info
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'put_density_on_atomic_grid', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: dir, handle, ipgf, iset, iso, iso_proc, &
                                                            ispin, maxso, n, na, nr, nset, nsgfi, &
                                                            nsoi, nspins, sgfi, starti
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf, nsgf_set
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: so, work1, work2
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: dso
      REAL(dp), DIMENSION(:, :), POINTER                 :: dgr1, dgr2, ga, gr, ri_sphi_so, zet
      REAL(dp), DIMENSION(:, :, :), POINTER              :: dga1, dga2, rhoa, rhob
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: ri_dcoeff_so
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (grid_atom, ri_basis, qs_kind_set, lmax, npgf, zet, nsgf_set, ri_sphi_so)
      NULLIFY (lmin, first_sgf, rhoa, rhob, ga, gr, dgr1, dgr2, dga1, dga2, ri_dcoeff_so)

      CALL timeset(routineN, handle)

!  Strategy: it makes sense to evaluate the spherical orbital on the grid (because of symmetry)
!            From there, one can directly contract into sgf using ri_sphi_so and then take the weight
!            The spherical orbital were precomputed and split in a purely radial and a purely
!            angular part. The full values on each grid point are obtain through gemm

!  Generalities
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(atom_kind), basis_set=ri_basis, basis_type="RI_XAS")
      CALL get_gto_basis_set(ri_basis, lmax=lmax, npgf=npgf, zet=zet, nset=nset, nsgf_set=nsgf_set, &
                             first_sgf=first_sgf, lmin=lmin, maxso=maxso)

!  Get the grid and the info we need from it
      grid_atom => xas_atom_env%grid_atom_set(atom_kind)%grid_atom
      na = grid_atom%ng_sphere
      nr = grid_atom%nr
      n = na*nr
      nspins = xas_atom_env%nspins
      ri_sphi_so => xas_atom_env%ri_sphi_so(atom_kind)%array

!  Point to the rho_set densities
      rhoa => rho_set%rhoa
      rhob => rho_set%rhob
      rhoa = 0.0_dp; rhob = 0.0_dp; 
      IF (do_gga) THEN
         DO dir = 1, 3
            rho_set%drhoa(dir)%array = 0.0_dp
            rho_set%drhob(dir)%array = 0.0_dp
         END DO
      END IF

!  Point to the precomputed SO
      ga => xas_atom_env%ga(atom_kind)%array
      gr => xas_atom_env%gr(atom_kind)%array
      IF (do_gga) THEN
         dga1 => xas_atom_env%dga1(atom_kind)%array
         dga2 => xas_atom_env%dga2(atom_kind)%array
         dgr1 => xas_atom_env%dgr1(atom_kind)%array
         dgr2 => xas_atom_env%dgr2(atom_kind)%array
      END IF

!  Need to express the ri_dcoeffs in terms of so (and not sgf)
      ALLOCATE (ri_dcoeff_so(nspins))
      DO ispin = 1, nspins
         ALLOCATE (ri_dcoeff_so(ispin)%array(nset*maxso))
         ri_dcoeff_so(ispin)%array = 0.0_dp

         !for a given so, loop over sgf and sum
         DO iset = 1, nset
            sgfi = first_sgf(1, iset)
            nsoi = npgf(iset)*nsoset(lmax(iset))
            nsgfi = nsgf_set(iset)
            ALLOCATE (work1(nsoi, 1), work2(nsgfi, 1))

            work2(:, 1) = ri_dcoeff(ispin)%array(sgfi:sgfi + nsgfi - 1)
            CALL dgemm('N', 'N', nsoi, 1, nsgfi, 1.0_dp, ri_sphi_so(1:nsoi, sgfi:sgfi + nsgfi - 1), &
                       nsoi, work2, nsgfi, 0.0_dp, work1, nsoi)
            ri_dcoeff_so(ispin)%array((iset - 1)*maxso + 1:(iset - 1)*maxso + nsoi) = work1(:, 1)

            DEALLOCATE (work1, work2)
         END DO

      END DO

      !allocate space to store the spherical orbitals on the grid
      ALLOCATE (so(na, nr))
      IF (do_gga) ALLOCATE (dso(na, nr, 3))

!  Loop over the spherical orbitals on this proc
      DO iso_proc = 1, batch_info%nso_proc(atom_kind)
         iset = batch_info%so_proc_info(atom_kind)%array(1, iso_proc)
         ipgf = batch_info%so_proc_info(atom_kind)%array(2, iso_proc)
         iso = batch_info%so_proc_info(atom_kind)%array(3, iso_proc)
         IF (iso < 0) CYCLE

         starti = (iset - 1)*maxso + (ipgf - 1)*nsoset(lmax(iset))

         !the spherical orbital on the grid
         CALL dgemm('N', 'T', na, nr, 1, 1.0_dp, ga(:, iso_proc:iso_proc), na, &
                    gr(:, iso_proc:iso_proc), nr, 0.0_dp, so(:, :), na)

         !the gradient on the grid
         IF (do_gga) THEN

            DO dir = 1, 3
               CALL dgemm('N', 'T', na, nr, 1, 1.0_dp, dga1(:, iso_proc:iso_proc, dir), na, &
                          dgr1(:, iso_proc:iso_proc), nr, 0.0_dp, dso(:, :, dir), na)
               CALL dgemm('N', 'T', na, nr, 1, 1.0_dp, dga2(:, iso_proc:iso_proc, dir), na, &
                          dgr2(:, iso_proc:iso_proc), nr, 1.0_dp, dso(:, :, dir), na)
            END DO
         END IF

         !put the so on the grid with the approriate coefficients and sum
         CALL daxpy(n, ri_dcoeff_so(1)%array(starti + iso), so(:, :), 1, rhoa(:, :, 1), 1)

         IF (nspins == 2) THEN
            CALL daxpy(n, ri_dcoeff_so(2)%array(starti + iso), so(:, :), 1, rhob(:, :, 1), 1)
         END IF

         IF (do_gga) THEN

            !put the gradient of the so on the grid with correspond RI coeff
            DO dir = 1, 3
               CALL daxpy(n, ri_dcoeff_so(1)%array(starti + iso), dso(:, :, dir), &
                          1, rho_set%drhoa(dir)%array(:, :, 1), 1)

               IF (nspins == 2) THEN
                  CALL daxpy(n, ri_dcoeff_so(2)%array(starti + iso), dso(:, :, dir), &
                             1, rho_set%drhob(dir)%array(:, :, 1), 1)
               END IF
            END DO !dir
         END IF !do_gga

      END DO

! Treat spin restricted case (=> copy alpha into beta)
      IF (nspins == 1) THEN
         CALL dcopy(n, rhoa(:, :, 1), 1, rhob(:, :, 1), 1)

         IF (do_gga) THEN
            DO dir = 1, 3
               CALL dcopy(n, rho_set%drhoa(dir)%array(:, :, 1), 1, rho_set%drhob(dir)%array(:, :, 1), 1)
            END DO
         END IF
      END IF

! Note: sum over procs is done outside

!  clean-up
      DO ispin = 1, nspins
         DEALLOCATE (ri_dcoeff_so(ispin)%array)
      END DO
      DEALLOCATE (ri_dcoeff_so)

      CALL timestop(handle)

   END SUBROUTINE put_density_on_atomic_grid

! **************************************************************************************************
!> \brief Adds the density of a given source atom with source kind (with ri_dcoeff) on the atomic
!>        grid belonging to another target atom of target kind. The evaluations of the basis
!>        function first requires the evaluation of the x,y,z coordinates on each grid point of
!>        target atom wrt to the position of source atom
!> \param rho_set where the densities are stored
!> \param ri_dcoeff the arrays containing the RI density coefficient of source_iat, for each spin
!> \param source_iat the index of the source atom
!> \param source_ikind the kind of the source atom
!> \param target_iat the index of the target atom
!> \param target_ikind the kind of the target atom
!> \param sr starting r index for the local grid
!> \param er ending r index for the local grid
!> \param do_gga whether the gradient of the density is needed
!> \param xas_atom_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE put_density_on_other_grid(rho_set, ri_dcoeff, source_iat, source_ikind, target_iat, &
                                        target_ikind, sr, er, do_gga, xas_atom_env, qs_env)

      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      TYPE(cp_1d_r_p_type), DIMENSION(:)                 :: ri_dcoeff
      INTEGER, INTENT(IN)                                :: source_iat, source_ikind, target_iat, &
                                                            target_ikind, sr, er
      LOGICAL, INTENT(IN)                                :: do_gga
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'put_density_on_other_grid', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: dir, handle, ia, ico, ipgf, ir, iset, &
                                                            isgf, lx, ly, lz, n, na, nr, nset, &
                                                            nspins, sgfi, start
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf, nsgf_set
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf
      REAL(dp)                                           :: rmom
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: sgf
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: co, dsgf, pos
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :, :)       :: dco
      REAL(dp), DIMENSION(3)                             :: rs, rst, rt
      REAL(dp), DIMENSION(:, :), POINTER                 :: ri_sphi, zet
      REAL(dp), DIMENSION(:, :, :), POINTER              :: rhoa, rhob
      TYPE(cell_type), POINTER                           :: cell
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, ri_basis, lmax, npgf, nsgf_set, lmin, zet, first_sgf, grid_atom)
      NULLIFY (harmonics, rhoa, rhob, particle_set, cell, ri_sphi)

      !Same logic as the  put_density_on_own_grid routine. Loop over orbitals, put them on the grid,
      !contract into sgf and daxpy with coeff. Notable difference: use cartesian orbitals instead of
      !spherical, since the center of the grid is not the origin and thus, spherical symmetry can't
      !be exploited so well

      CALL timeset(routineN, handle)

      !Generalities
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, cell=cell)
      !want basis of the source atom
      CALL get_qs_kind(qs_kind_set(source_ikind), basis_set=ri_basis, basis_type="RI_XAS")
      CALL get_gto_basis_set(ri_basis, lmax=lmax, npgf=npgf, zet=zet, nset=nset, nsgf_set=nsgf_set, &
                             first_sgf=first_sgf, lmin=lmin, sphi=ri_sphi)

      ! Want the grid and harmonics of the target atom
      grid_atom => xas_atom_env%grid_atom_set(target_ikind)%grid_atom
      harmonics => xas_atom_env%harmonics_atom_set(target_ikind)%harmonics_atom
      na = grid_atom%ng_sphere
      nr = er - sr + 1
      n = na*nr
      nspins = xas_atom_env%nspins

      !  Point to the rho_set densities
      rhoa => rho_set%rhoa
      rhob => rho_set%rhob

      !  Need the source-target position vector
      rs = pbc(particle_set(source_iat)%r, cell)
      rt = pbc(particle_set(target_iat)%r, cell)
      rst = pbc(rs, rt, cell)

      ! Precompute the positions on the target grid
      ALLOCATE (pos(na, sr:er, 4))
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(NONE), &
!$OMP SHARED(na,sr,er,pos,harmonics,grid_atom,rst), &
!$OMP PRIVATE(ia,ir)
      DO ir = sr, er
         DO ia = 1, na
            pos(ia, ir, 1:3) = harmonics%a(:, ia)*grid_atom%rad(ir) + rst
            pos(ia, ir, 4) = pos(ia, ir, 1)**2 + pos(ia, ir, 2)**2 + pos(ia, ir, 3)**2
         END DO
      END DO
!$OMP END PARALLEL DO

      ! Loop over the cartesian gaussian functions and evaluate them
      DO iset = 1, nset

         !allocate space to store the cartesian orbtial on the grid
         ALLOCATE (co(na, sr:er, npgf(iset)*ncoset(lmax(iset))))
         IF (do_gga) ALLOCATE (dco(na, sr:er, 3, npgf(iset)*ncoset(lmax(iset))))

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP SHARED(co,npgf,ncoset,lmax,lmin,indco,pos,zet,iset,na,sr,er,do_gga,dco), &
!$OMP PRIVATE(ipgf,start,ico,lx,ly,lz,ia,ir,rmom)

!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
         DO ir = sr, er
            DO ia = 1, na
               co(ia, ir, :) = 0.0_dp
               IF (do_gga) THEN
                  dco(ia, ir, :, :) = 0.0_dp
               END IF
            END DO
         END DO
!$OMP END DO NOWAIT

         DO ipgf = 1, npgf(iset)
            start = (ipgf - 1)*ncoset(lmax(iset))

            !loop over the cartesian orbitals
            DO ico = ncoset(lmin(iset) - 1) + 1, ncoset(lmax(iset))
               lx = indco(1, ico)
               ly = indco(2, ico)
               lz = indco(3, ico)

               ! compute g = x**lx * y**ly * z**lz * exp(-zet * r**2)
!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
               DO ir = sr, er
                  DO ia = 1, na
                     rmom = EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                     IF (lx /= 0) rmom = rmom*pos(ia, ir, 1)**lx
                     IF (ly /= 0) rmom = rmom*pos(ia, ir, 2)**ly
                     IF (lz /= 0) rmom = rmom*pos(ia, ir, 3)**lz
                     co(ia, ir, start + ico) = rmom
                     !co(ia, ir, start + ico) = pos(ia, ir, 1)**lx*pos(ia, ir, 2)**ly*pos(ia, ir, 3)**lz &
                     !                          *EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                  END DO
               END DO
!$OMP END DO NOWAIT

               IF (do_gga) THEN
                  !the gradient: dg_x = lx*x**(lx-1) * y**ly * z**lz * exp(-zet * r**2)
                  !                     -2*zet* x**(lx+1) * y**ly * z**lz * exp(-zet * r**2)
                  !                   = (lx*x**(lx-1) - 2*zet*x**(lx+1)) * y**ly * z**lz * exp(-zet * r**2)

                  !x direction, special case if lx == 0
                  IF (lx == 0) THEN
!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
                     DO ir = sr, er
                        DO ia = 1, na
                           rmom = -2.0_dp*pos(ia, ir, 1)*zet(ipgf, iset)*EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                           IF (ly /= 0) rmom = rmom*pos(ia, ir, 2)**ly
                           IF (lz /= 0) rmom = rmom*pos(ia, ir, 3)**lz
                           dco(ia, ir, 1, start + ico) = rmom
!                          dco(ia, ir, 1, start + ico) = -2.0_dp*pos(ia, ir, 1)*zet(ipgf, iset) &
!                                                        *pos(ia, ir, 2)**ly*pos(ia, ir, 3)**lz &
!                                                        *EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                        END DO
                     END DO
!$OMP END DO NOWAIT
                  ELSE
!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
                     DO ir = sr, er
                        DO ia = 1, na
                           IF (lx /= 1) THEN
                              rmom = (lx*pos(ia, ir, 1)**(lx - 1) - 2.0_dp*pos(ia, ir, 1)**(lx + 1)* &
                                      zet(ipgf, iset))*EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                           ELSE
                              rmom = (1.0_dp - 2.0_dp*pos(ia, ir, 1)**2*zet(ipgf, iset))* &
                                     EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                           END IF
                           IF (ly /= 0) rmom = rmom*pos(ia, ir, 2)**ly
                           IF (lz /= 0) rmom = rmom*pos(ia, ir, 3)**lz
                           dco(ia, ir, 1, start + ico) = rmom
!                          dco(ia, ir, 1, start + ico) = (lx*pos(ia, ir, 1)**(lx - 1) &
!                                                         - 2.0_dp*pos(ia, ir, 1)**(lx + 1)*zet(ipgf, iset)) &
!                                                        *pos(ia, ir, 2)**ly*pos(ia, ir, 3)**lz &
!                                                        *EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                        END DO
                     END DO
!$OMP END DO NOWAIT
                  END IF !lx == 0

                  !y direction, special case if ly == 0
                  IF (ly == 0) THEN
!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
                     DO ir = sr, er
                        DO ia = 1, na
                           rmom = -2.0_dp*pos(ia, ir, 2)*zet(ipgf, iset)*EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                           IF (lx /= 0) rmom = rmom*pos(ia, ir, 1)**lx
                           IF (lz /= 0) rmom = rmom*pos(ia, ir, 3)**lz
                           dco(ia, ir, 2, start + ico) = rmom
!                          dco(ia, ir, 2, start + ico) = -2.0_dp*pos(ia, ir, 2)*zet(ipgf, iset) &
!                                                        *pos(ia, ir, 1)**lx*pos(ia, ir, 3)**lz &
!                                                        *EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                        END DO
                     END DO
!$OMP END DO NOWAIT
                  ELSE
!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
                     DO ir = sr, er
                        DO ia = 1, na
                           IF (ly /= 1) THEN
                              rmom = (ly*pos(ia, ir, 2)**(ly - 1) - 2.0_dp*pos(ia, ir, 2)**(ly + 1)*zet(ipgf, iset)) &
                                     *EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                           ELSE
                              rmom = (1.0_dp - 2.0_dp*pos(ia, ir, 2)**2*zet(ipgf, iset)) &
                                     *EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                           END IF
                           IF (lx /= 0) rmom = rmom*pos(ia, ir, 1)**lx
                           IF (lz /= 0) rmom = rmom*pos(ia, ir, 3)**lz
                           dco(ia, ir, 2, start + ico) = rmom
!                          dco(ia, ir, 2, start + ico) = (ly*pos(ia, ir, 2)**(ly - 1) &
!                                                         - 2.0_dp*pos(ia, ir, 2)**(ly + 1)*zet(ipgf, iset)) &
!                                                        *pos(ia, ir, 1)**lx*pos(ia, ir, 3)**lz &
!                                                        *EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                        END DO
                     END DO
!$OMP END DO NOWAIT
                  END IF !ly == 0

                  !z direction, special case if lz == 0
                  IF (lz == 0) THEN
!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
                     DO ir = sr, er
                        DO ia = 1, na
                           rmom = -2.0_dp*pos(ia, ir, 3)*zet(ipgf, iset)*EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                           IF (lx /= 0) rmom = rmom*pos(ia, ir, 1)**lx
                           IF (ly /= 0) rmom = rmom*pos(ia, ir, 2)**ly
                           dco(ia, ir, 3, start + ico) = rmom
!                          dco(ia, ir, 3, start + ico) = -2.0_dp*pos(ia, ir, 3)*zet(ipgf, iset) &
!                                                        *pos(ia, ir, 1)**lx*pos(ia, ir, 2)**ly &
!                                                        *EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                        END DO
                     END DO
!$OMP END DO NOWAIT
                  ELSE
!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
                     DO ir = sr, er
                        DO ia = 1, na
                           IF (lz /= 1) THEN
                              rmom = (lz*pos(ia, ir, 3)**(lz - 1) - 2.0_dp*pos(ia, ir, 3)**(lz + 1)* &
                                      zet(ipgf, iset))*EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                           ELSE
                              rmom = (1.0_dp - 2.0_dp*pos(ia, ir, 3)**2*zet(ipgf, iset))* &
                                     EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                           END IF
                           IF (lx /= 0) rmom = rmom*pos(ia, ir, 1)**lx
                           IF (ly /= 0) rmom = rmom*pos(ia, ir, 2)**ly
                           dco(ia, ir, 3, start + ico) = rmom
!                          dco(ia, ir, 3, start + ico) = (lz*pos(ia, ir, 3)**(lz - 1) &
!                                                         - 2.0_dp*pos(ia, ir, 3)**(lz + 1)*zet(ipgf, iset)) &
!                                                        *pos(ia, ir, 1)**lx*pos(ia, ir, 2)**ly &
!                                                        *EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                        END DO
                     END DO
!$OMP END DO NOWAIT
                  END IF !lz == 0

               END IF !gga

            END DO !ico
         END DO !ipgf

!$OMP END PARALLEL

         !contract the co into sgf
         ALLOCATE (sgf(na, sr:er))
         IF (do_gga) ALLOCATE (dsgf(na, sr:er, 3))
         sgfi = first_sgf(1, iset) - 1

         DO isgf = 1, nsgf_set(iset)
            sgf = 0.0_dp
            IF (do_gga) dsgf = 0.0_dp

            DO ipgf = 1, npgf(iset)
               start = (ipgf - 1)*ncoset(lmax(iset))
               DO ico = ncoset(lmin(iset) - 1) + 1, ncoset(lmax(iset))
                  CALL daxpy(n, ri_sphi(start + ico, sgfi + isgf), co(:, sr:er, start + ico), 1, sgf(:, sr:er), 1)
               END DO !ico
            END DO !ipgf

            !add the density to the grid
            CALL daxpy(n, ri_dcoeff(1)%array(sgfi + isgf), sgf(:, sr:er), 1, rhoa(:, sr:er, 1), 1)

            IF (nspins == 2) THEN
               CALL daxpy(n, ri_dcoeff(2)%array(sgfi + isgf), sgf(:, sr:er), 1, rhob(:, sr:er, 1), 1)
            END IF

            !deal with the gradient
            IF (do_gga) THEN

               DO ipgf = 1, npgf(iset)
                  start = (ipgf - 1)*ncoset(lmax(iset))
                  DO ico = ncoset(lmin(iset) - 1) + 1, ncoset(lmax(iset))
                     DO dir = 1, 3
                        CALL daxpy(n, ri_sphi(start + ico, sgfi + isgf), dco(:, sr:er, dir, start + ico), &
                                   1, dsgf(:, sr:er, dir), 1)
                     END DO
                  END DO !ico
               END DO !ipgf

               DO dir = 1, 3
                  CALL daxpy(n, ri_dcoeff(1)%array(sgfi + isgf), dsgf(:, sr:er, dir), 1, &
                             rho_set%drhoa(dir)%array(:, sr:er, 1), 1)

                  IF (nspins == 2) THEN
                     CALL daxpy(n, ri_dcoeff(2)%array(sgfi + isgf), dsgf(:, sr:er, dir), 1, &
                                rho_set%drhob(dir)%array(:, sr:er, 1), 1)
                  END IF
               END DO
            END IF !do_gga

         END DO !isgf

         DEALLOCATE (co, sgf)
         IF (do_gga) DEALLOCATE (dco, dsgf)
      END DO !iset

      !Treat spin-restricted case (copy alpha into beta)
      IF (nspins == 1) THEN
         CALL dcopy(n, rhoa(:, sr:er, 1), 1, rhob(:, sr:er, 1), 1)

         IF (do_gga) THEN
            DO dir = 1, 3
               CALL dcopy(n, rho_set%drhoa(dir)%array(:, sr:er, 1), 1, rho_set%drhob(dir)%array(:, sr:er, 1), 1)
            END DO
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE put_density_on_other_grid

! **************************************************************************************************
!> \brief Computes the norm of the density gradient on the atomic grid
!> \param rho_set ...
!> \param atom_kind ...
!> \param xas_atom_env ...
!> \note GGA is assumed
! **************************************************************************************************
   SUBROUTINE compute_norm_drho(rho_set, atom_kind, xas_atom_env)

      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      INTEGER, INTENT(IN)                                :: atom_kind
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_norm_drho', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: dir, ia, ir, n, na, nr, nspins

      na = xas_atom_env%grid_atom_set(atom_kind)%grid_atom%ng_sphere
      nr = xas_atom_env%grid_atom_set(atom_kind)%grid_atom%nr
      n = na*nr
      nspins = xas_atom_env%nspins

      rho_set%norm_drhoa = 0.0_dp
      rho_set%norm_drhob = 0.0_dp
      rho_set%norm_drho = 0.0_dp

      DO dir = 1, 3
         DO ir = 1, nr
            DO ia = 1, na
               rho_set%norm_drhoa(ia, ir, 1) = rho_set%norm_drhoa(ia, ir, 1) &
                                               + rho_set%drhoa(dir)%array(ia, ir, 1)**2
            END DO !ia
         END DO !ir
      END DO !dir
      rho_set%norm_drhoa = SQRT(rho_set%norm_drhoa)

      IF (nspins == 1) THEN
         !spin-restricted
         CALL dcopy(n, rho_set%norm_drhoa(:, :, 1), 1, rho_set%norm_drhob(:, :, 1), 1)
      ELSE
         DO dir = 1, 3
            DO ir = 1, nr
               DO ia = 1, na
                  rho_set%norm_drhob(ia, ir, 1) = rho_set%norm_drhob(ia, ir, 1) &
                                                  + rho_set%drhob(dir)%array(ia, ir, 1)**2
               END DO
            END DO
         END DO
         rho_set%norm_drhob = SQRT(rho_set%norm_drhob)
      END IF

      DO dir = 1, 3
         DO ir = 1, nr
            DO ia = 1, na
               rho_set%norm_drho(ia, ir, 1) = rho_set%norm_drho(ia, ir, 1) + &
                                              (rho_set%drhoa(dir)%array(ia, ir, 1) + &
                                               rho_set%drhob(dir)%array(ia, ir, 1))**2
            END DO
         END DO
      END DO
      rho_set%norm_drho = SQRT(rho_set%norm_drho)

   END SUBROUTINE compute_norm_drho

! **************************************************************************************************
!> \brief Precomputes the spherical orbitals of the RI basis on the excited atom grids
!> \param do_gga whether the gradient needs to be computed for GGA or not
!> \param batch_info the parallelization info to complete with so distribution info
!> \param xas_atom_env ...
!> \param qs_env ...
!> \note the functions are split in a purely angular part of size na and a purely radial part of
!>       size nr. The full function on the grid can simply be obtained with dgemm and we save space
! **************************************************************************************************
   SUBROUTINE precompute_so_dso(do_gga, batch_info, xas_atom_env, qs_env)

      LOGICAL, INTENT(IN)                                :: do_gga
      TYPE(batch_info_type)                              :: batch_info
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'precompute_so_dso', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: bo(2), dir, handle, ikind, ipgf, iset, &
                                                            iso, iso_proc, l, maxso, n, na, nkind, &
                                                            nr, nset, nso_proc, nsotot, starti
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf, nsgf_set
      INTEGER, DIMENSION(:, :), POINTER                  :: so_proc_info
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: rexp
      REAL(dp), DIMENSION(:, :), POINTER                 :: dgr1, dgr2, ga, gr, slm, zet
      REAL(dp), DIMENSION(:, :, :), POINTER              :: dga1, dga2, dslm_dxyz
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, harmonics, grid_atom, slm, dslm_dxyz, ri_basis, lmax, lmin, npgf)
      NULLIFY (nsgf_set, zet, para_env, so_proc_info)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, para_env=para_env)
      nkind = SIZE(qs_kind_set)

      ALLOCATE (batch_info%so_proc_info(nkind))
      ALLOCATE (batch_info%nso_proc(nkind))
      ALLOCATE (batch_info%so_bo(2, nkind))

      DO ikind = 1, nkind

         NULLIFY (xas_atom_env%ga(ikind)%array)
         NULLIFY (xas_atom_env%gr(ikind)%array)
         NULLIFY (xas_atom_env%dga1(ikind)%array)
         NULLIFY (xas_atom_env%dga2(ikind)%array)
         NULLIFY (xas_atom_env%dgr1(ikind)%array)
         NULLIFY (xas_atom_env%dgr2(ikind)%array)

         NULLIFY (batch_info%so_proc_info(ikind)%array)

         IF (.NOT. ANY(xas_atom_env%excited_kinds == ikind)) CYCLE

         !grid info
         harmonics => xas_atom_env%harmonics_atom_set(ikind)%harmonics_atom
         grid_atom => xas_atom_env%grid_atom_set(ikind)%grid_atom

         na = grid_atom%ng_sphere
         nr = grid_atom%nr
         n = na*nr

         slm => harmonics%slm
         dslm_dxyz => harmonics%dslm_dxyz

         !basis info
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=ri_basis, basis_type="RI_XAS")
         CALL get_gto_basis_set(ri_basis, lmax=lmax, npgf=npgf, zet=zet, nset=nset, &
                                nsgf_set=nsgf_set, lmin=lmin, maxso=maxso)
         nsotot = maxso*nset

         !we split all so among the processors of the batch
         bo = get_limit(nsotot, batch_info%batch_size, batch_info%ipe)
         nso_proc = bo(2) - bo(1) + 1
         batch_info%so_bo(:, ikind) = bo
         batch_info%nso_proc(ikind) = nso_proc

         !store info about the so's set, pgf and index
         ALLOCATE (batch_info%so_proc_info(ikind)%array(3, nso_proc))
         so_proc_info => batch_info%so_proc_info(ikind)%array
         so_proc_info = -1 !default is -1 => set so value to zero
         DO iset = 1, nset
            DO ipgf = 1, npgf(iset)
               starti = (iset - 1)*maxso + (ipgf - 1)*nsoset(lmax(iset))
               DO iso = nsoset(lmin(iset) - 1) + 1, nsoset(lmax(iset))

                  !only consider so that are on this proc
                  IF (starti + iso < bo(1) .OR. starti + iso > bo(2)) CYCLE
                  iso_proc = starti + iso - bo(1) + 1
                  so_proc_info(1, iso_proc) = iset
                  so_proc_info(2, iso_proc) = ipgf
                  so_proc_info(3, iso_proc) = iso

               END DO
            END DO
         END DO

         !Put the gaussians and their gradient as purely angular or radial arrays
         ALLOCATE (xas_atom_env%ga(ikind)%array(na, nso_proc))
         ALLOCATE (xas_atom_env%gr(ikind)%array(nr, nso_proc))
         xas_atom_env%ga(ikind)%array = 0.0_dp; xas_atom_env%gr(ikind)%array = 0.0_dp
         IF (do_gga) THEN
            ALLOCATE (xas_atom_env%dga1(ikind)%array(na, nso_proc, 3))
            ALLOCATE (xas_atom_env%dgr1(ikind)%array(nr, nso_proc))
            ALLOCATE (xas_atom_env%dga2(ikind)%array(na, nso_proc, 3))
            ALLOCATE (xas_atom_env%dgr2(ikind)%array(nr, nso_proc))
            xas_atom_env%dga1(ikind)%array = 0.0_dp; xas_atom_env%dgr1(ikind)%array = 0.0_dp
            xas_atom_env%dga2(ikind)%array = 0.0_dp; xas_atom_env%dgr2(ikind)%array = 0.0_dp
         END IF

         ga => xas_atom_env%ga(ikind)%array
         gr => xas_atom_env%gr(ikind)%array
         dga1 => xas_atom_env%dga1(ikind)%array
         dga2 => xas_atom_env%dga2(ikind)%array
         dgr1 => xas_atom_env%dgr1(ikind)%array
         dgr2 => xas_atom_env%dgr2(ikind)%array

         ALLOCATE (rexp(nr))

         DO iso_proc = 1, nso_proc
            iset = so_proc_info(1, iso_proc)
            ipgf = so_proc_info(2, iso_proc)
            iso = so_proc_info(3, iso_proc)
            IF (iso < 0) CYCLE

            l = indso(1, iso)

            !radial part of the gaussian
            rexp(1:nr) = EXP(-zet(ipgf, iset)*grid_atom%rad2(1:nr))
            gr(1:nr, iso_proc) = grid_atom%rad(1:nr)**l*rexp(1:nr)

            !angular part of the gaussian
            ga(1:na, iso_proc) = slm(1:na, iso)

            IF (do_gga) THEN
               !radial part of the gradient => same in all three direction
               dgr1(1:nr, iso_proc) = grid_atom%rad(1:nr)**(l - 1)*rexp(1:nr)
               dgr2(1:nr, iso_proc) = -2.0_dp*zet(ipgf, iset)*grid_atom%rad(1:nr)**(l + 1)*rexp(1:nr)

               !angular part of the gradient
               DO dir = 1, 3
                  dga1(1:na, iso_proc, dir) = dslm_dxyz(dir, 1:na, iso)
                  dga2(1:na, iso_proc, dir) = harmonics%a(dir, 1:na)*slm(1:na, iso)
               END DO
            END IF

         END DO !iso_proc

         DEALLOCATE (rexp)
      END DO !ikind

      CALL timestop(handle)

   END SUBROUTINE precompute_so_dso

! **************************************************************************************************
!> \brief Integrate the xc kernel as a function of r on the atomic grids for the RI_XAS basis
!> \param int_fxc the global array containing the (P|fxc|Q) integrals, for all spin configurations
!> \param xas_atom_env ...
!> \param xas_tdp_control ...
!> \param qs_env ...
!> \note Note that if closed-shell, alpha-alpha term and beta-beta terms are the same
!>       Store the (P|fxc|Q) integrals on the processor they were computed on
!>       int_fxc(1)%matrix is alpha-alpha, 2: alpha-beta, 3: beta-beta
! **************************************************************************************************
   SUBROUTINE integrate_fxc_atoms(int_fxc, xas_atom_env, xas_tdp_control, qs_env)

      TYPE(cp_2d_r_p_type), DIMENSION(:, :), POINTER     :: int_fxc
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_fxc_atoms', &
         routineP = moduleN//':'//routineN

      INTEGER :: batch_size, dir, er, ex_bo(2), handle, i, iatom, ibatch, iex, ikind, inb, ipe, &
         mepos, na, natom, nb, nb_bo(2), nbatch, nbk, nex_atom, nr, num_pe, sr, subcomm
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(:), POINTER                     :: exat_neighbors
      LOGICAL                                            :: do_gga, do_sc, do_sf
      TYPE(batch_info_type)                              :: batch_info
      TYPE(cp_1d_r_p_type), DIMENSION(:, :), POINTER     :: ri_dcoeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: input, xc_functionals
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: rho_set

      NULLIFY (particle_set, qs_kind_set, dft_control, deriv_set, para_env, exat_neighbors)
      NULLIFY (rho_set, input, xc_functionals)

      CALL timeset(routineN, handle)

!  Initialize
      CALL get_qs_env(qs_env, particle_set=particle_set, qs_kind_set=qs_kind_set, natom=natom, &
                      dft_control=dft_control, input=input, para_env=para_env)
      ALLOCATE (int_fxc(natom, 4))
      DO iatom = 1, natom
         DO i = 1, 4
            NULLIFY (int_fxc(iatom, i)%array)
         END DO
      END DO
      nex_atom = SIZE(xas_atom_env%excited_atoms)
      !spin conserving in the general sense here
      do_sc = xas_tdp_control%do_spin_cons .OR. xas_tdp_control%do_singlet .OR. xas_tdp_control%do_triplet
      do_sf = xas_tdp_control%do_spin_flip

!  Get some info on the functionals
      xc_functionals => section_vals_get_subs_vals(input, "DFT%XAS_TDP%KERNEL%XC_FUNCTIONAL")
      ! ask for lsd in any case
      needs = xc_functionals_get_needs(xc_functionals, lsd=.TRUE., add_basic_components=.TRUE.)
      do_gga = needs%drho_spin !because either LDA or GGA, and the former does not need gradient

!  Distribute the excited atoms over batches of processors
!  Then, the spherical orbital of the RI basis are distributed among the procs of the batch, making
!  the GGA integration very efficient
      num_pe = para_env%num_pe
      mepos = para_env%mepos

      !create a batch_info_type
      CALL get_proc_batch_sizes(batch_size, nbatch, nex_atom, num_pe)

      !the batch index
      ibatch = mepos/batch_size
      !the proc index within the batch
      ipe = MODULO(mepos, batch_size)

      batch_info%batch_size = batch_size
      batch_info%nbatch = nbatch
      batch_info%ibatch = ibatch
      batch_info%ipe = ipe

      !create a subcommunicator for this batch
      CALL mp_comm_split_direct(para_env%group, subcomm, ibatch)
      NULLIFY (batch_info%para_env)
      CALL cp_para_env_create(batch_info%para_env, subcomm)

!  Precompute the spherical orbital of the RI basis (and maybe their gradient) on the grids of the
!  excited atoms. Needed for the GGA integration and to actually put the density on the grid
      CALL precompute_so_dso(do_gga, batch_info, xas_atom_env, qs_env)

      !distribute the excted atoms over the batches
      ex_bo = get_limit(nex_atom, nbatch, ibatch)

!  Looping over the excited atoms
      DO iex = ex_bo(1), ex_bo(2)

         iatom = xas_atom_env%excited_atoms(iex)
         ikind = particle_set(iatom)%atomic_kind%kind_number
         exat_neighbors => xas_atom_env%exat_neighbors(iex)%array
         ri_dcoeff => xas_atom_env%ri_dcoeff(:, :, iex)

!     General grid/basis info
         na = xas_atom_env%grid_atom_set(ikind)%grid_atom%ng_sphere
         nr = xas_atom_env%grid_atom_set(ikind)%grid_atom%nr

!     Creating a xc_rho_set to store the density and dset for the kernel
         bounds(1:2, 1:3) = 1
         bounds(2, 1) = na
         bounds(2, 2) = nr

         CALL xc_rho_set_create(rho_set=rho_set, local_bounds=bounds, &
                                rho_cutoff=dft_control%qs_control%eps_rho_rspace, &
                                drho_cutoff=dft_control%qs_control%eps_rho_rspace)
         CALL xc_dset_create(deriv_set, local_bounds=bounds)

         ! allocate internals of the rho_set
         CALL xc_rho_set_atom_update(rho_set, needs, nspins=2, bo=bounds)

!     Put the density, and possibly its gradient,  on the grid (for this atom)
         CALL put_density_on_atomic_grid(rho_set, ri_dcoeff(iatom, :), ikind, &
                                         do_gga, batch_info, xas_atom_env, qs_env)

!     Take the neighboring atom contributions to the density (and gradient)
!     distribute the grid among the procs (for best load balance)
         nb_bo = get_limit(nr, batch_size, ipe)
         sr = nb_bo(1); er = nb_bo(2)
         DO inb = 1, SIZE(exat_neighbors)

            nb = exat_neighbors(inb)
            nbk = particle_set(nb)%atomic_kind%kind_number
            CALL put_density_on_other_grid(rho_set, ri_dcoeff(nb, :), nb, nbk, iatom, &
                                           ikind, sr, er, do_gga, xas_atom_env, qs_env)

         END DO

         ! make sure contributions from different procs are summed up
         CALL mp_sum(rho_set%rhoa, batch_info%para_env%group)
         CALL mp_sum(rho_set%rhob, batch_info%para_env%group)
         IF (do_gga) THEN
            DO dir = 1, 3
               CALL mp_sum(rho_set%drhoa(dir)%array, batch_info%para_env%group)
               CALL mp_sum(rho_set%drhob(dir)%array, batch_info%para_env%group)
            END DO
         END IF

!     In case of GGA, also need the norm of the density gradient
         IF (do_gga) CALL compute_norm_drho(rho_set, ikind, xas_atom_env)

!     Compute the required derivatives
         CALL xc_functionals_eval(xc_functionals, lsd=.TRUE., rho_set=rho_set, deriv_set=deriv_set, &
                                  deriv_order=2)

         !spin-conserving (LDA part)
         IF (do_sc) THEN
            CALL integrate_sc_fxc(int_fxc, iatom, ikind, deriv_set, xas_atom_env, qs_env)
         END IF

         !spin-flip (LDA part)
         IF (do_sf) THEN
            CALL integrate_sf_fxc(int_fxc, iatom, ikind, rho_set, deriv_set, xas_atom_env, qs_env)
         END IF

         !Gradient correction (note: spin-flip only keeps the lda part, aka ALDA0)
         IF (do_gga .AND. do_sc) THEN
            CALL integrate_gga_fxc(int_fxc, iatom, ikind, batch_info, rho_set, deriv_set, &
                                   xas_atom_env, qs_env)
         END IF

!     Clean-up
         CALL xc_dset_release(deriv_set)
         CALL xc_rho_set_release(rho_set)
      END DO !iex

      CALL release_batch_info(batch_info)

      !Not necessary to sync, but makes sure that any load inbalance is reported here
      CALL mp_sync(para_env%group)

      CALL timestop(handle)

   END SUBROUTINE integrate_fxc_atoms

! **************************************************************************************************
!> \brief Integrate the gradient correction part of the xc kernel on the atomic grid
!> \param int_fxc the array containing the (P|fxc|Q) integrals
!> \param iatom the index of the current excited atom
!> \param ikind the index of the current excited kind
!> \param batch_info how the so are distributed over the processor batch
!> \param rho_set the variable contatinind the density and its gradient
!> \param deriv_set the functional derivatives
!> \param xas_atom_env ...
!> \param qs_env ...
!> \note Ignored in case of pure LDA, added on top of the LDA kernel in case of GGA
! **************************************************************************************************
   SUBROUTINE integrate_gga_fxc(int_fxc, iatom, ikind, batch_info, rho_set, deriv_set, &
                                xas_atom_env, qs_env)

      TYPE(cp_2d_r_p_type), DIMENSION(:, :), POINTER     :: int_fxc
      INTEGER, INTENT(IN)                                :: iatom, ikind
      TYPE(batch_info_type)                              :: batch_info
      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_gga_fxc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: bo(2), dir, handle, i, ia, ir, jpgf, &
                                                            jset, jso, l, maxso, na, nr, nset, &
                                                            nsgf, nsoi, nsotot, startj, ub
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: rexp
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: int_sgf, res, so, work
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: dso
      REAL(dp), DIMENSION(:, :), POINTER                 :: dgr1, dgr2, ga, gr, ri_sphi_so, weight, &
                                                            zet
      REAL(dp), DIMENSION(:, :, :), POINTER              :: dga1, dga2
      TYPE(cp_2d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: int_so, vxc
      TYPE(cp_3d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: vxg
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (grid_atom, ri_basis, qs_kind_set, ga, gr, dgr1, dgr2, lmax, lmin, npgf)
      NULLIFY (weight, ri_sphi_so, dga1, dga2, para_env, harmonics, zet)

      !Strategy: we need to compute <phi_i|fxc|phij>, most of existing application of the 2nd
      !          functional derivative involve the response density, and the expression of the
      !          integral int (fxc*n^1) is well known. We substitute the spherical orbital phi_j
      !          in place of n^1 in the formula and thus perform the first integration. Then
      !          we obtain something in the form int (fxc*phi_j) = Vxc - div (Vxg) that we can
      !          put on the grid and treat like a potential. The second integration is done by
      !          using the divergence theorem and numerical integration:
      !          <phi_i|fxc|phi_j> = int phi_i*(Vxc - div(Vxg)) = int phi_i*Vxc + grad(phi_i).Vxg
      !          Note the sign change and the dot product.

      CALL timeset(routineN, handle)

      !If closed shell, only compute f_aa and f_ab (ub = 2)
      ub = 2
      IF (xas_atom_env%nspins == 2) ub = 3

      !Get the necessary grid info
      harmonics => xas_atom_env%harmonics_atom_set(ikind)%harmonics_atom
      grid_atom => xas_atom_env%grid_atom_set(ikind)%grid_atom
      na = grid_atom%ng_sphere
      nr = grid_atom%nr
      weight => grid_atom%weight

      !get the ri_basis info
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, para_env=para_env)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=ri_basis, basis_type="RI_XAS")

      CALL get_gto_basis_set(gto_basis_set=ri_basis, lmax=lmax, lmin=lmin, nsgf=nsgf, &
                             maxso=maxso, npgf=npgf, nset=nset, zet=zet)
      nsotot = nset*maxso

      !Point to the precomputed so
      ga => xas_atom_env%ga(ikind)%array
      gr => xas_atom_env%gr(ikind)%array
      dgr1 => xas_atom_env%dgr1(ikind)%array
      dgr2 => xas_atom_env%dgr2(ikind)%array
      dga1 => xas_atom_env%dga1(ikind)%array
      dga2 => xas_atom_env%dga2(ikind)%array

      !Before integration, wanna pre-divide all relevant derivastives by the nrom of the gradient
      CALL divide_by_norm_drho(deriv_set, rho_set, lsd=.TRUE.)

      !Wanna integrate <phi_i|fxc|phi_j>, start looping over phi_j and do the first integration, then
      !collect vxc and vxg and loop over phi_i for the second integration
      !Note: we do not use the CG coefficients because they are only useful when there is a product
      !      of Gaussians, which is not really the case here
      !Note: the spherical orbitals for phi_i are distributed among the prcos of the current batch

      ALLOCATE (so(na, nr))
      ALLOCATE (dso(na, nr, 3))
      ALLOCATE (rexp(nr))

      ALLOCATE (vxc(ub))
      ALLOCATE (vxg(ub))
      ALLOCATE (int_so(ub))
      DO i = 1, ub
         ALLOCATE (vxc(i)%array(na, nr))
         ALLOCATE (vxg(i)%array(na, nr, 3))
         ALLOCATE (int_so(i)%array(nsotot, nsotot))
         vxc(i)%array = 0.0_dp; vxg(i)%array = 0.0_dp; int_so(i)%array = 0.0_dp
      END DO

      DO jset = 1, nset
         DO jpgf = 1, npgf(jset)
            startj = (jset - 1)*maxso + (jpgf - 1)*nsoset(lmax(jset))
            DO jso = nsoset(lmin(jset) - 1) + 1, nsoset(lmax(jset))
               l = indso(1, jso)

               !put the so phi_j and its gradient on the grid
               !more efficient to recompute it rather than mp_bcast each chunk

               rexp(1:nr) = EXP(-zet(jpgf, jset)*grid_atom%rad2(1:nr))
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE), &
!$OMP SHARED(nr,na,so,dso,grid_atom,l,rexp,harmonics,jso,zet,jset,jpgf), &
!$OMP PRIVATE(ir,ia,dir)
               DO ir = 1, nr
                  DO ia = 1, na

                     !so
                     so(ia, ir) = grid_atom%rad(ir)**l*rexp(ir)*harmonics%slm(ia, jso)

                     !dso
                     dso(ia, ir, :) = 0.0_dp
                     DO dir = 1, 3
                        dso(ia, ir, dir) = dso(ia, ir, dir) &
                                           + grid_atom%rad(ir)**(l - 1)*rexp(ir)*harmonics%dslm_dxyz(dir, ia, jso) &
                                           - 2.0_dp*zet(jpgf, jset)*grid_atom%rad(ir)**(l + 1)*rexp(ir) &
                                           *harmonics%a(dir, ia)*harmonics%slm(ia, jso)
                     END DO
                  END DO
               END DO
!$OMP END PARALLEL DO

               !Perform the first integration (analytically)
               CALL get_vxc_vxg(vxc, vxg, so, dso, na, nr, rho_set, deriv_set, weight)

               !For a given phi_j, compute the second integration with all phi_i at once
               !=> allows for efficient gemm to take place, especially since so are distributed
               nsoi = batch_info%nso_proc(ikind)
               bo = batch_info%so_bo(:, ikind)
               ALLOCATE (res(nsoi, nsoi), work(na, nsoi))
               res = 0.0_dp; work = 0.0_dp

               DO i = 1, ub

                  !integrate so*Vxc and store in the int_so
                  CALL dgemm('N', 'N', na, nsoi, nr, 1.0_dp, vxc(i)%array(:, :), na, &
                             gr(:, 1:nsoi), nr, 0.0_dp, work, na)
                  CALL dgemm('T', 'N', nsoi, nsoi, na, 1.0_dp, work, na, &
                             ga(:, 1:nsoi), na, 0.0_dp, res, nsoi)
                  int_so(i)%array(bo(1):bo(2), startj + jso) = get_diag(res)

                  DO dir = 1, 3

                     ! integrate and sum up Vxg*dso
                     CALL dgemm('N', 'N', na, nsoi, nr, 1.0_dp, vxg(i)%array(:, :, dir), na, &
                                dgr1(:, 1:nsoi), nr, 0.0_dp, work, na)
                     CALL dgemm('T', 'N', nsoi, nsoi, na, 1.0_dp, work, na, &
                                dga1(:, 1:nsoi, dir), na, 0.0_dp, res, nsoi)
                     CALL daxpy(nsoi, 1.0_dp, get_diag(res), 1, int_so(i)%array(bo(1):bo(2), startj + jso), 1)

                     CALL dgemm('N', 'N', na, nsoi, nr, 1.0_dp, vxg(i)%array(:, :, dir), na, &
                                dgr2(:, 1:nsoi), nr, 0.0_dp, work, na)
                     CALL dgemm('T', 'N', nsoi, nsoi, na, 1.0_dp, work, na, &
                                dga2(:, 1:nsoi, dir), na, 0.0_dp, res, nsoi)
                     CALL daxpy(nsoi, 1.0_dp, get_diag(res), 1, int_so(i)%array(bo(1):bo(2), startj + jso), 1)

                  END DO

               END DO !i
               DEALLOCATE (res, work)

            END DO !jso
         END DO !jpgf
      END DO !jset

      !Contract into sgf and add to already computed LDA part of int_fxc
      ri_sphi_so => xas_atom_env%ri_sphi_so(ikind)%array
      ALLOCATE (int_sgf(nsgf, nsgf))
      DO i = 1, ub
         CALL mp_sum(int_so(i)%array, batch_info%para_env%group)
         CALL contract_so2sgf(int_sgf, int_so(i)%array, ri_basis, ri_sphi_so)
         CALL daxpy(nsgf*nsgf, 1.0_dp, int_sgf, 1, int_fxc(iatom, i)%array, 1)
      END DO

      !Clean-up
      DO i = 1, ub
         DEALLOCATE (vxc(i)%array)
         DEALLOCATE (vxg(i)%array)
         DEALLOCATE (int_so(i)%array)
      END DO
      DEALLOCATE (vxc, vxg, int_so)

      CALL timestop(handle)

   END SUBROUTINE integrate_gga_fxc

! **************************************************************************************************
!> \brief Computes the first integration of the GGA part of <phi_i|fxc|phi_j>, i.e. int fxc*phi_j.
!>        The result is of the form Vxc - div(Vxg). Up to 3 results are returned, correspoinding to
!>        f_aa, f_ab and (if open-shell) f_bb
!> \param vxc ...
!> \param vxg ...
!> \param so the spherical orbital on the grid
!> \param dso the derivative of the spherical orbital on the grid
!> \param na ...
!> \param nr ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param weight the grid weight
!> \note This method is extremely similar to xc_calc_2nd_deriv of xc.F, but because it is a special
!>       case that can be further optimized and because the interface of the original routine does
!>       not fit this code, it has been re-written (no pw, no rho1_set but just the so, etc...)
! **************************************************************************************************
   SUBROUTINE get_vxc_vxg(vxc, vxg, so, dso, na, nr, rho_set, deriv_set, weight)

      TYPE(cp_2d_r_p_type), DIMENSION(:)                 :: vxc
      TYPE(cp_3d_r_p_type), DIMENSION(:)                 :: vxg
      REAL(dp), DIMENSION(:, :)                          :: so
      REAL(dp), DIMENSION(:, :, :)                       :: dso
      INTEGER, INTENT(IN)                                :: na, nr
      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      REAL(dp), DIMENSION(:, :), POINTER                 :: weight

      CHARACTER(len=*), PARAMETER :: routineN = 'get_vxc_vxg', routineP = moduleN//':'//routineN

      INTEGER                                            :: dir, handle, i, ia, ir, ub
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: dot_proda, dot_prodb, tmp
      REAL(dp), DIMENSION(:, :, :), POINTER              :: d1e, d2e, norm_drhoa, norm_drhob
      TYPE(xc_derivative_type), POINTER                  :: deriv

      NULLIFY (norm_drhoa, norm_drhob, d2e, d1e, deriv)

      CALL timeset(routineN, handle)

      !Note: this routines follows the order of the terms in equaiton (A.7) of Thomas Chassaing
      !      thesis, except for the pure LDA terms that are dropped. The n^(1)_a and n^(1)_b
      !      response densities are replaced by the spherical orbital.
      !      The usual spin ordering is used: aa => 1, ab => 2 , bb => 3

      !point to the relevant components of rho_set
      ub = SIZE(vxc)
      norm_drhoa => rho_set%norm_drhoa
      norm_drhob => rho_set%norm_drhob

      !Some init
      DO i = 1, ub
         vxc(i)%array = 0.0_dp
         vxg(i)%array = 0.0_dp
      END DO

      ALLOCATE (tmp(na, nr), dot_proda(na, nr), dot_prodb(na, nr))
      dot_proda = 0.0_dp; dot_prodb = 0.0_dp

      !Strategy: most terms are either multiplied by drhoa or drhob => group those first and then
      !          multiply. Also most terms are multiplied by the dot product grad_n . grad_so, so
      !          precompute it as well

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP SHARED(dot_proda,dot_prodb,tmp,vxc,vxg,deriv_set,rho_set,na,nr,norm_drhoa,norm_drhob,dso, &
!$OMP        so,weight,ub), &
!$OMP PRIVATE(ia,ir,dir,deriv,d1e,d2e)

      !Precompute the very common dot products grad_na . grad_so and grad_nb . grad_so
      DO dir = 1, 3
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               dot_proda(ia, ir) = dot_proda(ia, ir) + rho_set%drhoa(dir)%array(ia, ir, 1)*dso(ia, ir, dir)
               dot_prodb(ia, ir) = dot_prodb(ia, ir) + rho_set%drhob(dir)%array(ia, ir, 1)*dso(ia, ir, dir)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT
      END DO !dir

      !Deal with f_aa

      !Vxc, first term
      deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(norm_drhoa)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            vxc(1)%array(ia, ir) = d2e(ia, ir, 1)*dot_proda(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxc, second term
      deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(norm_drho)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            vxc(1)%array(ia, ir) = vxc(1)%array(ia, ir) + d2e(ia, ir, 1)*dot_prodb(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxc, take the grid weight into acocunt
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            vxc(1)%array(ia, ir) = vxc(1)%array(ia, ir)*weight(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxg, first term (to be multiplied by drhoa)
      deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(norm_drhoa)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            tmp(ia, ir) = d2e(ia, ir, 1)*so(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxg, second term (to be multiplied by drhoa)
      deriv => xc_dset_get_derivative(deriv_set, "(norm_drhoa)(norm_drho)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            tmp(ia, ir) = tmp(ia, ir) + d2e(ia, ir, 1)*dot_prodb(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxg, third term (to be multiplied by drhoa)
      deriv => xc_dset_get_derivative(deriv_set, "(norm_drhoa)(norm_drhoa)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            tmp(ia, ir) = tmp(ia, ir) + d2e(ia, ir, 1)*dot_proda(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxg, fourth term (to be multiplied by drhoa)
      deriv => xc_dset_get_derivative(deriv_set, "(norm_drhoa)")
      CALL xc_derivative_get(deriv, deriv_data=d1e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            tmp(ia, ir) = tmp(ia, ir) - d1e(ia, ir, 1)*dot_proda(ia, ir) &
                          /MAX(norm_drhoa(ia, ir, 1), rho_set%drho_cutoff)**2
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !put tmp*drhoa in Vxg (so that we can reuse it for drhob terms)
      DO dir = 1, 3
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               vxg(1)%array(ia, ir, dir) = tmp(ia, ir)*rho_set%drhoa(dir)%array(ia, ir, 1)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT
      END DO !dir

      !Vxg, fifth term (to be multiplied by drhob)
      deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(norm_drho)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            tmp(ia, ir) = d2e(ia, ir, 1)*so(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxg, sixth term (to be multiplied by drhob)
      deriv => xc_dset_get_derivative(deriv_set, "(norm_drhoa)(norm_drho)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            tmp(ia, ir) = tmp(ia, ir) + d2e(ia, ir, 1)*dot_proda(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxg, seventh term (to be multiplied by drhob)
      deriv => xc_dset_get_derivative(deriv_set, "(norm_drho)(norm_drho)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            tmp(ia, ir) = tmp(ia, ir) + d2e(ia, ir, 1)*dot_prodb(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !put tmp*drhob in Vxg
      DO dir = 1, 3
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               vxg(1)%array(ia, ir, dir) = vxg(1)%array(ia, ir, dir) + tmp(ia, ir)*rho_set%drhob(dir)%array(ia, ir, 1)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT
      END DO !dir

      !Vxg, last term
      deriv => xc_dset_get_derivative(deriv_set, "(norm_drhoa)")
      CALL xc_derivative_get(deriv, deriv_data=d1e)
      DO dir = 1, 3
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               vxg(1)%array(ia, ir, dir) = vxg(1)%array(ia, ir, dir) + d1e(ia, ir, 1)*dso(ia, ir, dir)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT
      END DO !dir

      !Vxg, take the grid weight into account
      DO dir = 1, 3
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               vxg(1)%array(ia, ir, dir) = vxg(1)%array(ia, ir, dir)*weight(ia, ir)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT
      END DO !dir

      !Deal with fab

      !Vxc, first term
      deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(norm_drhob)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            vxc(2)%array(ia, ir) = d2e(ia, ir, 1)*dot_prodb(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxc, second term
      deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(norm_drho)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            vxc(2)%array(ia, ir) = vxc(2)%array(ia, ir) + d2e(ia, ir, 1)*dot_proda(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxc, take the grid weight into acocunt
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            vxc(2)%array(ia, ir) = vxc(2)%array(ia, ir)*weight(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxg, first term (to be multiplied by drhoa)
      deriv => xc_dset_get_derivative(deriv_set, "(rhob)(norm_drhoa)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            tmp(ia, ir) = d2e(ia, ir, 1)*so(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxg, second term (to be multiplied by drhoa)
      deriv => xc_dset_get_derivative(deriv_set, "(norm_drhoa)(norm_drho)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            tmp(ia, ir) = tmp(ia, ir) + d2e(ia, ir, 1)*dot_proda(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxg, third term (to be multiplied by drhoa)
      deriv => xc_dset_get_derivative(deriv_set, "(norm_drhoa)(norm_drhob)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            tmp(ia, ir) = tmp(ia, ir) + d2e(ia, ir, 1)*dot_prodb(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !put tmp*drhoa in Vxg
      DO dir = 1, 3
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               vxg(2)%array(ia, ir, dir) = tmp(ia, ir)*rho_set%drhoa(dir)%array(ia, ir, 1)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT
      END DO !dir

      !Vxg, fourth term (to be multiplied by drhob)
      deriv => xc_dset_get_derivative(deriv_set, "(rhob)(norm_drho)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            tmp(ia, ir) = d2e(ia, ir, 1)*so(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxg, fifth term (to be multiplied by drhob)
      deriv => xc_dset_get_derivative(deriv_set, "(norm_drho)(norm_drhob)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            tmp(ia, ir) = tmp(ia, ir) + d2e(ia, ir, 1)*dot_prodb(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !Vxg, sixth term (to be multiplied by drhob)
      deriv => xc_dset_get_derivative(deriv_set, "(norm_drho)(norm_drho)")
      CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
      DO ir = 1, nr
         DO ia = 1, na
            tmp(ia, ir) = tmp(ia, ir) + d2e(ia, ir, 1)*dot_proda(ia, ir)
         END DO !ia
      END DO !ir
!$OMP END DO NOWAIT

      !put tmp*drhob in Vxg
      DO dir = 1, 3
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               vxg(2)%array(ia, ir, dir) = vxg(2)%array(ia, ir, dir) + tmp(ia, ir)*rho_set%drhob(dir)%array(ia, ir, 1)
            END DO
         END DO
!$OMP END DO NOWAIT
      END DO

      !Vxg, last term
      deriv => xc_dset_get_derivative(deriv_set, "(norm_drho)")
      CALL xc_derivative_get(deriv, deriv_data=d1e)
      DO dir = 1, 3
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               vxg(2)%array(ia, ir, dir) = vxg(2)%array(ia, ir, dir) + d1e(ia, ir, 1)*dso(ia, ir, dir)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT
      END DO !dir

      !Vxg, take the grid weight into account
      DO dir = 1, 3
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               vxg(2)%array(ia, ir, dir) = vxg(2)%array(ia, ir, dir)*weight(ia, ir)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT
      END DO !dir

      !Deal with f_bb, if so required
      IF (ub == 3) THEN

         !Vxc, first term
         deriv => xc_dset_get_derivative(deriv_set, "(rhob)(norm_drhob)")
         CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               vxc(3)%array(ia, ir) = d2e(ia, ir, 1)*dot_prodb(ia, ir)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT

         !Vxc, second term
         deriv => xc_dset_get_derivative(deriv_set, "(rhob)(norm_drho)")
         CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               vxc(3)%array(ia, ir) = vxc(3)%array(ia, ir) + d2e(ia, ir, 1)*dot_proda(ia, ir)
            END DO !i
         END DO !ir
!$OMP END DO NOWAIT

         !Vxc, take the grid weight into acocunt
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               vxc(3)%array(ia, ir) = vxc(3)%array(ia, ir)*weight(ia, ir)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT

         !Vxg, first term (to be multiplied by drhob)
         deriv => xc_dset_get_derivative(deriv_set, "(rhob)(norm_drhob)")
         CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               tmp(ia, ir) = d2e(ia, ir, 1)*so(ia, ir)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT

         !Vxg, second term (to be multiplied by drhob)
         deriv => xc_dset_get_derivative(deriv_set, "(norm_drhob)(norm_drho)")
         CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               tmp(ia, ir) = tmp(ia, ir) + d2e(ia, ir, 1)*dot_proda(ia, ir)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT

         !Vxg, third term (to be multiplied by drhob)
         deriv => xc_dset_get_derivative(deriv_set, "(norm_drhob)(norm_drhob)")
         CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               tmp(ia, ir) = tmp(ia, ir) + d2e(ia, ir, 1)*dot_prodb(ia, ir)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT

         !Vxg, fourth term (to be multiplied by drhob)
         deriv => xc_dset_get_derivative(deriv_set, "(norm_drhob)")
         CALL xc_derivative_get(deriv, deriv_data=d1e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               tmp(ia, ir) = tmp(ia, ir) - d1e(ia, ir, 1)*dot_prodb(ia, ir) &
                             /MAX(norm_drhob(ia, ir, 1), rho_set%drho_cutoff)**2
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT

         !put tmp*drhob in Vxg (so that we can reuse it for drhoa terms)
         DO dir = 1, 3
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
            DO ir = 1, nr
               DO ia = 1, na
                  vxg(3)%array(ia, ir, dir) = tmp(ia, ir)*rho_set%drhob(dir)%array(ia, ir, 1)
               END DO !ia
            END DO !ir
!$OMP END DO NOWAIT
         END DO !dir

         !Vxg, fifth term (to be multiplied by drhoa)
         deriv => xc_dset_get_derivative(deriv_set, "(rhob)(norm_drho)")
         CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               tmp(ia, ir) = d2e(ia, ir, 1)*so(ia, ir)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT

         !Vxg, sixth term (to be multiplied by drhoa)
         deriv => xc_dset_get_derivative(deriv_set, "(norm_drhob)(norm_drho)")
         CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               tmp(ia, ir) = tmp(ia, ir) + d2e(ia, ir, 1)*dot_prodb(ia, ir)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT

         !Vxg, seventh term (to be multiplied by drhoa)
         deriv => xc_dset_get_derivative(deriv_set, "(norm_drho)(norm_drho)")
         CALL xc_derivative_get(deriv, deriv_data=d2e)
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
         DO ir = 1, nr
            DO ia = 1, na
               tmp(ia, ir) = tmp(ia, ir) + d2e(ia, ir, 1)*dot_proda(ia, ir)
            END DO !ia
         END DO !ir
!$OMP END DO NOWAIT

         !put tmp*drhoa in Vxg
         DO dir = 1, 3
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
            DO ir = 1, nr
               DO ia = 1, na
                  vxg(3)%array(ia, ir, dir) = vxg(3)%array(ia, ir, dir) + &
                                              tmp(ia, ir)*rho_set%drhoa(dir)%array(ia, ir, 1)
               END DO !ia
            END DO !ir
!$OMP END DO NOWAIT
         END DO !dir

         !Vxg, last term
         deriv => xc_dset_get_derivative(deriv_set, "(norm_drhob)")
         CALL xc_derivative_get(deriv, deriv_data=d1e)
         DO dir = 1, 3
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
            DO ir = 1, nr
               DO ia = 1, na
                  vxg(3)%array(ia, ir, dir) = vxg(3)%array(ia, ir, dir) + d1e(ia, ir, 1)*dso(ia, ir, dir)
               END DO !ia
            END DO !ir
!$OMP END DO NOWAIT
         END DO !dir

         !Vxg, take the grid weight into account
         DO dir = 1, 3
!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
            DO ir = 1, nr
               DO ia = 1, na
                  vxg(3)%array(ia, ir, dir) = vxg(3)%array(ia, ir, dir)*weight(ia, ir)
               END DO !ia
            END DO !ir
!$OMP END DO NOWAIT
         END DO !dir

      END IF !f_bb

!$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE get_vxc_vxg

! **************************************************************************************************
!> \brief Integrate the fxc kernel in the spin-conserving case, be it closed- or open-shell
!> \param int_fxc the array containing the (P|fxc|Q) integrals
!> \param iatom the index of the current excited atom
!> \param ikind the index of the current excited kind
!> \param deriv_set the set of functional derivatives
!> \param xas_atom_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE integrate_sc_fxc(int_fxc, iatom, ikind, deriv_set, xas_atom_env, qs_env)

      TYPE(cp_2d_r_p_type), DIMENSION(:, :), POINTER     :: int_fxc
      INTEGER, INTENT(IN)                                :: iatom, ikind
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_sc_fxc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, maxso, na, nr, nset, nsotot, nspins, &
                                                            ri_nsgf, ub
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: fxc, int_so
      REAL(dp), DIMENSION(:, :), POINTER                 :: ri_sphi_so
      TYPE(cp_3d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: d2e
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xc_derivative_type), POINTER                  :: deriv

      NULLIFY (grid_atom, deriv, ri_basis, ri_sphi_so, qs_kind_set, dft_control)

      ! Initialization
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, dft_control=dft_control)
      grid_atom => xas_atom_env%grid_atom_set(ikind)%grid_atom
      na = grid_atom%ng_sphere
      nr = grid_atom%nr
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=ri_basis, basis_type="RI_XAS")
      CALL get_gto_basis_set(ri_basis, nset=nset, maxso=maxso, nsgf=ri_nsgf)
      nsotot = nset*maxso
      ri_sphi_so => xas_atom_env%ri_sphi_so(ikind)%array
      nspins = dft_control%nspins

      ! Get the second derivatives
      ALLOCATE (d2e(3))
      deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(rhoa)")
      CALL xc_derivative_get(deriv, deriv_data=d2e(1)%array)
      deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(rhob)")
      CALL xc_derivative_get(deriv, deriv_data=d2e(2)%array)
      deriv => xc_dset_get_derivative(deriv_set, "(rhob)(rhob)")
      CALL xc_derivative_get(deriv, deriv_data=d2e(3)%array)

      ! Allocate some work arrays
      ALLOCATE (fxc(na, nr))
      ALLOCATE (int_so(nsotot, nsotot))

      ! Integrate for all three derivatives, taking the grid weight into account
      ! If closed shell, do not need to integrate beta-beta as it is the same as alpha-alpha
      ub = 2; IF (nspins == 2) ub = 3
      DO i = 1, ub

         int_so = 0.0_dp
         fxc(:, :) = d2e(i)%array(:, :, 1)*grid_atom%weight(:, :)
         CALL integrate_so_prod(int_so, fxc, ikind, xas_atom_env, qs_env)

         !contract into sgf. Array allocated on current processor only
         ALLOCATE (int_fxc(iatom, i)%array(ri_nsgf, ri_nsgf))
         int_fxc(iatom, i)%array = 0.0_dp
         CALL contract_so2sgf(int_fxc(iatom, i)%array, int_so, ri_basis, ri_sphi_so)

      END DO

   END SUBROUTINE integrate_sc_fxc

! **************************************************************************************************
!> \brief Integrate the fxc kernel in the spin-flip case (open-shell assumed). The spin-flip LDA
!>        kernel reads: fxc = 1/(rhoa - rhob) * (dE/drhoa - dE/drhob)
!> \param int_fxc the array containing the (P|fxc|Q) integrals
!> \param iatom the index of the current excited atom
!> \param ikind the index of the current excited kind
!> \param rho_set the density in the atomic grid
!> \param deriv_set the set of functional derivatives
!> \param xas_atom_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE integrate_sf_fxc(int_fxc, iatom, ikind, rho_set, deriv_set, xas_atom_env, qs_env)

      TYPE(cp_2d_r_p_type), DIMENSION(:, :), POINTER     :: int_fxc
      INTEGER, INTENT(IN)                                :: iatom, ikind
      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_sf_fxc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ia, ir, maxso, na, nr, nset, nsotot, &
                                                            ri_nsgf
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: fxc, int_so
      REAL(dp), DIMENSION(:, :), POINTER                 :: ri_sphi_so
      REAL(dp), DIMENSION(:, :, :), POINTER              :: rhoa, rhob
      TYPE(cp_3d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: d1e, d2e
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xc_derivative_type), POINTER                  :: deriv

      NULLIFY (grid_atom, deriv, ri_basis, ri_sphi_so, qs_kind_set, rhoa, rhob, dft_control)

      ! Initialization
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, dft_control=dft_control)
      grid_atom => xas_atom_env%grid_atom_set(ikind)%grid_atom
      na = grid_atom%ng_sphere
      nr = grid_atom%nr
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=ri_basis, basis_type="RI_XAS")
      CALL get_gto_basis_set(ri_basis, nset=nset, maxso=maxso, nsgf=ri_nsgf)
      nsotot = nset*maxso
      ri_sphi_so => xas_atom_env%ri_sphi_so(ikind)%array
      rhoa => rho_set%rhoa
      rhob => rho_set%rhob

      ALLOCATE (d1e(2))
      deriv => xc_dset_get_derivative(deriv_set, "(rhoa)")
      CALL xc_derivative_get(deriv, deriv_data=d1e(1)%array)
      deriv => xc_dset_get_derivative(deriv_set, "(rhob)")
      CALL xc_derivative_get(deriv, deriv_data=d1e(2)%array)

      ! In case rhoa -> rhob, take the limit, which involves the (already computed) second derivatives
      ALLOCATE (d2e(3))
      deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(rhoa)")
      CALL xc_derivative_get(deriv, deriv_data=d2e(1)%array)
      deriv => xc_dset_get_derivative(deriv_set, "(rhoa)(rhob)")
      CALL xc_derivative_get(deriv, deriv_data=d2e(2)%array)
      deriv => xc_dset_get_derivative(deriv_set, "(rhob)(rhob)")
      CALL xc_derivative_get(deriv, deriv_data=d2e(3)%array)

      !Compute the kernel on the grid. Already take weight into acocunt there
      ALLOCATE (fxc(na, nr))
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(NONE), &
!$OMP SHARED(grid_atom,fxc,d1e,d2e,dft_control,na,nr,rhoa,rhob), &
!$OMP PRIVATE(ia,ir)
      DO ir = 1, nr
         DO ia = 1, na

            !Need to be careful not to divide by zero. Assume that if rhoa == rhob, then
            !take the limit fxc = 0.5* (f_aa + f_bb - 2f_ab)
            IF (ABS(rhoa(ia, ir, 1) - rhob(ia, ir, 1)) > dft_control%qs_control%eps_rho_rspace) THEN
               fxc(ia, ir) = grid_atom%weight(ia, ir)/(rhoa(ia, ir, 1) - rhob(ia, ir, 1)) &
                             *(d1e(1)%array(ia, ir, 1) - d1e(2)%array(ia, ir, 1))
            ELSE
               fxc(ia, ir) = 0.5_dp*grid_atom%weight(ia, ir)* &
                             (d2e(1)%array(ia, ir, 1) + d2e(3)%array(ia, ir, 1) - 2._dp*d2e(2)%array(ia, ir, 1))
            END IF

         END DO
      END DO
!$OMP END PARALLEL DO

      ! Integrate wrt to so
      ALLOCATE (int_so(nsotot, nsotot))
      int_so = 0.0_dp
      CALL integrate_so_prod(int_so, fxc, ikind, xas_atom_env, qs_env)

      ! Contract into sgf. Array located on current processor only
      ALLOCATE (int_fxc(iatom, 4)%array(ri_nsgf, ri_nsgf))
      int_fxc(iatom, 4)%array = 0.0_dp
      CALL contract_so2sgf(int_fxc(iatom, 4)%array, int_so, ri_basis, ri_sphi_so)

   END SUBROUTINE integrate_sf_fxc

! **************************************************************************************************
!> \brief Contract spherical orbitals to spherical Gaussians (so to sgf)
!> \param int_sgf the array with the sgf integrals
!> \param int_so the array with the so integrals (to contract)
!> \param basis the corresponding gto basis set
!> \param sphi_so the contraction coefficients for the s:
! **************************************************************************************************
   SUBROUTINE contract_so2sgf(int_sgf, int_so, basis, sphi_so)

      REAL(dp), DIMENSION(:, :)                          :: int_sgf, int_so
      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), DIMENSION(:, :)                          :: sphi_so

      CHARACTER(len=*), PARAMETER :: routineN = 'contract_so2sgf', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iset, jset, maxso, nset, nsoi, nsoj, &
                                                            sgfi, sgfj, starti, startj
      INTEGER, DIMENSION(:), POINTER                     :: lmax, npgf, nsgf_set
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf

      NULLIFY (nsgf_set, npgf, lmax, first_sgf)

      CALL get_gto_basis_set(basis, nset=nset, maxso=maxso, nsgf_set=nsgf_set, first_sgf=first_sgf, &
                             npgf=npgf, lmax=lmax)

      DO iset = 1, nset
         starti = (iset - 1)*maxso + 1
         nsoi = npgf(iset)*nsoset(lmax(iset))
         sgfi = first_sgf(1, iset)

         DO jset = 1, nset
            startj = (jset - 1)*maxso + 1
            nsoj = npgf(jset)*nsoset(lmax(jset))
            sgfj = first_sgf(1, jset)

            CALL ab_contract(int_sgf(sgfi:sgfi + nsgf_set(iset) - 1, sgfj:sgfj + nsgf_set(jset) - 1), &
                             int_so(starti:starti + nsoi - 1, startj:startj + nsoj - 1), &
                             sphi_so(:, sgfi:), sphi_so(:, sgfj:), nsoi, nsoj, &
                             nsgf_set(iset), nsgf_set(jset))
         END DO !jset
      END DO !iset

   END SUBROUTINE contract_so2sgf

! **************************************************************************************************
!> \brief Integrate the product of spherical gaussian orbitals with the xc kernel on the atomic grid
!> \param intso the integral in terms of spherical orbitals
!> \param fxc the xc kernel at each grid point
!> \param ikind the kind of the atom we integrate for
!> \param xas_atom_env ...
!> \param qs_env ...
!> \note Largely copied from gaVxcgb_noGC. Rewritten here because we need our own atomic grid,
!>       harmonics, basis set and we do not need the soft vxc. Could have tweaked the original, but
!>       it would have been messy. Also we do not need rho_atom (too big and fancy for us)
!>       We also go over the whole range of angular momentum l
! **************************************************************************************************
   SUBROUTINE integrate_so_prod(intso, fxc, ikind, xas_atom_env, qs_env)

      REAL(dp), DIMENSION(:, :), INTENT(INOUT)           :: intso
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: fxc
      INTEGER, INTENT(IN)                                :: ikind
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_so_prod', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, ia, ic, icg, ipgf1, ipgf2, iset1, iset2, iso, iso1, iso2, l, ld, lmax12, &
         lmin12, m1, m2, max_iso_not0, max_iso_not0_local, max_s_harm, maxl, maxso, n1, n2, na, &
         ngau1, ngau2, nngau1, nr, nset, size1
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cg_n_list
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cg_list
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: g1, g2
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: gfxcg, gg, matso
      REAL(dp), DIMENSION(:, :), POINTER                 :: zet
      REAL(dp), DIMENSION(:, :, :), POINTER              :: my_CG
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (grid_atom, harmonics, ri_basis, qs_kind_set, lmax, lmin, npgf, zet, my_CG)

!  Initialization
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=ri_basis, basis_type="RI_XAS")
      grid_atom => xas_atom_env%grid_atom_set(ikind)%grid_atom
      harmonics => xas_atom_env%harmonics_atom_set(ikind)%harmonics_atom

      CALL get_gto_basis_set(ri_basis, lmax=lmax, lmin=lmin, maxso=maxso, maxl=maxl, npgf=npgf, &
                             nset=nset, zet=zet)

      na = grid_atom%ng_sphere
      nr = grid_atom%nr
      my_CG => harmonics%my_CG
      max_iso_not0 = harmonics%max_iso_not0
      max_s_harm = harmonics%max_s_harm
      CPASSERT(2*maxl .LE. indso(1, max_iso_not0))

      ALLOCATE (g1(nr), g2(nr), gg(nr, 0:2*maxl))
      ALLOCATE (gfxcg(na, 0:2*maxl))
      ALLOCATE (matso(nsoset(maxl), nsoset(maxl)))
      ALLOCATE (cg_list(2, nsoset(maxl)**2, max_s_harm), cg_n_list(max_s_harm))

      g1 = 0.0_dp
      g2 = 0.0_dp
      m1 = 0
!  Loop over the product of so
      DO iset1 = 1, nset
         n1 = nsoset(lmax(iset1))
         m2 = 0
         DO iset2 = 1, nset
            CALL get_none0_cg_list(my_CG, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                   max_s_harm, lmax(iset1) + lmax(iset2), cg_list, cg_n_list, &
                                   max_iso_not0_local)
            CPASSERT(max_iso_not0_local .LE. max_iso_not0)

            n2 = nsoset(lmax(iset2))
            DO ipgf1 = 1, npgf(iset1)
               ngau1 = n1*(ipgf1 - 1) + m1
               size1 = nsoset(lmax(iset1)) - nsoset(lmin(iset1) - 1)
               nngau1 = nsoset(lmin(iset1) - 1) + ngau1

               g1(:) = EXP(-zet(ipgf1, iset1)*grid_atom%rad2(1:nr))
               DO ipgf2 = 1, npgf(iset2)
                  ngau2 = n2*(ipgf2 - 1) + m2

                  g2(:) = EXP(-zet(ipgf2, iset2)*grid_atom%rad2(1:nr))
                  lmin12 = lmin(iset1) + lmin(iset2)
                  lmax12 = lmax(iset1) + lmax(iset2)

                  !get the gaussian product
                  gg = 0.0_dp
                  IF (lmin12 == 0) THEN
                     gg(:, lmin12) = g1(:)*g2(:)
                  ELSE
                     gg(:, lmin12) = grid_atom%rad2l(1:nr, lmin12)*g1(:)*g2(:)
                  END IF

                  DO l = lmin12 + 1, lmax12
                     gg(:, l) = grid_atom%rad(1:nr)*gg(:, l - 1)
                  END DO

                  ld = lmax12 + 1
                  CALL dgemm('N', 'N', na, ld, nr, 1.0_dp, fxc(1:na, 1:nr), na, gg(:, 0:lmax12), &
                             nr, 0.0_dp, gfxcg(1:na, 0:lmax12), na)

                  !integrate
                  matso = 0.0_dp
                  DO iso = 1, max_iso_not0_local
                     DO icg = 1, cg_n_list(iso)
                        iso1 = cg_list(1, icg, iso)
                        iso2 = cg_list(2, icg, iso)
                        l = indso(1, iso1) + indso(1, iso2)

                        DO ia = 1, na
                           matso(iso1, iso2) = matso(iso1, iso2) + gfxcg(ia, l)* &
                                               my_CG(iso1, iso2, iso)*harmonics%slm(ia, iso)
                        END DO !ia
                     END DO !icg
                  END DO !iso

                  !write in integral matrix
                  DO ic = nsoset(lmin(iset2) - 1) + 1, nsoset(lmax(iset2))
                     iso1 = nsoset(lmin(iset1) - 1) + 1
                     iso2 = ngau2 + ic
                     CALL daxpy(size1, 1.0_dp, matso(iso1, ic), 1, intso(nngau1 + 1, iso2), 1)
                  END DO !ic

               END DO !ipgf2
            END DO ! ipgf1
            m2 = m2 + maxso
         END DO !iset2
         m1 = m1 + maxso
      END DO !iset1

      CALL timestop(handle)

   END SUBROUTINE integrate_so_prod

! **************************************************************************************************
!> \brief This routine computes the integral of a potential V wrt the derivative of the spherical
!>        orbitals, that is <df/dx|V|dg/dy> on the atomic grid.
!> \param intso the integral in terms of the spherical orbitals (well, their derivative)
!> \param V the potential (put on the grid and weighted) to integrate
!> \param ikind the atomic kind for which we integrate
!> \param xas_atom_env ...
!> \param qs_env ...
!> \note The atomic grids are taken fron xas_atom_env and the orbitals are the normal ones. Ok since
!>       the grid and spherical harmonics for those grids are at least as good as the GAPW ones
! **************************************************************************************************
   SUBROUTINE integrate_so_dxdy_prod(intso, V, ikind, xas_atom_env, qs_env)

      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: intso
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: V
      INTEGER, INTENT(IN)                                :: ikind
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_so_dxdy_prod', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ipgf, iset, iso, j, jpgf, jset, jso, &
                                                            k, l, maxso, na, nr, nset, starti, &
                                                            startj
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: fga, fgr, r1, r2, work
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: a1, a2
      REAL(dp), DIMENSION(:, :), POINTER                 :: slm, zet
      REAL(dp), DIMENSION(:, :, :), POINTER              :: dslm_dxyz
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (grid_atom, harmonics, basis, qs_kind_set, dslm_dxyz, slm, lmin, lmax, npgf, zet)

!  Getting what we need from the atom_env
      harmonics => xas_atom_env%harmonics_atom_set(ikind)%harmonics_atom
      grid_atom => xas_atom_env%grid_atom_set(ikind)%grid_atom

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis, basis_type="ORB")

      na = grid_atom%ng_sphere
      nr = grid_atom%nr

      slm => harmonics%slm
      dslm_dxyz => harmonics%dslm_dxyz

!  Getting what we need from the orbital basis
      CALL get_gto_basis_set(gto_basis_set=basis, lmax=lmax, lmin=lmin, &
                             maxso=maxso, npgf=npgf, nset=nset, zet=zet)

!  Separate the functions into purely r and purely angular parts, compute them all
!  and use matrix mutliplication for the integral. We use f for x derivative ang g for y

      ! Separating the functions. Note that the radial part is the same for x and y derivatives
      ALLOCATE (a1(na, nset*maxso, 3), a2(na, nset*maxso, 3))
      ALLOCATE (r1(nr, nset*maxso), r2(nr, nset*maxso))
      a1 = 0.0_dp; a2 = 0.0_dp
      r1 = 0.0_dp; r2 = 0.0_dp

      DO iset = 1, nset
         DO ipgf = 1, npgf(iset)
            starti = (iset - 1)*maxso + (ipgf - 1)*nsoset(lmax(iset))
            DO iso = nsoset(lmin(iset) - 1) + 1, nsoset(lmax(iset))
               l = indso(1, iso)

               ! The x derivative of the spherical orbital, divided in angular and radial parts
               ! Two of each are needed because d/dx(r^l Y_lm) * exp(-al*r^2) + r^l Y_lm * ! d/dx(exp-al*r^2)

               ! the purely radial part of d/dx(r^l Y_lm) * exp(-al*r^2) (same for y)
               r1(1:nr, starti + iso) = grid_atom%rad(1:nr)**(l - 1)*EXP(-zet(ipgf, iset)*grid_atom%rad2(1:nr))

               ! the purely radial part of r^l Y_lm * d/dx(exp-al*r^2) (same for y)
               r2(1:nr, starti + iso) = -2.0_dp*zet(ipgf, iset)*grid_atom%rad(1:nr)**(l + 1) &
                                        *EXP(-zet(ipgf, iset)*grid_atom%rad2(1:nr))

               DO i = 1, 3
                  ! the purely angular part of d/dx(r^l Y_lm) * exp(-al*r^2)
                  a1(1:na, starti + iso, i) = dslm_dxyz(i, 1:na, iso)

                  ! the purely angular part of r^l Y_lm * d/dx(exp-al*r^2)
                  a2(1:na, starti + iso, i) = harmonics%a(i, 1:na)*slm(1:na, iso)
               END DO

            END DO !iso
         END DO !ipgf
      END DO !iset

      ! Do the integration in terms of so using matrix products
      intso = 0.0_dp
      ALLOCATE (fga(na, 1))
      ALLOCATE (fgr(nr, 1))
      ALLOCATE (work(na, 1))
      fga = 0.0_dp; fgr = 0.0_dp; work = 0.0_dp

      DO iset = 1, nset
         DO jset = 1, nset
            DO ipgf = 1, npgf(iset)
               starti = (iset - 1)*maxso + (ipgf - 1)*nsoset(lmax(iset))
               DO jpgf = 1, npgf(jset)
                  startj = (jset - 1)*maxso + (jpgf - 1)*nsoset(lmax(jset))

                  DO i = 1, 3
                     j = MOD(i, 3) + 1
                     k = MOD(i + 1, 3) + 1

                     DO iso = nsoset(lmin(iset) - 1) + 1, nsoset(lmax(iset))
                        DO jso = nsoset(lmin(jset) - 1) + 1, nsoset(lmax(jset))

                           !Two component per function => 4 terms in total

                           ! take r1*a1(j) * V * r1*a1(k)
                           fgr(1:nr, 1) = r1(1:nr, starti + iso)*r1(1:nr, startj + jso)
                           fga(1:na, 1) = a1(1:na, starti + iso, j)*a1(1:na, startj + jso, k)

                           CALL dgemm('N', 'N', na, 1, nr, 1.0_dp, V, na, fgr, nr, 0.0_dp, work, na)
                           CALL dgemm('T', 'N', 1, 1, na, 1.0_dp, work, na, fga, na, 0.0_dp, &
                                      intso(starti + iso, startj + jso, i), 1)

                           ! add r1*a1(j) * V * r2*a2(k)
                           fgr(1:nr, 1) = r1(1:nr, starti + iso)*r2(1:nr, startj + jso)
                           fga(1:na, 1) = a1(1:na, starti + iso, j)*a2(1:na, startj + jso, k)

                           CALL dgemm('N', 'N', na, 1, nr, 1.0_dp, V, na, fgr, nr, 0.0_dp, work, na)
                           CALL dgemm('T', 'N', 1, 1, na, 1.0_dp, work, na, fga, na, 1.0_dp, &
                                      intso(starti + iso, startj + jso, i), 1)

                           ! add r2*a2(j) * V * r1*a1(k)
                           fgr(1:nr, 1) = r2(1:nr, starti + iso)*r1(1:nr, startj + jso)
                           fga(1:na, 1) = a2(1:na, starti + iso, j)*a1(1:na, startj + jso, k)

                           CALL dgemm('N', 'N', na, 1, nr, 1.0_dp, V, na, fgr, nr, 0.0_dp, work, na)
                           CALL dgemm('T', 'N', 1, 1, na, 1.0_dp, work, na, fga, na, 1.0_dp, &
                                      intso(starti + iso, startj + jso, i), 1)

                           ! add the last term: r2*a2(j) * V * r2*a2(k)
                           fgr(1:nr, 1) = r2(1:nr, starti + iso)*r2(1:nr, startj + jso)
                           fga(1:na, 1) = a2(1:na, starti + iso, j)*a2(1:na, startj + jso, k)

                           CALL dgemm('N', 'N', na, 1, nr, 1.0_dp, V, na, fgr, nr, 0.0_dp, work, na)
                           CALL dgemm('T', 'N', 1, 1, na, 1.0_dp, work, na, fga, na, 1.0_dp, &
                                      intso(starti + iso, startj + jso, i), 1)

                        END DO !jso
                     END DO !iso

                  END DO !i
               END DO !jpgf
            END DO !ipgf
         END DO !jset
      END DO !iset

      DO i = 1, 3
         intso(:, :, i) = intso(:, :, i) - TRANSPOSE(intso(:, :, i))
      END DO

   END SUBROUTINE integrate_so_dxdy_prod

! **************************************************************************************************
!> \brief Computes the SOC matrix elements with respect to the ORB basis set for each atomic kind
!>        and put them as the block diagonal of dbcsr_matrix
!> \param matrix_soc the matrix where the SOC is stored
!> \param xas_atom_env ...
!> \param qs_env ...
!> \note We compute: <da_dx|V\(4c^2-2V)|db_dy> - <da_dy|V\(4c^2-2V)|db_dx>, where V is a model
!>       potential on the atomic grid. What we get is purely imaginary
! **************************************************************************************************
   SUBROUTINE integrate_soc_atoms(matrix_soc, xas_atom_env, qs_env)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_soc
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_soc_atoms', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, handle, i, iat, ikind, ir, jat, &
                                                            maxso, na, nkind, nr, nset, nsgf
      REAL(dp)                                           :: zeff
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: Vr
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: V
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: intso
      REAL(dp), DIMENSION(:, :), POINTER                 :: sphi_so
      REAL(dp), DIMENSION(:, :, :), POINTER              :: intsgf
      TYPE(cp_3d_r_p_type), DIMENSION(:), POINTER        :: int_soc
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (int_soc, basis, qs_kind_set, sphi_so, matrix_s)
      NULLIFY (particle_set)

      CALL timeset(routineN, handle)

!  Initialization
      CALL get_qs_env(qs_env, nkind=nkind, qs_kind_set=qs_kind_set, matrix_s=matrix_s, &
                      particle_set=particle_set)
      ALLOCATE (int_soc(nkind))

!  Loop over the kinds to compute the integrals
      DO ikind = 1, nkind

         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis, basis_type="ORB", zeff=zeff)
         CALL get_gto_basis_set(basis, nset=nset, maxso=maxso)
         ALLOCATE (intso(nset*maxso, nset*maxso, 3))

         ! compute the model potential on the grid
         grid => xas_atom_env%grid_atom_set(ikind)%grid_atom
         nr = grid%nr
         na = grid%ng_sphere
         ALLOCATE (Vr(nr))
         !TODO: do we need potential contribution from neighboring atoms ?!?
         CALL calculate_model_potential(Vr, grid, zeff)

         !Compute V/(4c^2-2V) and weight it
         ALLOCATE (V(na, nr))
         V = 0.0_dp
         DO ir = 1, nr
            CALL daxpy(na, Vr(ir)/(4.0_dp*c_light_au**2 - 2.0_dp*Vr(ir)), grid%weight(1:na, ir), 1, &
                       V(1:na, ir), 1)
         END DO
         DEALLOCATE (Vr)

         ! compute the integral <da_dx|...|db_dy> in terms of so
         CALL integrate_so_dxdy_prod(intso, V, ikind, xas_atom_env, qs_env)
         DEALLOCATE (V)

         ! contract in terms of sgf
         CALL get_gto_basis_set(basis, nsgf=nsgf)
         ALLOCATE (int_soc(ikind)%array(nsgf, nsgf, 3))
         intsgf => int_soc(ikind)%array
         sphi_so => xas_atom_env%orb_sphi_so(ikind)%array
         intsgf = 0.0_dp

         DO i = 1, 3
            CALL contract_so2sgf(intsgf(:, :, i), intso(:, :, i), basis, sphi_so)
         END DO

         DEALLOCATE (intso)
      END DO !ikind

!  Build the matrix_soc based on the matrix_s (but anti-symmetric)
      DO i = 1, 3
         CALL dbcsr_create(matrix_soc(i)%matrix, name="SOC MATRIX", template=matrix_s(1)%matrix, &
                           matrix_type="A")
      END DO

!  Iterate over its diagonal blocks and fill=it
      CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iat, column=jat, blk=blk)
         IF (.NOT. iat == jat) CYCLE
         ikind = particle_set(iat)%atomic_kind%kind_number

         DO i = 1, 3
            CALL dbcsr_put_block(matrix_soc(i)%matrix, iat, iat, int_soc(ikind)%array(:, :, i))
         END DO

      END DO !iat
      CALL dbcsr_iterator_stop(iter)
      DO i = 1, 3
         CALL dbcsr_finalize(matrix_soc(i)%matrix)
      END DO

! Clean-up
      DO ikind = 1, nkind
         DEALLOCATE (int_soc(ikind)%array)
      END DO
      DEALLOCATE (int_soc)

      CALL timestop(handle)

   END SUBROUTINE integrate_soc_atoms

END MODULE xas_tdp_atom
