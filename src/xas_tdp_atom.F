
!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
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
                                              cp_2d_r_p_type,&
                                              cp_3d_r_p_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              qs_control_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_create, dbcsr_finalize, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_put_block
   USE generic_os_integrals,            ONLY: int_overlap_ab_os
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                              get_number_of_lebedev_grid,&
                                              init_lebedev_grids,&
                                              lebedev_grid
   USE mathconstants,                   ONLY: dfac,&
                                              pi
   USE mathlib,                         ONLY: invmat_symm
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_sum
   USE orbital_pointers,                ONLY: indco,&
                                              indso,&
                                              nco,&
                                              ncoset,&
                                              nso,&
                                              nsoset
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: a_fine
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
   USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                              neighbor_list_set_p_type
   USE qs_o3c_methods,                  ONLY: calculate_o3c_integrals,&
                                              contract12_o3c
   USE qs_o3c_types,                    ONLY: get_o3c_iterator_info,&
                                              init_o3c_container,&
                                              o3c_container_type,&
                                              o3c_iterate,&
                                              o3c_iterator_create,&
                                              o3c_iterator_release,&
                                              o3c_iterator_type,&
                                              release_o3c_container
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE spherical_harmonics,             ONLY: clebsch_gordon,&
                                              clebsch_gordon_deallocate,&
                                              clebsch_gordon_init
   USE util,                            ONLY: get_limit,&
                                              locate
   USE xas_tdp_types,                   ONLY: xas_atom_env_type,&
                                              xas_tdp_control_type,&
                                              xas_tdp_env_type
   USE xas_tdp_utils,                   ONLY: build_xas_tdp_3c_nl
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
!> \brief Initializes a xas_atom_env type given the qs_env
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

      ALLOCATE (xas_atom_env%grid_atom_set(nkind))
      ALLOCATE (xas_atom_env%harmonics_atom_set(nkind))
      ALLOCATE (xas_atom_env%ri_sphi_so(nkind))
      ALLOCATE (xas_atom_env%orb_sphi_so(nkind))
      ALLOCATE (xas_atom_env%ri_dcoeff(natom, nspins))

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
!> \param qs_env ...
!> \param do_xc Whether the xc kernel will ne computed on the atomic grids. If not, the harmonics
!>        are built for the orbital basis for all kinds.
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

      CHARACTER(len=2)                                   :: symbol
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
         DO iso1 = nsoset(lc1-1)+1, nsoset(lc1)
            l1 = indso(1, iso1)
            m1 = indso(2, iso1)
            DO lc2 = 0, maxlgto
               DO iso2 = nsoset(lc2-1)+1, nsoset(lc2)
                  l2 = indso(1, iso2)
                  m2 = indso(2, iso2)
                  CALL clebsch_gordon(l1, m1, l2, m2, rga)
                  IF (l1+l2 > llmax) THEN
                     l1l2 = llmax
                  ELSE
                     l1l2 = l1+l2
                  END IF
                  mp = m1+m2
                  mm = m1-m2
                  IF (m1*m2 < 0 .OR. (m1*m2 == 0 .AND. (m1 < 0 .OR. m2 < 0))) THEN
                     mp = -ABS(mp)
                     mm = -ABS(mm)
                  ELSE
                     mp = ABS(mp)
                     mm = ABS(mm)
                  END IF
                  DO lp = MOD(l1+l2, 2), l1l2, 2
                     il = lp/2+1
                     IF (ABS(mp) <= lp) THEN
                     IF (mp >= 0) THEN
                        iso = nsoset(lp-1)+lp+1+mp
                     ELSE
                        iso = nsoset(lp-1)+lp+1-ABS(mp)
                     END IF
                     my_CG(iso1, iso2, iso) = rga(il, 1)
                     ENDIF
                     IF (mp /= mm .AND. ABS(mm) <= lp) THEN
                     IF (mm >= 0) THEN
                        iso = nsoset(lp-1)+lp+1+mm
                     ELSE
                        iso = nsoset(lp-1)+lp+1-ABS(mm)
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
         CALL get_qs_kind(qs_kind_set(ikind), ngrid_rad=nr, ngrid_ang=na, element_symbol=symbol)

         !take the grid dimension given as input, if none, take the GAPW ones above
         IF (ANY(grid_info == symbol)) THEN
            DO igrid = 1, SIZE(grid_info, 1)
               IF (grid_info(igrid, 1) == symbol) THEN

                  !hack to convert string into integer
                  READ (grid_info(igrid, 2), *, iostat=stat) na
                  IF (stat .NE. 0) CPABORT("The 'na' value for the GRID keyword must be an integer")
                  READ (grid_info(igrid, 3), *, iostat=stat) nr
                  IF (stat .NE. 0) CPABORT("The 'nr' value for the GRID keyword must be an integer")
                  CYCLE
               END IF
            END DO
         END IF

         ll = get_number_of_lebedev_grid(n=na)
         na = lebedev_grid(ll)%n
         la = lebedev_grid(ll)%l
         grid_atom%ng_sphere = na
         grid_atom%nr = nr

         IF (llmax > la) THEN
            CPWARN("A Lebedev grid for XAS TDP is built for a too low angular momentum l")
            CPWARN("To fix this (l too low), increase the value of 'na' inder the GRID keyword.")
         END IF

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
      llmax_p1 = SIZE(grid_atom%rad2l, 2)-1

!     Find the index corresponding to the limiting radius (small ir => large radius)
      DO ir = 1, nr
         IF (grid_atom%rad(ir) < max_radius) THEN
            first_ir = ir
            EXIT
         END IF
      END DO
      new_nr = nr-first_ir+1

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
            start_s = (ipgf-1)*nsoset(lmax(iset))
            start_c = (ipgf-1)*ncoset(lmax(iset))

            DO l = lmin(iset), lmax(iset)
               DO iso = 1, nso(l)
                  DO ico = 1, nco(l)

                     lx = indco(1, ico+ncoset(l-1))
                     ly = indco(2, ico+ncoset(l-1))
                     lz = indco(3, ico+ncoset(l-1))

                     factor = orbtramat(l)%s2c(iso, ico) &
                              *SQRT(4.0_dp*pi/dfac(2*l+1)*dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1))

                     CALL daxpy(nsgf_set(iset), factor, &
                                sphi(start_c+ncoset(l-1)+ico, sgfi:sgfi+nsgf_set(iset)-1), 1, &
                                sphi_so(start_s+nsoset(l-1)+iso, sgfi:sgfi+nsgf_set(iset)-1), 1)

                  END DO !ico
               END DO !iso
            END DO !l
         END DO !ipgf
      END DO !iset

   END SUBROUTINE compute_sphi_so

! **************************************************************************************************
!> \brief Find the neighbor of each excited atom by looping over the orbital overlap matrix
!>        Also returns an array containing all excited atom index and all neighbors, uniquely
!> \param all_ri_atoms the indices of all atoms where an RI expansion will be done
!> \param all_ri_kinds the kinds of all RI atoms
!> \param xas_atom_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE find_exat_neighbors(all_ri_atoms, all_ri_kinds, xas_atom_env, qs_env)

      INTEGER, DIMENSION(:), POINTER                     :: all_ri_atoms, all_ri_kinds
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'find_exat_neighbors', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, i, iat, iblk, iex, jblk, mepos, &
                                                            natom, nb, nex
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inb, who_is_there
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: n_neighbors
      INTEGER, DIMENSION(:), POINTER                     :: excited_atoms
      TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER        :: exat_neighbors
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (matrix_s, excited_atoms, para_env, exat_neighbors, particle_set)

      ! Initialization
      excited_atoms => xas_atom_env%excited_atoms
      nex = SIZE(excited_atoms)
      CALL get_qs_env(qs_env, matrix_s=matrix_s, para_env=para_env, natom=natom, particle_set=particle_set)
      ALLOCATE (xas_atom_env%exat_neighbors(nex))
      exat_neighbors => xas_atom_env%exat_neighbors
      mepos = para_env%mepos

      ! First loop over S => count the number of neighbors
      ALLOCATE (n_neighbors(nex, 0:para_env%num_pe-1))
      n_neighbors = 0

      CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)

         !avoid having oneself as neighbor
         IF (iblk == jblk) CYCLE

         IF (ANY(excited_atoms == iblk)) THEN
            iex = locate(excited_atoms, iblk)
            n_neighbors(iex, mepos) = n_neighbors(iex, mepos)+1
         END IF
         IF (ANY(excited_atoms == jblk)) THEN
            iex = locate(excited_atoms, jblk)
            n_neighbors(iex, mepos) = n_neighbors(iex, mepos)+1
         END IF

      END DO !iter
      CALL dbcsr_iterator_stop(iter)
      CALL mp_sum(n_neighbors, para_env%group)

      ! Allocate the exat_neighbors array at the correct length
      DO iex = 1, nex
         ALLOCATE (exat_neighbors(iex)%array(SUM(n_neighbors(iex, :))))
         exat_neighbors(iex)%array = 0
      END DO

      ! Loop a second time over S, this time fill the neighbors details
      CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
      ALLOCATE (inb(nex))
      inb = 1
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iblk, column=jblk, blk=blk)
         IF (iblk == jblk) CYCLE

         IF (ANY(excited_atoms == iblk)) THEN
            iex = locate(excited_atoms, iblk)
            exat_neighbors(iex)%array(SUM(n_neighbors(iex, 0:mepos-1))+inb(iex)) = jblk
            inb(iex) = inb(iex)+1
         END IF
         IF (ANY(excited_atoms == jblk)) THEN
            iex = locate(excited_atoms, jblk)
            exat_neighbors(iex)%array(SUM(n_neighbors(iex, 0:mepos-1))+inb(iex)) = iblk
            inb(iex) = inb(iex)+1
         END IF

      END DO !iter
      CALL dbcsr_iterator_stop(iter)

      ! Get all the indices in a sigle array for RI expansion
      ALLOCATE (who_is_there(natom))
      who_is_there = 0

      ! Make sure that all process have the whole info
      DO iex = 1, nex
         CALL mp_sum(exat_neighbors(iex)%array, para_env%group)
         who_is_there(iex) = 1
         DO nb = 1, SIZE(exat_neighbors(iex)%array)
            who_is_there(exat_neighbors(iex)%array(nb)) = 1
         END DO
      END DO

      ALLOCATE (all_ri_atoms(natom), all_ri_kinds(natom))
      i = 0
      DO iat = 1, natom
         IF (who_is_there(iat) == 1) THEN
            i = i+1
            all_ri_atoms(i) = iat
            all_ri_kinds(i) = particle_set(iat)%atomic_kind%kind_number
         END IF
      END DO
      CALL reallocate(all_ri_atoms, 1, i)
      CALL reallocate(all_ri_kinds, 1, i)

   END SUBROUTINE find_exat_neighbors

! **************************************************************************************************
!> \brief Compute the coefficients to project the density on the RI_XAS basis set
!> \param xas_atom_env ...
!> \param qs_env ...
!> \note The density is n = sum_ab P_ab*phi_a*phi_b, the RI basis covers the products of orbital sgfs
!>       => n = sum_ab sum_cd P_ab (phi_a phi_b xi_c) S_cd^-1 xi_d
!>            = sum_d coeff_d xi_d , where xi are the RI basis func.
!>       The 3-center overlap is needed, as well as the inverse of the RI basis overlap. Note that
!>       the RI basis functions are only centered on excited atoms and their neighbors
! **************************************************************************************************
   SUBROUTINE calculate_density_coeffs(xas_atom_env, qs_env)

      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_density_coeffs', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iat, ikind, ilist, ispin, katom, &
                                                            kkind, mepos, natom, nkind, nsgf_kind, &
                                                            nspins, nthread
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nsgf
      INTEGER, DIMENSION(:), POINTER                     :: all_ri_atoms, all_ri_kinds
      REAL(dp)                                           :: dmax, factor
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: work
      REAL(dp), DIMENSION(:, :), POINTER                 :: tvec
      TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER        :: sinv
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_orb, basis_set_ri
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: ac_list, sab_orb
      TYPE(o3c_container_type), POINTER                  :: o3c
      TYPE(o3c_iterator_type)                            :: o3c_iterator
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho

      NULLIFY (qs_kind_set, basis_set_ri, basis_set_orb, ac_list, rho, rho_ao, sinv, ri_basis)
      NULLIFY (o3c, tvec, particle_set, para_env, all_ri_atoms)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, nkind=nkind, sab_orb=sab_orb, qs_kind_set=qs_kind_set, rho=rho, &
                      natom=natom, particle_set=particle_set)
      nspins = SIZE(xas_atom_env%ri_dcoeff, 2)

!  Get the neighbors of the excited atoms and get all the atoms where RI needs to be done
      CALL find_exat_neighbors(all_ri_atoms, all_ri_kinds, xas_atom_env, qs_env)

!  Create the needed neighbor list and basis set lists.
      ALLOCATE (basis_set_ri(nkind))
      ALLOCATE (basis_set_orb(nkind))
      CALL basis_set_list_setup(basis_set_ri, "RI_XAS", qs_kind_set)
      CALL basis_set_list_setup(basis_set_orb, "ORB", qs_kind_set)
      CALL build_xas_tdp_3c_nl(ac_list, basis_set_orb, basis_set_ri, "OVERLAP", all_ri_atoms, &
                               all_ri_kinds, qs_env)

!  Create the inverse RI overlap matrices.
      ALLOCATE (sinv(nkind))
      DO ikind = 1, nkind
         NULLIFY (sinv(ikind)%array)
         IF (.NOT. ANY(all_ri_kinds == ikind)) CYCLE
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=ri_basis, basis_type="RI_XAS", nsgf=nsgf_kind)
         ALLOCATE (sinv(ikind)%array(nsgf_kind, nsgf_kind))

         CALL int_overlap_ab_os(sab=sinv(ikind)%array, rab=(/0.0_dp, 0.0_dp, 0.0_dp/), fba=ri_basis, &
                                fbb=ri_basis, calculate_forces=.FALSE., debug=.FALSE., dmax=dmax)
         CALL invmat_symm(sinv(ikind)%array)
      END DO

!  Create the 3-center overlap integrals for the RI basis
      ALLOCATE (o3c)
      CALL init_o3c_container(o3c, nspins, basis_set_orb, basis_set_orb, basis_set_ri, sab_orb, ac_list)
      CALL calculate_o3c_integrals(o3c)

!  Contract the 3-center tensor with the density matrix
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CALL contract12_o3c(o3c, rho_ao)

      !keep in mind that double occupation is included in rho_ao in case of closed-shell
      factor = 0.5_dp; IF (nspins == 2) factor = 1.0_dp

!  Allocate space for the projected density coefficients.
      ALLOCATE (nsgf(natom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=nsgf, basis=basis_set_ri)
      DO ispin = 1, nspins
         DO iat = 1, natom
            NULLIFY (xas_atom_env%ri_dcoeff(iat, ispin)%array)
            IF (.NOT. ANY(all_ri_atoms == iat)) CYCLE
            ikind = particle_set(iat)%atomic_kind%kind_number
            ALLOCATE (xas_atom_env%ri_dcoeff(iat, ispin)%array(nsgf(iat)))
            xas_atom_env%ri_dcoeff(iat, ispin)%array = 0.0_dp
         END DO
      END DO

!  Iterate over the o3c to access the contracted tvec = sum_ab P_ab*(abc).
      nthread = 1
!$    nthread = omp_get_max_threads()
      CALL o3c_iterator_create(o3c, o3c_iterator, nthread=nthread)

!$OMP PARALLEL DEFAULT(NONE)&
!$OMP SHARED (o3c_iterator,xas_atom_env,nsgf,sinv,nspins,factor)&
!$OMP PRIVATE (mepos,katom,kkind,tvec,work)

      mepos = 0
!$    mepos = omp_get_thread_num()

      DO WHILE (o3c_iterate(o3c_iterator, mepos=mepos) == 0)
         CALL get_o3c_iterator_info(o3c_iterator, mepos=mepos, katom=katom, kkind=kkind, tvec=tvec)

!        The coefficient for atom katom of kkind are tvec*S^-1
         ALLOCATE (work(nsgf(katom), nspins))
         CALL dgemm('N', 'N', nsgf(katom), nspins, nsgf(katom), 1.0_dp, sinv(kkind)%array, nsgf(katom), &
                    tvec, nsgf(katom), 0.0_dp, work, nsgf(katom))
!$OMP CRITICAL(spin1)
         xas_atom_env%ri_dcoeff(katom, 1)%array(:) = xas_atom_env%ri_dcoeff(katom, 1)%array(:)+factor*work(:, 1)
!$OMP END CRITICAL(spin1)

         IF (nspins == 2) THEN
!$OMP CRITICAL(spin2)
            xas_atom_env%ri_dcoeff(katom, 2)%array(:) = xas_atom_env%ri_dcoeff(katom, 2)%array(:)+work(:, 2)
!$OMP END CRITICAL(spin2)
         END IF

         DEALLOCATE (work)
      END DO
!$OMP END PARALLEL
      CALL o3c_iterator_release(o3c_iterator)

      CALL get_qs_env(qs_env, para_env=para_env)
      DO ispin = 1, nspins
         DO iat = 1, natom
            IF (.NOT. ANY(all_ri_atoms == iat)) CYCLE
            CALL mp_sum(xas_atom_env%ri_dcoeff(iat, ispin)%array, para_env%group)
         END DO !iat
      END DO !ispin

!  some clean-up
      DO ilist = 1, SIZE(ac_list)
         CALL deallocate_neighbor_list_set(ac_list(ilist)%neighbor_list_set)
      END DO
      DO ikind = 1, nkind
         IF (ASSOCIATED(sinv(ikind)%array)) DEALLOCATE (sinv(ikind)%array)
      END DO
      DEALLOCATE (sinv)
      CALL release_o3c_container(o3c)
      DEALLOCATE (basis_set_ri, basis_set_orb, o3c, ac_list, all_ri_atoms, all_ri_kinds)

      CALL timestop(handle)

   END SUBROUTINE calculate_density_coeffs

! **************************************************************************************************
!> \brief Evaluates the density on a given atomic grid
!> \param rho_set where the densities are stored
!> \param atom_index the index of the atom in question
!> \param atom_kind the kind of the atom in question
!> \param xas_atom_env ...
!> \param qs_env ...
!> \note The density is expressed as n = sum_d coeff_d*xi_d. Knowing the coordinate of each grid
!>       grid point, one can simply evaluate xi_d(r)
! **************************************************************************************************
   SUBROUTINE put_density_on_atomic_grid(rho_set, atom_index, atom_kind, xas_atom_env, qs_env)

      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      INTEGER, INTENT(IN)                                :: atom_index, atom_kind
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'put_density_on_atomic_grid', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ipgf, ir, iset, isgf, iso, l, n, &
                                                            na, nr, nset, nspins, sgfi, start
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf, nsgf_set
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: sgf, so
      REAL(dp), DIMENSION(:, :), POINTER                 :: ri_sphi_so, zet
      REAL(dp), DIMENSION(:, :, :), POINTER              :: rhoa, rhob
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: ri_basis
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (grid_atom, ri_basis, qs_kind_set, harmonics, lmax, npgf, zet, nsgf_set, ri_sphi_so)
      NULLIFY (lmin, first_sgf, rhoa, rhob)

      CALL timeset(routineN, handle)

!  Strategy: it makes sense to evaluate the spherical orbital on the grid (because of symmetry)
!            From there, one can directly contract into sgf using ri_sphi_so and then take the weight

!  Generalities
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(atom_kind), basis_set=ri_basis, basis_type="RI_XAS")
      CALL get_gto_basis_set(ri_basis, lmax=lmax, npgf=npgf, zet=zet, nset=nset, nsgf_set=nsgf_set, &
                             first_sgf=first_sgf, lmin=lmin)

!  Get the grid, the harmonics and the info we need from it
      grid_atom => xas_atom_env%grid_atom_set(atom_kind)%grid_atom
      harmonics => xas_atom_env%harmonics_atom_set(atom_kind)%harmonics_atom
      CPASSERT(ASSOCIATED(grid_atom))
      CPASSERT(ASSOCIATED(harmonics))
      nr = grid_atom%nr
      na = grid_atom%ng_sphere
      n = na*nr
      nspins = SIZE(xas_atom_env%ri_dcoeff, 2)

!  Point to the rho_set densities
      rhoa => rho_set%rhoa
      rhob => rho_set%rhob
      rhoa = 0.0_dp; rhob = 0.0_dp; 
!  Loop over the sets and pgfs of the basis set
      DO iset = 1, nset

         !allocate space to store the spherical orbitals on the grid
         ALLOCATE (so(na, nr, npgf(iset)*nsoset(lmax(iset))))
         so = 0.0_dp

         DO ipgf = 1, npgf(iset)
            start = (ipgf-1)*nsoset(lmax(iset))

            !loop over the spherical gaussian orbitals
            DO iso = nsoset(lmin(iset)-1)+1, nsoset(lmax(iset))
               l = indso(1, iso)
               DO ir = 1, nr
                  !compute: g = r**l * Y_lm * exp(-zet* r**2)
                  CALL daxpy(na, EXP(-zet(ipgf, iset)*grid_atom%rad2(ir))*grid_atom%rad(ir)**l, &
                             harmonics%slm(1:na, iso), 1, so(1:na, ir, start+iso), 1)
               END DO
            END DO !iso

         END DO !ipgf

         !contract the so to get the sgf
         ALLOCATE (sgf(na, nr, nsgf_set(iset)))
         sgf = 0.0_dp
         sgfi = first_sgf(1, iset)-1

         DO isgf = 1, nsgf_set(iset)

            ri_sphi_so => xas_atom_env%ri_sphi_so(atom_kind)%array
            DO ipgf = 1, npgf(iset)
               start = (ipgf-1)*nsoset(lmax(iset))
               DO iso = nsoset(lmin(iset)-1)+1, nsoset(lmax(iset))
                  CALL daxpy(n, ri_sphi_so(start+iso, sgfi+isgf), so(:, :, start+iso), 1, sgf(:, :, isgf), 1)
               END DO !iso
            END DO !ipgf

            !put the sgf on the grid with the approriate coefficients and sum
            CALL daxpy(n, xas_atom_env%ri_dcoeff(atom_index, 1)%array(sgfi+isgf), sgf(:, :, isgf), &
                       1, rhoa(:, :, 1), 1)

            IF (nspins == 2) THEN
               CALL daxpy(n, xas_atom_env%ri_dcoeff(atom_index, 2)%array(sgfi+isgf), sgf(:, :, isgf), &
                          1, rhob(:, :, 1), 1)
            END IF

         END DO !isgf

         DEALLOCATE (so, sgf)
      END DO !iset

      ! Treat spin restricted case
      IF (nspins == 1) THEN
         CALL dcopy(n, rhoa(:, :, 1), 1, rhob(:, :, 1), 1)
      END IF

      CALL timestop(handle)

   END SUBROUTINE put_density_on_atomic_grid

! **************************************************************************************************
!> \brief Adds the density of a given source atom with source kind (with ri_dcoeff) on the atomic
!>        grid belonging to another target atom of target kind. The evaluations of the basis
!>        function first requires the evaluation of the x,y,z coordinates on each grid point of
!>        target atom wrt to the position of source atom
!> \param rho_set where the densities are stored
!> \param source_iat the index of the source atom
!> \param source_ikind the kind of the source atom
!> \param target_iat the index of the target atom
!> \param target_ikind the kind of the target atom
!> \param xas_atom_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE put_density_on_other_grid(rho_set, source_iat, source_ikind, target_iat, &
                                        target_ikind, xas_atom_env, qs_env)

      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      INTEGER, INTENT(IN)                                :: source_iat, source_ikind, target_iat, &
                                                            target_ikind
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'put_density_on_other_grid', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ia, ico, ipgf, ir, iset, isgf, &
                                                            lx, ly, lz, n, na, nr, nset, nspins, &
                                                            sgfi, start
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf, nsgf_set
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: co, pos, sgf
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
      nr = grid_atom%nr
      na = grid_atom%ng_sphere
      n = na*nr
      nspins = SIZE(xas_atom_env%ri_dcoeff, 2)

      !  Point to the rho_set densities
      rhoa => rho_set%rhoa
      rhob => rho_set%rhob

      !  Need the source-target position vector
      rs = pbc(particle_set(source_iat)%r, cell)
      rt = pbc(particle_set(target_iat)%r, cell)
      rst = pbc(rt-rs, cell)

      ! Precompute the positions on the target grid
      ALLOCATE (pos(na, nr, 4))
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(NONE), &
!$OMP SHARED(na,nr,pos,harmonics,grid_atom,rst), &
!$OMP PRIVATE(ia,ir)
      DO ir = 1, nr
         DO ia = 1, na
            pos(ia, ir, 1:3) = harmonics%a(:, ia)*grid_atom%rad(ir)+rst
            pos(ia, ir, 4) = pos(ia, ir, 1)**2+pos(ia, ir, 2)**2+pos(ia, ir, 3)**2
         END DO
      END DO
!$OMP END PARALLEL DO

      ! Loop over the cartesian gaussian functions and evaluate them
      DO iset = 1, nset

         !allocate space to store the cartesian orbtial on the grid
         ALLOCATE (co(na, nr, npgf(iset)*ncoset(lmax(iset))))

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP SHARED(co,npgf,ncoset,lmax,lmin,indco,pos,zet,iset,na,nr), &
!$OMP PRIVATE(ipgf,start,ico,lx,ly,lz,ia,ir)

!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
         DO ir = 1, nr
            DO ia = 1, na
               co(ia, ir, :) = 0.0_dp
            END DO
         END DO
!$OMP END DO NOWAIT

         DO ipgf = 1, npgf(iset)
            start = (ipgf-1)*ncoset(lmax(iset))

            !loop over the cartesian orbitals
            DO ico = ncoset(lmin(iset)-1)+1, ncoset(lmax(iset))
               lx = indco(1, ico)
               ly = indco(2, ico)
               lz = indco(3, ico)

               ! compute g = x**lx * y**ly * z**lz * exp(-zet * r**2)
               !co(:,:,start+ico) = pos(:,:,1)**lx * pos(:,:,2)**ly * pos(:,:,3)**lz &
               !                    * EXP(-zet(ipgf, iset) * pos(:,:,4))
!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
               DO ir = 1, nr
                  DO ia = 1, na
                     co(ia, ir, start+ico) = pos(ia, ir, 1)**lx*pos(ia, ir, 2)**ly*pos(ia, ir, 3)**lz &
                                             *EXP(-zet(ipgf, iset)*pos(ia, ir, 4))
                  END DO
               END DO
!$OMP END DO NOWAIT

            END DO !ico

         END DO !ipgf

!$OMP END PARALLEL

         !contract the co into sgf
         ALLOCATE (sgf(na, nr, nsgf_set(iset)))
         sgf = 0.0_dp
         sgfi = first_sgf(1, iset)-1

         DO isgf = 1, nsgf_set(iset)
            DO ipgf = 1, npgf(iset)
               start = (ipgf-1)*ncoset(lmax(iset))
               DO ico = ncoset(lmin(iset)-1)+1, ncoset(lmax(iset))
                  CALL daxpy(n, ri_sphi(start+ico, sgfi+isgf), co(:, :, start+ico), 1, sgf(:, :, isgf), 1)
               END DO !ico
            END DO !ipgf

            !add the density to the grid
            CALL daxpy(n, xas_atom_env%ri_dcoeff(source_iat, 1)%array(sgfi+isgf), sgf(:, :, isgf), 1, rhoa(:, :, 1), 1)

            IF (nspins == 2) THEN
               CALL daxpy(n, xas_atom_env%ri_dcoeff(source_iat, 2)%array(sgfi+isgf), sgf(:, :, isgf), 1, rhob(:, :, 1), 1)
            END IF
         END DO !isgf
         DEALLOCATE (co, sgf)
      END DO !iset

      IF (nspins == 1) THEN
         CALL dcopy(n, rhoa(:, :, 1), 1, rhob(:, :, 1), 1)
      END IF

      CALL timestop(handle)

   END SUBROUTINE put_density_on_other_grid

! **************************************************************************************************
!> \brief Integrate the xc kernel as a function of r on the atomic grids for the RI_XAS basis
!> \param int_fxc the global array containing the (P|fxc|Q) integrals, for all spin configurations
!> \param xas_atom_env ...
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

      INTEGER                                            :: bo(2), handle, i, iatom, iex, ikind, &
                                                            inb, mepos, na, natom, nb, nbk, &
                                                            nex_atom, nr, num_pe
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(:), POINTER                     :: exat_neighbors
      LOGICAL                                            :: do_sc, do_sf
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

!  Distribute the excited atoms over the processors
      num_pe = para_env%num_pe
      mepos = para_env%mepos
      bo = get_limit(nex_atom, num_pe, mepos)

!  Looping over the excited atoms
      DO iex = bo(1), bo(2)

         iatom = xas_atom_env%excited_atoms(iex)
         ikind = particle_set(iatom)%atomic_kind%kind_number
         exat_neighbors => xas_atom_env%exat_neighbors(iex)%array

!     General grid/basis info
         na = xas_atom_env%grid_atom_set(ikind)%grid_atom%ng_sphere
         nr = xas_atom_env%grid_atom_set(ikind)%grid_atom%nr

!     Get the xc functional subsection
         xc_functionals => section_vals_get_subs_vals(input, "DFT%XAS_TDP%KERNEL%XC_FUNCTIONAL")

!     Creating a xc_rho_set to store the density and dset for the kernel
         bounds(1:2, 1:3) = 1
         bounds(2, 1) = na
         bounds(2, 2) = nr
         CALL xc_rho_set_create(rho_set=rho_set, local_bounds=bounds, &
                                rho_cutoff=dft_control%qs_control%eps_rho_rspace)
         CALL xc_dset_create(deriv_set, local_bounds=bounds)
         ! ask for lsd in any case
         needs = xc_functionals_get_needs(xc_functionals, lsd=.TRUE., add_basic_components=.TRUE.)

         ! allocate internals of the rho_set
         CALL xc_rho_set_atom_update(rho_set, needs, nspins=2, bo=bounds)

!     Put the density on the grid (for this atom)
         CALL put_density_on_atomic_grid(rho_set, iatom, ikind, xas_atom_env, qs_env)

!     Take the neighboring atom contributions to the density
         DO inb = 1, SIZE(exat_neighbors)

            nb = exat_neighbors(inb)
            nbk = particle_set(nb)%atomic_kind%kind_number

            CALL put_density_on_other_grid(rho_set, nb, nbk, iatom, ikind, xas_atom_env, qs_env)
         END DO

         !spin-conserving
         IF (do_sc) THEN
            CALL integrate_sc_fxc(int_fxc, iatom, ikind, rho_set, deriv_set, xc_functionals, &
                                  xas_atom_env, qs_env)
         END IF

         !spin-flip
         IF (do_sf) THEN
            CALL integrate_sf_fxc(int_fxc, iatom, ikind, rho_set, deriv_set, xc_functionals, &
                                  xas_atom_env, qs_env)
         END IF

!     Clean-up
         CALL xc_dset_release(deriv_set)
         CALL xc_rho_set_release(rho_set)

      END DO !iex

      CALL timestop(handle)

   END SUBROUTINE integrate_fxc_atoms

! **************************************************************************************************
!> \brief Integrate the fxc kernel in the spin-conserving case, be it closed- or open-shell
!> \param int_fxc the array containing the (P|fxc|Q) integrals
!> \param iatom the index of the current excited atom
!> \param ikind the index of the current excited kind
!> \param rho_set the density in the atomic grid
!> \param deriv_set the set of functional derivatives
!> \param xc_functionals the input values for the xc functionals
!> \param xas_atom_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE integrate_sc_fxc(int_fxc, iatom, ikind, rho_set, deriv_set, xc_functionals, &
                               xas_atom_env, qs_env)

      TYPE(cp_2d_r_p_type), DIMENSION(:, :), POINTER     :: int_fxc
      INTEGER, INTENT(IN)                                :: iatom, ikind
      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(section_vals_type), POINTER                   :: xc_functionals
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

      ! Compute the kernel on the grid => only interested in second derivative
      CALL xc_functionals_eval(xc_functionals, lsd=.TRUE., rho_set=rho_set, deriv_set=deriv_set, &
                               deriv_order=-2)

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
!> \param xc_functionals the input values for the xc functionals
!> \param xas_atom_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE integrate_sf_fxc(int_fxc, iatom, ikind, rho_set, deriv_set, xc_functionals, &
                               xas_atom_env, qs_env)

      TYPE(cp_2d_r_p_type), DIMENSION(:, :), POINTER     :: int_fxc
      INTEGER, INTENT(IN)                                :: iatom, ikind
      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(section_vals_type), POINTER                   :: xc_functionals
      TYPE(xas_atom_env_type), POINTER                   :: xas_atom_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_sf_fxc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ia, ir, maxso, na, nr, nset, nsotot, &
                                                            ri_nsgf
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: fxc, int_so
      REAL(dp), DIMENSION(:, :), POINTER                 :: ri_sphi_so
      REAL(dp), DIMENSION(:, :, :), POINTER              :: rhoa, rhob
      TYPE(cp_3d_r_p_type), ALLOCATABLE, DIMENSION(:)    :: d1e
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

      ! Compute the first (only) derivatives on the grid
      CALL xc_functionals_eval(xc_functionals, lsd=.TRUE., rho_set=rho_set, deriv_set=deriv_set, &
                               deriv_order=-1)

      ALLOCATE (d1e(2))
      deriv => xc_dset_get_derivative(deriv_set, "(rhoa)")
      CALL xc_derivative_get(deriv, deriv_data=d1e(1)%array)
      deriv => xc_dset_get_derivative(deriv_set, "(rhob)")
      CALL xc_derivative_get(deriv, deriv_data=d1e(2)%array)

      !Compute the kernel on the grid. Already take weight into acocunt there
      ALLOCATE (fxc(na, nr))
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) DEFAULT(NONE), &
!$OMP SHARED(grid_atom,fxc,d1e,dft_control,na,nr,rhoa,rhob), &
!$OMP PRIVATE(ia,ir)
      DO ir = 1, nr
         DO ia = 1, na

            !Need to be careful not to devide by zero. Assume that if rhoa == rhob, then
            !dE/drhoa == dE/drhob and we send everything to zero
            IF (ABS(rhoa(ia, ir, 1)-rhob(ia, ir, 1)) > dft_control%qs_control%eps_rho_rspace) THEN
               fxc(ia, ir) = grid_atom%weight(ia, ir)/(rhoa(ia, ir, 1)-rhob(ia, ir, 1)) &
                             *(d1e(1)%array(ia, ir, 1)-d1e(2)%array(ia, ir, 1))
            ELSE
               fxc(ia, ir) = 0.0_dp
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
      CALL contract_so2sgf(int_fxc(iatom, 4)%array, int_so, ri_basis, ri_sphi_so)

   END SUBROUTINE integrate_sf_fxc

! **************************************************************************************************
!> \brief Contract spherical orbitals to spherical Gaussians (so to sgf)
!> \param int_sgf the arary with the sgf integrals
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
         starti = (iset-1)*maxso+1
         nsoi = npgf(iset)*nsoset(lmax(iset))
         sgfi = first_sgf(1, iset)

         DO jset = 1, nset
            startj = (jset-1)*maxso+1
            nsoj = npgf(jset)*nsoset(lmax(jset))
            sgfj = first_sgf(1, jset)

            CALL ab_contract(int_sgf(sgfi:sgfi+nsgf_set(iset)-1, sgfj:sgfj+nsgf_set(jset)-1), &
                             int_so(starti:starti+nsoi-1, startj:startj+nsoj-1), &
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

      nr = grid_atom%nr
      na = grid_atom%ng_sphere
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
                                   max_s_harm, lmax(iset1)+lmax(iset2), cg_list, cg_n_list, &
                                   max_iso_not0_local)
            CPASSERT(max_iso_not0_local .LE. max_iso_not0)

            n2 = nsoset(lmax(iset2))
            DO ipgf1 = 1, npgf(iset1)
               ngau1 = n1*(ipgf1-1)+m1
               size1 = nsoset(lmax(iset1))-nsoset(lmin(iset1)-1)
               nngau1 = nsoset(lmin(iset1)-1)+ngau1

               g1(1:nr) = EXP(-zet(ipgf1, iset1)*grid_atom%rad2(1:nr))
               DO ipgf2 = 1, npgf(iset2)
                  ngau2 = n2*(ipgf2-1)+m2

                  g2(1:nr) = EXP(-zet(ipgf2, iset2)*grid_atom%rad2(1:nr))
                  lmin12 = lmin(iset1)+lmin(iset2)
                  lmax12 = lmax(iset1)+lmax(iset2)

                  !get the gaussian product
                  gg = 0.0_dp
                  IF (lmin12 == 0) THEN
                     gg(1:nr, lmin12) = g1(1:nr)*g2(1:nr)
                  ELSE
                     gg(1:nr, lmin12) = grid_atom%rad2l(1:nr, lmin12)*g1(1:nr)*g2(1:nr)
                  END IF

                  DO l = lmin12+1, lmax12
                     gg(1:nr, l) = grid_atom%rad(1:nr)*gg(:, l-1)
                  END DO

                  ld = lmax12+1
                  CALL dgemm('N', 'N', na, ld, nr, 1.0_dp, fxc(1:na, 1:nr), na, gg(1:nr, 0:lmax12), &
                             nr, 0.0_dp, gfxcg(1:na, 0:lmax12), na)

                  !integrate
                  matso = 0.0_dp
                  DO iso = 1, max_iso_not0_local
                     DO icg = 1, cg_n_list(iso)
                        iso1 = cg_list(1, icg, iso)
                        iso2 = cg_list(2, icg, iso)
                        l = indso(1, iso1)+indso(1, iso2)

                        DO ia = 1, na
                           matso(iso1, iso2) = matso(iso1, iso2)+gfxcg(ia, l)* &
                                               my_CG(iso1, iso2, iso)*harmonics%slm(ia, iso)
                        END DO !ia
                     END DO !icg
                  END DO !iso

                  !write in integral matrix
                  DO ic = nsoset(lmin(iset2)-1)+1, nsoset(lmax(iset2))
                     iso1 = nsoset(lmin(iset1)-1)+1
                     iso2 = ngau2+ic
                     CALL daxpy(size1, 1.0_dp, matso(iso1, ic), 1, intso(nngau1+1, iso2), 1)
                  END DO !ic

               END DO !ipgf2
            END DO ! ipgf1
            m2 = m2+maxso
         END DO !iset2
         m1 = m1+maxso
      END DO !iset1

      CALL timestop(handle)

   END SUBROUTINE integrate_so_prod

! **************************************************************************************************
!> \brief This routine computes the integral of a potential V wrt the derivitive of the spherical
!>        orbitals, that is <df/dx|V|dg/dy> on the atomic grid.
!> \param intso the integral in terms of the spherical orbitals (well, their derivitive)
!> \param V the potential (put on the grid and wieghted) to integrate
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
            starti = (iset-1)*maxso+(ipgf-1)*nsoset(lmax(iset))
            DO iso = nsoset(lmin(iset)-1)+1, nsoset(lmax(iset))
               l = indso(1, iso)

               ! The x derivitive of the spherical orbital, divided in angular and radial parts
               ! Two of each are needed because d/dx(r^l Y_lm) * exp(-al*r^2) + r^l Y_lm * ! d/dx(exp-al*r^2)

               ! the purely radial part of d/dx(r^l Y_lm) * exp(-al*r^2) (same for y)
               r1(1:nr, starti+iso) = grid_atom%rad(1:nr)**(l-1)*EXP(-zet(ipgf, iset)*grid_atom%rad2(1:nr))

               ! the purely radial part of r^l Y_lm * d/dx(exp-al*r^2) (same for y)
               r2(1:nr, starti+iso) = -2.0_dp*zet(ipgf, iset)*grid_atom%rad(1:nr)**(l+1) &
                                      *EXP(-zet(ipgf, iset)*grid_atom%rad2(1:nr))

               DO i = 1, 3
                  ! the purely angular part of d/dx(r^l Y_lm) * exp(-al*r^2)
                  a1(1:na, starti+iso, i) = dslm_dxyz(i, 1:na, iso)

                  ! the purely angular part of r^l Y_lm * d/dx(exp-al*r^2)
                  a2(1:na, starti+iso, i) = harmonics%a(i, 1:na)*slm(1:na, iso)
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
               starti = (iset-1)*maxso+(ipgf-1)*nsoset(lmax(iset))
               DO jpgf = 1, npgf(jset)
                  startj = (jset-1)*maxso+(jpgf-1)*nsoset(lmax(jset))

                  DO i = 1, 3
                     j = MOD(i, 3)+1
                     k = MOD(i+1, 3)+1

                     DO iso = nsoset(lmin(iset)-1)+1, nsoset(lmax(iset))
                        DO jso = nsoset(lmin(jset)-1)+1, nsoset(lmax(jset))

                           !Two component per function => 4 terms in total

                           ! take r1*a1(j) * V * r1*a1(k)
                           fgr(1:nr, 1) = r1(1:nr, starti+iso)*r1(1:nr, startj+jso)
                           fga(1:na, 1) = a1(1:na, starti+iso, j)*a1(1:na, startj+jso, k)

                           CALL dgemm('N', 'N', na, 1, nr, 1.0_dp, V, na, fgr, nr, 0.0_dp, work, na)
                           CALL dgemm('T', 'N', 1, 1, na, 1.0_dp, work, na, fga, na, 0.0_dp, &
                                      intso(starti+iso, startj+jso, i), 1)

                           ! add r1*a1(j) * V * r2*a2(k)
                           fgr(1:nr, 1) = r1(1:nr, starti+iso)*r2(1:nr, startj+jso)
                           fga(1:na, 1) = a1(1:na, starti+iso, j)*a2(1:na, startj+jso, k)

                           CALL dgemm('N', 'N', na, 1, nr, 1.0_dp, V, na, fgr, nr, 0.0_dp, work, na)
                           CALL dgemm('T', 'N', 1, 1, na, 1.0_dp, work, na, fga, na, 1.0_dp, &
                                      intso(starti+iso, startj+jso, i), 1)

                           ! add r2*a2(j) * V * r1*a1(k)
                           fgr(1:nr, 1) = r2(1:nr, starti+iso)*r1(1:nr, startj+jso)
                           fga(1:na, 1) = a2(1:na, starti+iso, j)*a1(1:na, startj+jso, k)

                           CALL dgemm('N', 'N', na, 1, nr, 1.0_dp, V, na, fgr, nr, 0.0_dp, work, na)
                           CALL dgemm('T', 'N', 1, 1, na, 1.0_dp, work, na, fga, na, 1.0_dp, &
                                      intso(starti+iso, startj+jso, i), 1)

                           ! add the last term: r2*a2(j) * V * r2*a2(k)
                           fgr(1:nr, 1) = r2(1:nr, starti+iso)*r2(1:nr, startj+jso)
                           fga(1:na, 1) = a2(1:na, starti+iso, j)*a2(1:na, startj+jso, k)

                           CALL dgemm('N', 'N', na, 1, nr, 1.0_dp, V, na, fgr, nr, 0.0_dp, work, na)
                           CALL dgemm('T', 'N', 1, 1, na, 1.0_dp, work, na, fga, na, 1.0_dp, &
                                      intso(starti+iso, startj+jso, i), 1)

                        END DO !jso
                     END DO !iso

                  END DO !i
               END DO !jpgf
            END DO !ipgf
         END DO !jset
      END DO !iset

      DO i = 1, 3
         intso(:, :, i) = intso(:, :, i)-TRANSPOSE(intso(:, :, i))
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
      REAL(dp)                                           :: cspeed, zeff
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
      cspeed = 1.0_dp/a_fine

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
         CALL calculate_model_potential(Vr, grid, zeff)

         !Compute V/(4c^2-2V) and weight it
         ALLOCATE (V(na, nr))
         V = 0.0_dp
         DO ir = 1, nr
            CALL daxpy(na, Vr(ir)/(4.0_dp*cspeed**2-2.0_dp*Vr(ir)), grid%weight(1:na, ir), 1, &
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
