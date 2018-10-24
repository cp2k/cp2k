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
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: gto_basis_set_type, gto_basis_set_p_type, &
                                              get_gto_basis_set, set_gto_basis_set
   USE cell_types,                      ONLY: cell_type
   USE cp_array_utils,                  ONLY: cp_2d_r_p_type, cp_3d_r_p_type
   USE cp_control_types,                ONLY: qs_control_type, dft_control_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type, dbcsr_distribution_type,&
                                              dbcsr_create, dbcsr_print, dbcsr_put_block, &
                                              dbcsr_finalize, dbcsr_release, dbcsr_get_block_p, &
                                              dbcsr_release_p, dbcsr_iterator_type, &  
                                              dbcsr_iterator_blocks_left, dbcsr_iterator_start, &
                                              dbcsr_iterator_stop, dbcsr_iterator_next_block, &
                                              dbcsr_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE generic_os_integrals,            ONLY: int_overlap_ab_os
   USE input_section_types,             ONLY: section_vals_val_get
   USE kinds,                           ONLY: dp, default_string_length
   USE lebedev,                         ONLY: deallocate_lebedev_grids, get_number_of_lebedev_grid,&
                                              init_lebedev_grids, lebedev_grid
   USE mathconstants,                   ONLY: dfac, pi, gamma1
   USE mathlib,                         ONLY: invmat_symm
   USE memory_utilities,                ONLY: reallocate
   USE molecule_types,                  ONLY: molecule_type
   USE orbital_pointers,                ONLY: nsoset, indso, ncoset, nco, nso, indco
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: a_fine
   USE qs_environment_types,            ONLY: qs_environment_type, get_qs_env
   USE qs_grid_atom,                    ONLY: grid_atom_type, allocate_grid_atom, &
                                              create_grid_atom
   USE qs_harmonics_atom,               ONLY: create_harmonics_atom, get_maxl_CG, &
                                              harmonics_atom_type, allocate_harmonics_atom, &
                                              get_none0_cg_list
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: qs_kind_type, get_qs_kind, get_qs_kind_set
   USE qs_neighbor_lists,               ONLY: build_neighbor_lists, local_atoms_type, &
                                              pair_radius_setup, atom2d_cleanup, atom2d_build
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type, &
                                              deallocate_neighbor_list_set
   USE qs_o3c_methods,                  ONLY: calculate_o3c_integrals, contract12_o3c
   USE qs_o3c_types,                    ONLY: init_o3c_container, o3c_container_type, &
                                              release_o3c_container, o3c_iterate, &
                                              o3c_iterator_create, o3c_iterator_release, &
                                              o3c_iterator_type, get_o3c_iterator_info
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_types,                    ONLY: qs_rho_type, qs_rho_get
   USE spherical_harmonics,             ONLY: clebsch_gordon, clebsch_gordon_deallocate, &
                                              clebsch_gordon_init
   USE xas_tdp_types,                   ONLY: xas_atom_env_type
   USE xc_functionals_utilities,        ONLY: set_util, calc_rs_pw

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "xas_tdp_atom"

   PUBLIC :: init_xas_atom_env, integrate_fxc_atoms, integrate_soc_atoms

CONTAINS

! **************************************************************************************************
!> \brief Initializes a xas_atom_env type given the qs_env
!> \param xas_atom_env the xas_atom_env to initialize
!> \param excited_atoms the indices of the excited atoms
!> \param excited_kinds the indices of the excited kinds
!> \param grid_info the raw input giving the grid dimensions
!> \param qs_env ...
!> \param do_fxc whether the xc-kernel will be computed on the grid
! **************************************************************************************************
   SUBROUTINE init_xas_atom_env(xas_atom_env, excited_atoms, excited_kinds, grid_info, qs_env, do_fxc)

      TYPE(xas_atom_env_type), POINTER                :: xas_atom_env
      TYPE(qs_environment_type), POINTER              :: qs_env
      INTEGER, DIMENSION(:), POINTER                  :: excited_atoms
      INTEGER, DIMENSION(:), POINTER                  :: excited_kinds
      CHARACTER(len=default_string_length), &
                     DIMENSION(:,:), POINTER          :: grid_info
      LOGICAL, INTENT(IN)                             :: do_fxc

      CHARACTER(len=*), PARAMETER :: routineN = 'init_xas_atom_env', routineP = moduleN//":"//routineN

      INTEGER                                         :: nkind, natom, nex_kinds, handle, &
                                                         nex_atoms, ikind
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
                                             POINTER  :: basis_set_ri 
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
   
      NULLIFY(qs_kind_set, basis_set_ri, particle_set) 

      CALL timeset(routineN, handle)

!  Initializing the type
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, natom=natom, particle_set=particle_set)

      nkind = SIZE(qs_kind_set)   
      nex_kinds = SIZE(excited_kinds)
      nex_atoms = SIZE(excited_atoms)

      ALLOCATE(xas_atom_env%excited_atoms(nex_atoms))
      ALLOCATE(xas_atom_env%excited_kinds(nex_kinds))
      ALLOCATE(xas_atom_env%grid_atom_set(nkind))
      ALLOCATE(xas_atom_env%harmonics_atom_set(nkind))
      ALLOCATE(xas_atom_env%ri_sphi_so(nkind))
      ALLOCATE(xas_atom_env%orb_sphi_so(nkind))
      ALLOCATE(xas_atom_env%ri_blk_size(natom))
      ALLOCATE(xas_atom_env%ri_dcoeff(natom))

      xas_atom_env%excited_atoms = excited_atoms
      xas_atom_env%excited_kinds = excited_kinds

!  Allocate and initialize the atomic grids and harmonics
      CALL init_xas_atom_grid_harmo(xas_atom_env, grid_info, qs_env, do_fxc)

!  Create the block size for the integral matrix based on the size of the RI_XAS basis.
      ALLOCATE(basis_set_ri(nkind))
      CALL basis_set_list_setup(basis_set_ri, "RI_XAS", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=xas_atom_env%ri_blk_size, basis=basis_set_ri)
      DEALLOCATE(basis_set_ri)

!  Compute the contraction coefficients for spherical orbitals
      DO ikind = 1,nkind 
         CALL compute_sphi_so(ikind, "ORB", xas_atom_env%orb_sphi_so(ikind)%array, qs_env)
         IF (ANY(excited_kinds == ikind) .AND. do_fxc) THEN
            CALL compute_sphi_so(ikind, "RI_XAS", xas_atom_env%ri_sphi_so(ikind)%array, qs_env)
         END IF
      END DO !ikind

!  Compute the coefficients to expand the density in the RI_XAS basis, if requested
      IF (do_fxc) CALL calculate_density_coeffs(xas_atom_env, qs_env)

      CALL timestop(handle)

   END SUBROUTINE init_xas_atom_env

! **************************************************************************************************
!> \brief Initializes the atomic grids and harmonics for the RI atomic calculations
!> \param xas_atom_env ...
!> \param qs_env ...
!> \param do_fxc Whether the xc kernel will ne computed on the atomic grids. If not, the harmonics
!>        are built for the orbital basis for all kinds.
!> \note Largely inspired by init_rho_atom subroutine 
! **************************************************************************************************
   SUBROUTINE init_xas_atom_grid_harmo(xas_atom_env, grid_info, qs_env, do_fxc)

      TYPE(xas_atom_env_type), POINTER                :: xas_atom_env
      CHARACTER(len=default_string_length), &
                           DIMENSION(:,:), POINTER    :: grid_info
      TYPE(qs_environment_type), POINTER              :: qs_env
      LOGICAL, INTENT(IN)                             :: do_fxc

      CHARACTER(len=*), PARAMETER :: routineN = "init_xas_atom_grid_harmo", &
                                     routineP = moduleN//":"//routineN

      REAL(dp), DIMENSION(:, :, :), POINTER           :: my_CG
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)          :: rga
      INTEGER                                         :: llmax, maxlgto, max_s_harm, max_s_set, &
                                                         lcleb, lc1, iso1, l1, m1, lc2, iso2, &
                                                         l2, m2, l1l2, mp, mm, lp, il, iso, &
                                                         ikind, quadrature, ll, na, la, nr, &
                                                         maxs, maxl, igrid, stat
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(gto_basis_set_type), POINTER               :: tmp_basis
      TYPE(grid_atom_type), POINTER                   :: grid_atom
      TYPE(harmonics_atom_type), POINTER              :: harmonics
      TYPE(qs_control_type), POINTER                  :: qs_control
      TYPE(dft_control_type), POINTER                 :: dft_control
      CHARACTER(len=2)                                :: symbol

      NULLIFY(my_CG, qs_kind_set, tmp_basis, grid_atom, harmonics, qs_control, dft_control)

!  Initialization of some integer for the CG coeff generation
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, dft_control=dft_control)
      IF(do_fxc) THEN 
         CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto, basis_type="RI_XAS")
      ELSE
         CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto, basis_type="ORB")
      END IF
      qs_control => dft_control%qs_control

      !maximum expansion
      llmax = 2*maxlgto
      max_s_harm = nsoset(llmax)
      max_s_set = nsoset(maxlgto)
      lcleb = MAX(llmax, 2*maxlgto, 1)

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
      
      DO ikind = 1,SIZE(xas_atom_env%grid_atom_set)

!        Allocate the grid and the harmonics for this kind
         CALL allocate_grid_atom(xas_atom_env%grid_atom_set(ikind)%grid_atom)
         CALL allocate_harmonics_atom(xas_atom_env%harmonics_atom_set(ikind)%harmonics_atom)

         NULLIFY(grid_atom, harmonics)
         grid_atom => xas_atom_env%grid_atom_set(ikind)%grid_atom         
         harmonics => xas_atom_env%harmonics_atom_set(ikind)%harmonics_atom 

!        Initialize some integers
         CALL get_qs_kind(qs_kind_set(ikind), ngrid_rad=nr, ngrid_ang=na, element_symbol=symbol)

         !take the grid dimension given as input, if none, take the GAPW ones above
         IF (ANY(grid_info == symbol)) THEN
            DO igrid = 1, SIZE(grid_info,1)
               IF (grid_info(igrid,1) == symbol) THEN

                  !hack to convert string into integer
                  READ(grid_info(igrid,2), *, iostat=stat) na
                  IF (stat .NE. 0) CPABORT("The 'na' value for the GRID keyword must be an integer")
                  READ(grid_info(igrid,3), *, iostat=stat) nr
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
         IF (ANY(xas_atom_env%excited_atoms == ikind) .AND. do_fxc) THEN
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=tmp_basis, basis_type="RI_XAS")
         ELSE
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=tmp_basis, basis_type="ORB")
         END IF
         CALL get_gto_basis_set(gto_basis_set=tmp_basis, maxl=maxl)

         CALL create_grid_atom(grid_atom, nr, na, llmax, ll, quadrature)

         maxs = nsoset(maxl)
         CALL create_harmonics_atom(harmonics, &
                                    my_CG, na, llmax, maxs, max_s_harm, ll, grid_atom%wa, &
                                    grid_atom%azi, grid_atom%pol) 
         CALL get_maxl_CG(harmonics, tmp_basis, llmax, max_s_harm)

      END DO

      CALL deallocate_lebedev_grids()
      DEALLOCATE(my_CG)

   END SUBROUTINE init_xas_atom_grid_harmo
   
! **************************************************************************************************
!> \brief Computes the contraction coefficients to go from spherical orbitals to sgf for a given 
!>        atomic kind and a given basis type. 
!> \param ikind the kind for which we compute the coefficients
!> \param basis_type the type of basis for which we compute
!> \param sphi_so where the new contraction coefficients are stored (not yet allocated)
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_sphi_so(ikind, basis_type, sphi_so, qs_env)

      INTEGER, INTENT(IN)                             :: ikind
      CHARACTER(len=*), INTENT(IN)                    :: basis_type
      REAL(dp), DIMENSION(:,:), POINTER               :: sphi_so
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "compute_sphi_so", routineP = moduleN//":"//routineN

      INTEGER                                         :: nset, maxso, iset, sgfi, ipgf, start_s, &
                                                         start_c, l, iso, ico, lx, ly, lz
      REAL(dp)                                        :: factor
      TYPE(gto_basis_set_type), POINTER               :: basis 
      INTEGER, DIMENSION(:), POINTER                  :: lmax, npgf, nsgf_set, lmin
      INTEGER, DIMENSION(:,:), POINTER                :: first_sgf
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      REAL(dp), DIMENSION(:,:), POINTER               :: sphi

      NULLIFY(basis, lmax, lmin, npgf, nsgf_set, qs_kind_set, first_sgf, sphi)

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis, basis_type=basis_type)
      CALL get_gto_basis_set(basis, lmax=lmax, nset=nset, npgf=npgf, maxso=maxso, lmin=lmin, &
                             nsgf_set=nsgf_set, sphi=sphi, first_sgf=first_sgf)

      ALLOCATE(sphi_so(maxso,SUM(nsgf_set)))
      sphi_so = 0.0_dp

      DO iset = 1, nset
            sgfi = first_sgf(1, iset)
            DO ipgf = 1,npgf(iset)
               start_s = (ipgf-1)*nsoset(lmax(iset))
               start_c = (ipgf-1)*ncoset(lmax(iset))

               DO l = lmin(iset), lmax(iset)
                  DO iso = 1, nso(l)
                     DO ico = 1, nco(l)

                        lx = indco(1, ico+ncoset(l-1))
                        ly = indco(2, ico+ncoset(l-1))
                        lz = indco(3, ico+ncoset(l-1))

                        factor = orbtramat(l)%s2c(iso,ico) &
                               * SQRT(4.0_dp*pi/dfac(2*l+1)*dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1))

                        CALL daxpy(nsgf_set(iset), factor, &
                                   sphi(start_c+ncoset(l-1)+ico,sgfi:sgfi+nsgf_set(iset)-1), 1, &
                                   sphi_so(start_s+nsoset(l-1)+iso, sgfi:sgfi+nsgf_set(iset)-1), 1)

                     END DO !ico
                  END DO !iso
               END DO !l
            END DO !ipgf
         END DO !iset

   END SUBROUTINE compute_sphi_so


! **************************************************************************************************
!> \brief Builds a neighbor lists set taylored for 3-center integral within XAS TDP, such that only
!>        excited atoms are taken into account for the list_b
!> \param ac_list the neighbor list ready for 3-center integrals
!> \param basis_a basis set list for atom a
!> \param basis_c basis set list for atom c
!> \param op_type to indicate whther the list should be built with overlap or Coulomb in mind
!> \param excited_atoms the indices of the excited atoms to consider
!> \param excited_kinds the indices of the excited kinds to consider
!> \param qs_env ...
!> \note Based on setup_neighbor_list with added features
! **************************************************************************************************
   SUBROUTINE build_xas_tdp_3c_nl(ac_list, basis_a, basis_c, op_type, excited_atoms, &
                                  excited_kinds, qs_env)

      TYPE(neighbor_list_set_p_type), DIMENSION(:), &                                                
                                             POINTER  :: ac_list
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
                                             POINTER  :: basis_a, basis_c
      CHARACTER(len=*), INTENT(IN)                    :: op_type
      INTEGER, DIMENSION(:), POINTER                  :: excited_atoms, excited_kinds
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN =  "build_xas_tdp_3c_nl", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: nkind, ikind
      REAL(dp)                                        :: subcells
      LOGICAL, DIMENSION(:), ALLOCATABLE              :: a_present, c_present
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: a_radius, c_radius
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER   :: atomic_kind_set
      TYPE(cell_type), POINTER                        :: cell
      TYPE(distribution_1d_type), POINTER             :: distribution_1d
      TYPE(distribution_2d_type), POINTER             :: distribution_2d
      TYPE(local_atoms_type), DIMENSION(:), &
                                          ALLOCATABLE :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER      :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set

      NULLIFY(atomic_kind_set, distribution_1d, distribution_2d, molecule_set, particle_set, cell)
      
      IF (op_type == "COULOMB ") CPABORT("3-center Coulomb neighbor list NYI")

!  Initialization
      CALL get_qs_env(qs_env, nkind=nkind)
      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)

      ALLOCATE(a_present(nkind), c_present(nkind))
      a_present = .FALSE.
      c_present = .FALSE.
      ALLOCATE(a_radius(nkind), c_radius(nkind))
      a_radius = 0.0_dp
      c_radius = 0.0_dp

!  Set up the radii
      DO ikind = 1,nkind
         !orbital basis set
         IF (ASSOCIATED(basis_a(ikind)%gto_basis_set)) THEN
            a_present(ikind) = .TRUE.
            CALL get_gto_basis_set(basis_a(ikind)%gto_basis_set, kind_radius=a_radius(ikind))
         END IF
         !RI_XAS basis set, take only excited kinds
         IF (ASSOCIATED(basis_c(ikind)%gto_basis_set) .AND. ANY(excited_kinds == ikind)) THEN
            c_present(ikind) = .TRUE.
            CALL get_gto_basis_set(basis_c(ikind)%gto_basis_set, kind_radius=c_radius(ikind))
         END IF
      END DO !ikind

      ALLOCATE(pair_radius(nkind, nkind))
      pair_radius = 0.0_dp
      CALL pair_radius_setup(a_present, c_present, a_radius, c_radius, pair_radius)

!  Actually setup the list
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, cell=cell, &
                      distribution_2d=distribution_2d, local_particles=distribution_1d, &
                      particle_set=particle_set, molecule_set=molecule_set)

      ALLOCATE(atom2d(nkind))
      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                        molecule_set, .FALSE., particle_set)

      CALL build_neighbor_lists(ac_list, particle_set, atom2d, cell, pair_radius, subcells, &
                                operator_type = "ABC", atomb_to_keep=excited_atoms, &
                                nlname="XAS_TDP_3c_nl")

!  Clean-up
      CALL atom2d_cleanup(atom2d) 

   END SUBROUTINE build_xas_tdp_3c_nl

! **************************************************************************************************
!> \brief Compute the coefficients to project the density on the RI_XAS basis set
!> \param xas_atom_env ...
!> \param qs_env ...
!> \note The density is n = sum_ab P_ab*phi_a*phi_b, the RI basis covers the products of orbital sgfs
!>       => n = sum_ab sum_cd P_ab (phi_a phi_b xi_c) S_cd^-1 xi_d
!>            = sum_d coeff_d xi_d , where xi are the RI basis func.
!>       The 3-center overlap is needed, as well as the inverse of the RI basis overlap. Note that 
!>       the RI basis functions are only centered on atoms of excited kinds 
! **************************************************************************************************
   SUBROUTINE calculate_density_coeffs(xas_atom_env, qs_env)

      TYPE(xas_atom_env_type), POINTER                :: xas_atom_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "calculate_density_coeffs", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: nkind, ilist, katom, kkind, ikind, natom, &
                                                         iat, nsgf_kind
      REAL(dp)                                        :: dmax
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(gto_basis_set_p_type), DIMENSION(:), &
                                             POINTER  :: basis_set_ri, basis_set_orb
      TYPE(gto_basis_set_type), POINTER               :: ri_basis
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &                                                
                                             POINTER  :: ac_list, sab_orb
      TYPE(o3c_container_type), POINTER               :: o3c
      TYPE(qs_rho_type), POINTER                      :: rho
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: rho_ao
      TYPE(o3c_iterator_type)                         :: o3c_iterator
      REAL(dp), DIMENSION(:,:), POINTER               :: tvec
      TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER     :: sinv
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set 
      INTEGER, DIMENSION(:), POINTER                  :: nsgf
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: work

      NULLIFY(qs_kind_set, basis_set_ri, basis_set_orb, ac_list, rho, rho_ao, sinv, ri_basis)
      NULLIFY(o3c, tvec, particle_set, nsgf)
      CPASSERT(ASSOCIATED(xas_atom_env))

      CALL get_qs_env(qs_env, nkind=nkind, sab_orb=sab_orb, qs_kind_set=qs_kind_set, rho=rho, &
                      natom=natom, particle_set=particle_set)

!  TODO: make sure that incomplete RI basis are delt with correctly
!  Create the needed neighbor list and basis set lists.
      ALLOCATE(basis_set_ri(nkind))
      ALLOCATE(basis_set_orb(nkind))
      CALL basis_set_list_setup(basis_set_ri, "RI_XAS", qs_kind_set)
      CALL basis_set_list_setup(basis_set_orb, "ORB", qs_kind_set)
      CALL build_xas_tdp_3c_nl(ac_list, basis_set_orb, basis_set_ri, "OVERLAP", &
                               xas_atom_env%excited_atoms, xas_atom_env%excited_kinds, qs_env)

!  Create the inverse RI overlap matrices. 
      ALLOCATE(sinv(nkind))
      DO ikind = 1,nkind
         NULLIFY(sinv(ikind)%array)
         IF (.NOT. ANY(xas_atom_env%excited_kinds == ikind)) CYCLE
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=ri_basis, basis_type="RI_XAS",nsgf=nsgf_kind)
         ALLOCATE(sinv(ikind)%array(nsgf_kind,nsgf_kind))
         ALLOCATE(work(nsgf_kind, nsgf_kind))

         CALL int_overlap_ab_os(sab=work, rab=(/0.0_dp,0.0_dp,0.0_dp/), fba=ri_basis, &
                                fbb=ri_basis, calculate_forces=.FALSE., debug=.FALSE.,  dmax=dmax)
         CALL invmat_symm(work)
         sinv(ikind)%array(:,:) = work(:,:)
         DEALLOCATE(work)
      END DO

!  Create the 3-center overlap integrals for the RI basis
      ALLOCATE(o3c)
      CALL init_o3c_container(o3c, 1, basis_set_orb, basis_set_orb, basis_set_ri, sab_orb, ac_list)
      CALL calculate_o3c_integrals(o3c)

!  Contract the 3-center tensor with the density matrix
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CALL contract12_o3c(o3c, rho_ao)

!  Allocate space for the projected density coefficients. 
      nsgf => xas_atom_env%ri_blk_size
      DO iat = 1,natom
         IF (.NOT. ANY(xas_atom_env%excited_atoms == iat)) CYCLE
         ikind = particle_set(iat)%atomic_kind%kind_number
         ALLOCATE(xas_atom_env%ri_dcoeff(iat)%array(nsgf(iat)))
         xas_atom_env%ri_dcoeff(iat)%array = 0.0_dp
      END DO

!  Iterate over the o3c to access the contracted tvec = sum_ab P_ab*(abc). 
!  TODO: careful when parallelizing
      CALL o3c_iterator_create(o3c, o3c_iterator) 
      DO WHILE (o3c_iterate(o3c_iterator) == 0)
         CALL get_o3c_iterator_info(o3c_iterator, katom=katom, kkind=kkind, tvec=tvec)
         IF (.NOT. ANY(xas_atom_env%excited_atoms == katom)) CYCLE

!        The coefficient for atom katom of kkind are tvec*S^-1
         ALLOCATE(work(nsgf(katom),1))
         CALL dgemm('N', 'N', nsgf(katom), 1, nsgf(katom), 1.0_dp, sinv(kkind)%array, nsgf(katom), &
                    tvec, nsgf(katom), 0.0_dp, work, nsgf(katom))
         xas_atom_env%ri_dcoeff(katom)%array(:) = xas_atom_env%ri_dcoeff(katom)%array(:) + work(:,1)
         DEALLOCATE(work)

      END DO
      CALL o3c_iterator_release(o3c_iterator)

!  some clean-up
      DO ilist = 1,SIZE(ac_list)
         CALL deallocate_neighbor_list_set(ac_list(ilist)%neighbor_list_set)
      END DO
      DO ikind = 1,nkind
         IF (ASSOCIATED(sinv(ikind)%array)) DEALLOCATE(sinv(ikind)%array)
      END DO
      DEALLOCATE(sinv)
      CALL release_o3c_container(o3c)
      DEALLOCATE(basis_set_ri, basis_set_orb, o3c, ac_list)

   END SUBROUTINE calculate_density_coeffs

! **************************************************************************************************
!> \brief Evaluates the density on a given atomic grid
!> \param rho the array containing the projectd density at each grid point
!> \param atom_index the index of the atom in question
!> \param atom_kind the kind of the atom in question
!> \param xas_atom_env ...
!> \param qs_env ...
!> \note The density is expressed as n = sum_d coeff_d*xi_d. Knowing the coordinate of each grid
!>       grid point, one can simply evaluate xi_d(r) 
! **************************************************************************************************
   SUBROUTINE put_density_on_atomic_grid(rho, atom_index, atom_kind, xas_atom_env, qs_env)

      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: rho
      INTEGER, INTENT(IN)                             :: atom_index, atom_kind
      TYPE(xas_atom_env_type), POINTER                :: xas_atom_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "put_density_on_atomic_grid", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: na, nr, iset, ipgf, l, iso, ia, nset, &
                                                         isgf, start, sgfi, n
      TYPE(grid_atom_type), POINTER                   :: grid_atom
      TYPE(harmonics_atom_type), POINTER              :: harmonics
      TYPE(gto_basis_set_type), POINTER               :: ri_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      INTEGER, DIMENSION(:), POINTER                  :: lmax, npgf, nsgf_set, lmin
      REAL(dp), DIMENSION(:,:), POINTER               :: zet, ri_sphi_so
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE         :: so, sgf
      INTEGER, DIMENSION(:,:), POINTER                :: first_sgf

      NULLIFY(grid_atom, ri_basis, qs_kind_set, harmonics, lmax, npgf, zet, nsgf_set, ri_sphi_so)
      NULLIFY(lmin, first_sgf)

!  Strategy: it makes sense to evaluate the spherical orbital on the grid (because of symmetry)
!            From there, one can directly contract into sgf using ri_sphi_so and then take the weight

!  Generalities
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(atom_kind), basis_set=ri_basis, basis_type="RI_XAS")
      CALL get_gto_basis_set(ri_basis, lmax=lmax, npgf=npgf, zet=zet, nset=nset, nsgf_set=nsgf_set,&
                             first_sgf=first_sgf, lmin=lmin)

!  Get the grid, the harmonics and the info we need from it
      grid_atom => xas_atom_env%grid_atom_set(atom_kind)%grid_atom
      harmonics => xas_atom_env%harmonics_atom_set(atom_kind)%harmonics_atom
      CPASSERT(ASSOCIATED(grid_atom))
      CPASSERT(ASSOCIATED(harmonics))
      nr = grid_atom%nr
      na = grid_atom%ng_sphere
      n = na*nr

!  Loop over the sets and pgfs of the basis set
      DO iset = 1,nset

         !allocate space to store the spherical orbitals on the grid
         ALLOCATE(so(npgf(iset)*nsoset(lmax(iset)),na,nr))
         so = 0.0_dp

         DO ipgf = 1,npgf(iset)
            start = (ipgf-1)*nsoset(lmax(iset))

            !loop over the spherical gaussian orbitals 
            DO iso = nsoset(lmin(iset)-1)+1, nsoset(lmax(iset))
               l = indso(1, iso)
               DO ia = 1,na
                  !compute: g = r**l * Y_lm * exp(-zet* r**2)
                  so(start+iso,ia,1:nr) = EXP(-zet(ipgf, iset)*grid_atom%rad2(1:nr)) & 
                                            *grid_atom%rad(1:nr)**l *harmonics%slm(ia,iso)
               END DO !ia
            END DO !iso

         END DO !ipgf

         !contract the so to get the sgf
         ALLOCATE(sgf(nsgf_set(iset),na,nr))
         sgf = 0.0_dp
         sgfi = first_sgf(1,iset) - 1

         DO isgf = 1,nsgf_set(iset)

            ri_sphi_so => xas_atom_env%ri_sphi_so(atom_kind)%array
            DO iso = 1,npgf(iset)*nsoset(lmax(iset))
               CALL daxpy(n, ri_sphi_so(iso,sgfi+isgf), so(iso,:,:), 1, sgf(isgf,:,:), 1)
            END DO !iso

            !put the sgf on the grid with the approriate coefficients and sum
            CALL daxpy(n, xas_atom_env%ri_dcoeff(atom_index)%array(sgfi+isgf), sgf(isgf,:,:), &
                       1, rho(:,:), 1)

         END DO !isgf

         DEALLOCATE(so, sgf)
      END DO !iset

   END SUBROUTINE put_density_on_atomic_grid

! **************************************************************************************************
!> \brief Integrate the xc kernel as a function of r on the atomic grids for the RI_XAS basis
!> \param int_mat the dbcsr_matrix in the usual block format where the integrals are stored
!> \param xas_atom_env ...
!> \param qs_env ...
!> \param singlet whether we want singlet or triplet excitations
! **************************************************************************************************
   SUBROUTINE integrate_fxc_atoms(int_mat, xas_atom_env, qs_env, singlet)

      TYPE(dbcsr_type), POINTER                       :: int_mat
      TYPE(xas_atom_env_type), POINTER                :: xas_atom_env  
      TYPE(qs_environment_type), POINTER              :: qs_env   
      LOGICAL, INTENT(IN)                             :: singlet

      CHARACTER(len=*), PARAMETER :: routineN = "integrate_fxc_atoms", &
                                     routineP = moduleN//":"//routineN

      TYPE(dbcsr_distribution_type), POINTER          :: dbcsr_dist
      REAL(dp), DIMENSION(:, :), ALLOCATABLE          :: fxc, int_so, int_sgf, rho
      INTEGER                                         :: natom, iat, ikind, na, nr, nsotot, nset, &
                                                         maxso, iset, jset, nsgf, sgfi, jat, blk, &
                                                         sgfj, n, nsoi, nsoj, starti, startj
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set 
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      TYPE(grid_atom_type), POINTER                   :: grid_atom 
      TYPE(harmonics_atom_type), POINTER              :: harmonics
      TYPE(gto_basis_set_type), POINTER               :: ri_basis 
      REAL(dp), DIMENSION(:,:), POINTER               :: ri_sphi_so
      INTEGER, DIMENSION(:), POINTER                  :: nsgf_set, npgf, lmax
      INTEGER, DIMENSION(:,:), POINTER                :: first_sgf
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: rs
      TYPE(dft_control_type), POINTER                 :: dft_control
      TYPE(dbcsr_iterator_type)                       :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s

      NULLIFY(dbcsr_dist, particle_set, qs_kind_set, grid_atom, harmonics, dft_control)
      NULLIFY(ri_basis, nsgf_set, first_sgf, lmax, npgf, ri_sphi_so, matrix_s)

!  Create the integral matrices based on the usual dbcsr_dist and the custom xas_atom_env block sizes
      CALL get_qs_env(qs_env, dbcsr_dist=dbcsr_dist, matrix_s=matrix_s)
      CALL dbcsr_create(matrix=int_mat, name="RI atomic Fxc integrals", matrix_type='S', &         
                        dist=dbcsr_dist, row_blk_size=xas_atom_env%ri_blk_size, &                        
                        col_blk_size=xas_atom_env%ri_blk_size)

      CALL get_qs_env(qs_env, particle_set=particle_set, qs_kind_set=qs_kind_set, natom=natom, &
                      dft_control=dft_control)

!  Loop over the atoms of interest. (iterate over diagonal of matrix_s cuz same distribution and
!                                    filled diagonal blocks) 
      CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix) 
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iat, column=jat, blk=blk)

         IF (.NOT. iat == jat) CYCLE
         IF (.NOT. ANY(xas_atom_env%excited_atoms == iat)) CYCLE
         ikind = particle_set(iat)%atomic_kind%kind_number

         grid_atom => xas_atom_env%grid_atom_set(ikind)%grid_atom
         harmonics => xas_atom_env%harmonics_atom_set(ikind)%harmonics_atom
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=ri_basis, basis_type="RI_XAS")

!     Allocate the density, kernel and integral arrays
         na = grid_atom%ng_sphere
         nr = grid_atom%nr
         CALL get_gto_basis_set(ri_basis, nset=nset, maxso=maxso)
         nsotot = nset*maxso 

         ALLOCATE(rho(na, nr)) 
         ALLOCATE(fxc(na, nr))
         ALLOCATE(int_so(nsotot, nsotot))
         int_so = 0.0_dp

!     Put the density on the grid
         CALL put_density_on_atomic_grid(rho, iat, ikind, xas_atom_env, qs_env)

!     Compute the fxc kernel on the atomic grids. Use calc_rs_pw because it takes 3d arrays as arg
         CALL set_util(dft_control%qs_control%eps_rho_rspace)
         n = na*nr
         ALLOCATE(rs(n))

         CALL calc_rs_pw(rho, rs, n)
         fxc = 0.0_dp

         IF (singlet) THEN
            CALL pade_kernel_singlet(n, rho, rs, fxc, dft_control%qs_control%eps_rho_rspace)
         ELSE
            CALL pade_kernel_triplet(n, rho, rs, fxc, dft_control%qs_control%eps_rho_rspace)
         END IF

!     Take the grid weights into account
         fxc(1:na,1:nr) = fxc(1:na,1:nr)*grid_atom%weight(1:na,1:nr)

!     Integrate
         CALL integrate_so_prod(int_so, fxc, ikind, xas_atom_env, qs_env) 

         DEALLOCATE(rs, fxc, rho)

!     Contract the integrals into sgf using the ri_sphi_so arrays
         CALL get_gto_basis_set(ri_basis, nsgf=nsgf, nsgf_set=nsgf_set, first_sgf=first_sgf, &
                                npgf=npgf, lmax=lmax)

         ALLOCATE(int_sgf(nsgf,nsgf))
         int_sgf = 0.0_dp
         ri_sphi_so => xas_atom_env%ri_sphi_so(ikind)%array

         DO iset = 1, nset
            starti = (iset-1)*maxso+1
            nsoi = npgf(iset)*nsoset(lmax(iset))
            sgfi = first_sgf(1, iset)

            DO jset= 1, nset
               startj = (jset-1)*maxso+1
               nsoj = npgf(jset)*nsoset(lmax(jset))
               sgfj = first_sgf(1, jset)

               CALL ab_contract(int_sgf(sgfi:sgfi+nsgf_set(iset)-1, sgfj:sgfj+nsgf_set(jset)-1), & 
                                int_so(starti:starti+nsoi-1, startj:startj+nsoj-1), &
                                ri_sphi_so(:,sgfi:), ri_sphi_so(:,sgfj:), nsoi, nsoj, &
                                nsgf_set(iset), nsgf_set(jset))
            END DO !jset
         END DO !iset

!     Put the integrals for this atom as blocks in the dbcsr matrices
         CALL dbcsr_put_block(int_mat, iat, iat, int_sgf)

!     Clean-up 
         DEALLOCATE(int_sgf, int_so)

      END DO !iat
      CALL dbcsr_iterator_stop(iter)

!  Finalize the matrix
      CALL dbcsr_finalize(int_mat)

   END SUBROUTINE integrate_fxc_atoms


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

      REAL(dp), DIMENSION(:,:), INTENT(INOUT)         :: intso
      REAL(dp), DIMENSION(:,:), INTENT(IN)            :: fxc
      INTEGER, INTENT(IN)                             :: ikind
      TYPE(xas_atom_env_type), POINTER                :: xas_atom_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "integrate_so_prod", &
                                     routineP = moduleN//":"//routineN 

      INTEGER                                         :: handle, nset, maxso, maxl, nr, na, &
                                                         max_iso_not0, max_s_harm, ld, ia, nngau1, &
                                                         m1, iset1, iset2, n1, m2, n2, ngau1, ngau2,&
                                                         max_iso_not0_local, size1, ipgf1, ipgf2, &
                                                         lmin12, lmax12, iso, icg, iso1, iso2, l, ic
      TYPE(grid_atom_type), POINTER                   :: grid_atom
      TYPE(harmonics_atom_type), POINTER              :: harmonics
      TYPE(gto_basis_set_type), POINTER               :: ri_basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      INTEGER, DIMENSION(:), POINTER                  :: lmax, lmin, npgf
      REAL(dp), DIMENSION(:,:), POINTER               :: zet
      REAL(dp), DIMENSION(:,:,:), POINTER             :: my_CG
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: g1, g2
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: gg, gfxcg, matso
      INTEGER, DIMENSION(:), ALLOCATABLE              :: cg_n_list
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE          :: cg_list
      
      CALL timeset(routineN, handle)

      NULLIFY(grid_atom, harmonics, ri_basis, qs_kind_set, lmax, lmin, npgf, zet, my_CG)

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

      ALLOCATE(g1(nr), g2(nr), gg(nr, 0:2*maxl))
      ALLOCATE(gfxcg(na, 0:2*maxl))
      ALLOCATE(matso(nsoset(maxl), nsoset(maxl)))
      ALLOCATE(cg_list(2, nsoset(maxl)**2, max_s_harm), cg_n_list(max_s_harm))

      g1 = 0.0_dp
      g2 = 0.0_dp
      m1 = 0
!  Loop over the product of so
      DO iset1 = 1,nset
         n1 = nsoset(lmax(iset1))
         m2 = 0
         DO iset2 = 1,nset
            CALL get_none0_cg_list(my_CG, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                   max_s_harm, lmax(iset1)+lmax(iset2), cg_list, cg_n_list, &
                                   max_iso_not0_local)
            CPASSERT(max_iso_not0_local .LE. max_iso_not0)

            n2 = nsoset(lmax(iset2))
            DO ipgf1 = 1,npgf(iset1)
               ngau1 = n1*(ipgf1-1)+m1
               size1 = nsoset(lmax(iset1))-nsoset(lmin(iset1)-1)
               nngau1 = nsoset(lmin(iset1)-1)+ngau1

               g1(1:nr) = EXP(-zet(ipgf1, iset1)*grid_atom%rad2(1:nr))
               DO ipgf2 = 1,npgf(iset2)
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
                           matso(iso1, iso2) = matso(iso1, iso2) + gfxcg(ia,l) *&
                                               my_CG(iso1, iso2, iso) * harmonics%slm(ia, iso)
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
!> \brief The actual pade kernel evaluation in the singlet case. Inspired from xc_pade.F
!> \param n the number of grid points
!> \param rho the array of the values of the density on the grid
!> \param rs the array of the values of rs on the grid
!> \param fxc the array of the values of fxc on the grid
!> \param eps_rho the density cutoff
!> \note TODO: OMP parallelize as in xc_pade.F
!>             In the singlet case, take d2exc_d2na + d2exc_dnadnb, can be shown to be the formula
!>             below, only dependent on n (spin-unpolarized)
! **************************************************************************************************
   SUBROUTINE pade_kernel_singlet(n, rho, rs, fxc, eps_rho)

      REAL(dp), DIMENSION(*), INTENT(IN)              :: rho, rs
      REAL(dp), DIMENSION(*), INTENT(INOUT)           :: fxc
      REAL(dp), INTENT(IN)                            :: eps_rho
      INTEGER                                         :: n

      CHARACTER(len=*), PARAMETER :: routineN = "pade_kernel_singlet", routineP=moduleN//":"//routineN 

      REAL(dp), PARAMETER  :: f13 = 1.0_dp/3.0_dp, &
                              f49 = 4.0_dp/9.0_dp, &
                              two = 2.0_dp, &
                              three = 3.0_dp, &
                              four = 4.0_dp, &
                              six = 6.0_dp, &
                              twelve = 12.0_dp

      REAL(dp), PARAMETER  :: a0 = 0.4581652932831429E+0_dp, &                                      
                              a1 = 0.2217058676663745E+1_dp, &
                              a2 = 0.7405551735357053E+0_dp, & 
                              a3 = 0.1968227878617998E-1_dp, &
                              b1 = 1.0000000000000000E+0_dp, &
                              b2 = 0.4504130959426697E+1_dp, &
                              b3 = 0.1110667363742916E+1_dp, & 
                              b4 = 0.2359291751427506E-1_dp 

      REAL(dp)                                        :: drs, d2rs, p, dp, d2p, q, dq, d2q, df, d2f
      INTEGER                                         :: ip

!     the functional has the form:               e = -f(rs) = -p(rs)/q(rs)
!     the first derivative wrt to rho (n):       de_dn = -df_drs*drs_dn
!     the seconde derivative wrt rho :           d2e_d2n = -d2f_d2rs*drs_dn**2 - df_drs*d2rs_d2n

!     rs has the form:                           rs = (3/4*pi)**1/3 * n**(-1/3)
!     the first derivative wrt n:                drs = -1/3 * rs * n**(-1)
!     the second derivative wrt n:               d2rs = 4/9 * rs * n**(-2)

!     We need p,q and the derivatives:           p = a0 + a1*rs + a2*rs**2 + a3*rs**3  
!                                                dp = a1 + 2*a2*rs + 3*a3*rs**2                     
!                                                d2p = 2*a2 + 6*a3*rs

!                                                q = b1*rs + b2*rs**2 + b3*rs**3 + b4*rs**4
!                                                dq = b1 + 2*b2*rs + 3*b3*rs**2 + 4*b4*rs**3
!                                                d2q = 2*b2 + 6*b3*rs + 12*b4*rs**2

!     Finally the derivatives of f wrt rs:       df = (dp*q - p*dq)/q**2
!                                                d2f = -2*dp*dq/q**2 + d2p/q 
!                                                      + p*( 2*dq**2/q**3 - d2q/q**2)

      DO ip = 1,n
         IF (rho(ip) > eps_rho) THEN

            p = a0+(a1+(a2+a3*rs(ip))*rs(ip))*rs(ip)
            q = (b1+(b2+(b3+b4*rs(ip))*rs(ip))*rs(ip))*rs(ip)  

            dp = a1 + two*a2*rs(ip) + three*a3*rs(ip)**2
            dq = b1 + two*b2*rs(ip) + three*b3*rs(ip)**2 + four*b4*rs(ip)**3

            d2p = two*a2 + six*a3*rs(ip)
            d2q = two*b2 + six*b3*rs(ip) + twelve*b4*rs(ip)**2

            df = (dp*q - p*dq)/q**2
            d2f = -two*dp*dq/q**2 + d2p/q + p*( two*dq**2/q**3 - d2q/q**2)

            drs = -f13*rs(ip)/rho(ip)
            d2rs = f49*rs(ip)/rho(ip)**2

            !the factor 2 comes from the sum of the 2 derivatives
            fxc(ip) = -two*(d2f*drs**2 + df*d2rs)

         END IF
      END DO

   END SUBROUTINE pade_kernel_singlet

! **************************************************************************************************
!> \brief The actual pade kernel evaluation in the triplet case. Inspired from xc_pade.F
!> \param n the number of grid points
!> \param rho the array of the values of the density on the grid                                           
!> \param rs the array of the values of rs on the grid                                                     
!> \param fxc the array of the values of fxc on the grid                                                   
!> \param eps_rho the density cutoff                                                                       
!> \note TODO: OMP parallelize as in xc_pade.F                                                       
!>             In the triplet case, take d2exc_d2na - d2exc_dnadnb, can be shown to be the formula   
!>             below, only dependent on n (spin-unpolarized) 
! **************************************************************************************************
   SUBROUTINE pade_kernel_triplet(n, rho, rs, fxc, eps_rho)
      
      REAL(dp), DIMENSION(*), INTENT(IN)              :: rho, rs  
      REAL(dp), DIMENSION(*), INTENT(INOUT)           :: fxc  
      REAL(dp), INTENT(IN)                            :: eps_rho
      INTEGER                                         :: n

      CHARACTER(len=*), PARAMETER :: routineN = "pade_kernel_triplet", routineP=moduleN//":"//routineN

      REAL(dp), PARAMETER  :: a0 = 0.4581652932831429E+0_dp, &
                              a1 = 0.2217058676663745E+1_dp, & 
                              a2 = 0.7405551735357053E+0_dp, &
                              a3 = 0.1968227878617998E-1_dp, &
                              b1 = 1.0000000000000000E+0_dp, &
                              b2 = 0.4504130959426697E+1_dp, &
                              b3 = 0.1110667363742916E+1_dp, &
                              b4 = 0.2359291751427506E-1_dp

      REAL(dp), PARAMETER  :: da0 = 0.119086804055547E+0_dp, &
                              da1 = 0.6157402568883345E+0_dp, &
                              da2 = 0.1574201515892867E+0_dp, & 
                              da3 = 0.3532336663397157E-2_dp, &
                              db1 = 0.0000000000000000E+0_dp, &
                              db2 = 0.2673612973836267E+0_dp, &
                              db3 = 0.2052004607777787E+0_dp, &
                              db4 = 0.4200005045691381E-2_dp 

      REAL(dp), PARAMETER  :: d2f_d2xi = 4.0_dp/9.0_dp * 1.0_dp/(2.0_dp**(1.0_dp/3.0_dp)-1.0_dp), &
                              two = 2.0_dp

      REAL(dp)                                        :: p, dp, q, dq, de_df
      INTEGER                                         :: ip

!     We can express the functional as a function of fx =>        e = e(f) = -p(f)/q(f)
!     In spin polarized cases, one can show that we only need:    de_df
!                                                         and:    d2f_d2xi

!     The first derivitive of the functional has the form:        de_df = - (dp*q - p*dq)/q^2

!     We need p,q and their derivative wrt to f (at xi=0):        p = a0 + a1*rs + a2*rs^2 +a3*rs^3
!                                                                 dp = da0 +da1*rs +da2*rs^2 +da3*rs^3 

!                                                                 q = b1*rs + b2*rs^2 + b3*rs^3 +b4*rs^4
!                                                                 dq = db1*rs +db2*rs^2 +db3*rs^3 +db4*rs^4

!     The double derivative of f wrt xi (at xi=0):                d2f_d2xi = 4/9 * 1/(2^1/3 - 1)

!     In the end, we compute:                                     2/n^2 * de_df * d2f_d2xi

      DO ip = 1,n
         IF (rho(ip) > eps_rho) THEN 

            p = a0+(a1+(a2+a3*rs(ip))*rs(ip))*rs(ip)
            q = (b1+(b2+(b3+b4*rs(ip))*rs(ip))*rs(ip))*rs(ip)   

            dp = da0+(da1+(da2+da3*rs(ip))*rs(ip))*rs(ip)
            dq = (db1+(db2+(db3+db4*rs(ip))*rs(ip))*rs(ip))*rs(ip)

            de_df = -(dp*q - p*dq)/q**2

            fxc(ip) = two/(rho(ip))**2 * de_df *d2f_d2xi

         END IF
      END DO

   END SUBROUTINE pade_kernel_triplet

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
      
      REAL(dp), DIMENSION(:,:,:), INTENT(INOUT)       :: intso
      REAL(dp), DIMENSION(:,:), INTENT(IN)            :: V
      INTEGER, INTENT(IN)                             :: ikind
      TYPE(xas_atom_env_type), POINTER                :: xas_atom_env
      TYPE(qs_environment_type), POINTER              :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = "integrate_so_dxdy_prod", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: na, nr, iso, jso, maxso, nset, iset,&
                                                         jset, ipgf, jpgf, l, starti, startj, i, j,k
      TYPE(grid_atom_type), POINTER                   :: grid_atom
      TYPE(harmonics_atom_type), POINTER              :: harmonics
      TYPE(gto_basis_set_type), POINTER               :: basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      REAL(dp), DIMENSION(:, :), POINTER              :: slm
      REAL(dp), DIMENSION(:,:,:), POINTER             :: dslm_dxyz
      INTEGER, DIMENSION(:), POINTER                  :: lmax, lmin, npgf
      REAL(dp), DIMENSION(:, :), POINTER              :: zet
      REAL(dp), ALLOCATABLE, DIMENSION(:,:)           :: r1,  work, fga, fgr, r2
      REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)         :: a1, a2

      NULLIFY(grid_atom, harmonics, basis, qs_kind_set, dslm_dxyz, slm, lmin, lmax, npgf, zet)

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
                             maxso=maxso, npgf=npgf, nset=nset, zet=zet )


!  Separate the functions into purely r and purely angular parts, compute them all
!  and use matrix mutliplication for the integral. We use f for x derivative ang g for y

   ! Separating the functions. Note that the radial part is the same for x and y derivatives
      ALLOCATE(a1(na,nset*maxso,3), a2(na,nset*maxso,3))
      ALLOCATE(r1(nr,nset*maxso), r2(nr,nset*maxso))
      a1 = 0.0_dp; a2 = 0.0_dp
      r1 = 0.0_dp; r2 = 0.0_dp

      DO iset = 1,nset
         DO ipgf = 1,npgf(iset)
            starti = (iset-1)*maxso + (ipgf-1)*nsoset(lmax(iset))
            DO iso = nsoset(lmin(iset)-1)+1,nsoset(lmax(iset))
               l = indso(1,iso)

               ! The x derivitive of the spherical orbital, divided in angular and radial parts
               ! Two of each are needed because d/dx(r^l Y_lm) * exp(-al*r^2) + r^l Y_lm * ! d/dx(exp-al*r^2)

               ! the purely radial part of d/dx(r^l Y_lm) * exp(-al*r^2) (same for y)
               r1(1:nr,starti+iso) = grid_atom%rad(1:nr)**(l-1)*EXP(-zet(ipgf,iset)*grid_atom%rad2(1:nr))

               ! the purely radial part of r^l Y_lm * d/dx(exp-al*r^2) (same for y)
               r2(1:nr,starti+iso) = -2.0_dp*zet(ipgf,iset)*grid_atom%rad(1:nr)**(l+1) &
                                      *EXP(-zet(ipgf, iset)*grid_atom%rad2(1:nr))

               DO i = 1,3
                  ! the purely angular part of d/dx(r^l Y_lm) * exp(-al*r^2)
                  a1(1:na,starti+iso,i) = dslm_dxyz(i,1:na, iso)

                  ! the purely angular part of r^l Y_lm * d/dx(exp-al*r^2)
                  a2(1:na,starti+iso,i) = harmonics%a(i,1:na)*slm(1:na, iso)
               END DO

            END DO !iso
         END DO !ipgf
      END DO !iset

   ! Do the integration in terms of so using matrix products
      intso = 0.0_dp
      ALLOCATE(fga(na,1))
      ALLOCATE(fgr(nr,1))
      ALLOCATE(work(na,1))
      fga = 0.0_dp; fgr = 0.0_dp; work = 0.0_dp

      DO iset = 1,nset
         DO jset = 1,nset
            DO ipgf = 1,npgf(iset)
               starti = (iset-1)*maxso + (ipgf-1)*nsoset(lmax(iset))
               DO jpgf = 1,npgf(jset)
                  startj = (jset-1)*maxso + (jpgf-1)*nsoset(lmax(jset))

                  DO i = 1,3
                     j = MOD(i,3)+1
                     k = MOD(i+1,3)+1

                     DO iso = nsoset(lmin(iset)-1)+1,nsoset(lmax(iset))
                        DO jso = nsoset(lmin(jset)-1)+1,nsoset(lmax(jset))

                           !Two component per function => 4 terms in total

                           ! take r1*a1(j) * V * r1*a1(k)
                           fgr(1:nr,1) = r1(1:nr,starti+iso)*r1(1:nr,startj+jso)
                           fga(1:na,1) = a1(1:na,starti+iso,j)*a1(1:na,startj+jso,k)

                           CALL dgemm('N', 'N', na, 1, nr, 1.0_dp, V, na, fgr, nr, 0.0_dp, work, na)
                           CALL dgemm('T', 'N', 1, 1, na, 1.0_dp, work, na, fga, na, 0.0_dp, &
                                      intso(starti+iso, startj+jso, i), 1)

                           ! add r1*a1(j) * V * r2*a2(k)
                           fgr(1:nr,1) = r1(1:nr,starti+iso)*r2(1:nr,startj+jso)
                           fga(1:na,1) = a1(1:na,starti+iso,j)*a2(1:na,startj+jso,k)

                           CALL dgemm('N', 'N', na, 1, nr, 1.0_dp, V, na, fgr, nr, 0.0_dp, work, na)
                           CALL dgemm('T', 'N', 1, 1, na, 1.0_dp, work, na, fga, na, 1.0_dp, &
                                      intso(starti+iso, startj+jso, i), 1)

                           ! add r2*a2(j) * V * r1*a1(k)
                           fgr(1:nr,1) = r2(1:nr,starti+iso)*r1(1:nr,startj+jso) 
                           fga(1:na,1) = a2(1:na,starti+iso,j)*a1(1:na,startj+jso,k)

                           CALL dgemm('N', 'N', na, 1, nr, 1.0_dp, V, na, fgr, nr, 0.0_dp, work, na)
                           CALL dgemm('T', 'N', 1, 1, na, 1.0_dp, work, na, fga, na, 1.0_dp, &
                                     intso(starti+iso, startj+jso, i), 1)

                           ! add the last term: r2*a2(j) * V * r2*a2(k)
                           fgr(1:nr,1) = r2(1:nr,starti+iso)*r2(1:nr,startj+jso)
                           fga(1:na,1) = a2(1:na,starti+iso,j)*a2(1:na,startj+jso,k)

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

      DO i = 1,3
         intso(:,:,i) = intso(:,:,i) - TRANSPOSE(intso(:,:,i))
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

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_soc
      TYPE(xas_atom_env_type), POINTER                :: xas_atom_env
      TYPE(qs_environment_type), POINTER              :: qs_env
      
      CHARACTER(len=*), PARAMETER :: routineN = "integrate_soc_atoms", &
                                     routineP = moduleN//":"//routineN

      INTEGER                                         :: ikind, nkind, nset, maxso, na, nr, ir,&
                                                         iset, jset, nsoi, nsoj, sgfi, sgfj, nsgf, &
                                                         starti, startj, iat, jat, blk, i
      TYPE(cp_3d_r_p_type), DIMENSION(:), POINTER     :: int_soc
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE         :: intso
      TYPE(gto_basis_set_type), POINTER               :: basis
      TYPE(qs_kind_type), DIMENSION(:), POINTER       :: qs_kind_set
      INTEGER, DIMENSION(:), POINTER                  :: npgf, lmax, nsgf_set
      INTEGER, DIMENSION(:,:), POINTER                :: first_sgf
      REAL(dp), DIMENSION(:,:), POINTER               :: sphi_so
      REAL(dp), DIMENSION(:,:,:), POINTER             :: intsgf
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(dbcsr_iterator_type)                       :: iter
      TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
      REAL(dp)                                        :: cspeed, zeff
      REAL(dp), DIMENSION(:,:), ALLOCATABLE           :: V
      REAL(dp), DIMENSION(:), ALLOCATABLE             :: Vr
      TYPE(grid_atom_type), POINTER                   :: grid
   
      NULLIFY(int_soc, basis, qs_kind_set, npgf, lmax, nsgf_set, first_sgf, sphi_so, matrix_s)
      NULLIFY(particle_set)

!  Initialization
      CALL get_qs_env(qs_env, nkind=nkind, qs_kind_set=qs_kind_set, matrix_s=matrix_s, &
                      particle_set=particle_set)
      ALLOCATE(int_soc(nkind))
      cspeed = 1.0_dp/a_fine

!  Loop over the kinds to compute the integrals
      DO ikind = 1,nkind
         
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis, basis_type="ORB", zeff=zeff)
         CALL get_gto_basis_set(basis, nset=nset, maxso=maxso)
         ALLOCATE(intso(nset*maxso, nset*maxso, 3))

         ! compute the model potential on the grid
         grid => xas_atom_env%grid_atom_set(ikind)%grid_atom
         nr = grid%nr 
         na = grid%ng_sphere
         ALLOCATE(Vr(nr))
         CALL calculate_model_potential(Vr, grid, zeff)

         !Compute V/(4c^2-2V) and weight it
         ALLOCATE(V(na,nr))
         DO ir = 1,nr
            V(1:na,ir) = Vr(ir)/(4.0_dp*cspeed**2-2.0_dp*Vr(ir)) *grid%weight(1:na,ir)
         END DO
         DEALLOCATE(Vr)

         ! compute the integral <da_dx|...|db_dy> in terms of so
         CALL integrate_so_dxdy_prod(intso, V, ikind, xas_atom_env, qs_env)
         DEALLOCATE(V)

         ! contract in terms of sgf
         CALL get_gto_basis_set(basis, nsgf=nsgf, nsgf_set=nsgf_set, first_sgf=first_sgf, npgf=npgf, &
                            lmax=lmax)
         ALLOCATE(int_soc(ikind)%array(nsgf,nsgf,3)) 
         intsgf => int_soc(ikind)%array
         sphi_so => xas_atom_env%orb_sphi_so(ikind)%array
         intsgf = 0.0_dp

         DO iset = 1,nset
            starti = (iset-1)*maxso+1
            nsoi = npgf(iset)*nsoset(lmax(iset))
            sgfi = first_sgf(1, iset)

            DO jset= 1, nset
               startj = (jset-1)*maxso+1
               nsoj = npgf(jset)*nsoset(lmax(jset))
               sgfj = first_sgf(1, jset)

               DO i = 1,3
                  CALL ab_contract(intsgf(sgfi:sgfi+nsgf_set(iset)-1, sgfj:sgfj+nsgf_set(jset)-1,i), &
                                   intso(starti:starti+nsoi-1, startj:startj+nsoj-1,i), &
                                   sphi_so(:,sgfi:), sphi_so(:,sgfj:), nsoi, nsoj, &
                                   nsgf_set(iset), nsgf_set(jset))
               END DO
            END DO !jset
         END DO !iset

         DEALLOCATE(intso)
      END DO !ikind

!  Build the matrix_soc based on the matrix_s (but anti-symmetric)
      DO i = 1,3
         CALL dbcsr_create(matrix_soc(i)%matrix, name="SOC MATRIX", template=matrix_s(1)%matrix, &
                           matrix_type="A")
      END DO

!  Iterate over its diagonal blocks and fill=it
      CALL dbcsr_iterator_start(iter, matrix_s(1)%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iat, column=jat, blk=blk)
         IF (.NOT. iat == jat) CYCLE
         ikind = particle_set(iat)%atomic_kind%kind_number

         DO i = 1,3
            CALL dbcsr_put_block(matrix_soc(i)%matrix, iat, iat, int_soc(ikind)%array(:,:,i))
         END DO
      
      END DO !iat
      CALL dbcsr_iterator_stop(iter)
      DO i = 1,3
         CALL dbcsr_finalize(matrix_soc(i)%matrix)
      END DO

! Clean-up
      DO ikind = 1,nkind
         DEALLOCATE(int_soc(ikind)%array)
      END DO
      DEALLOCATE(int_soc)

   END SUBROUTINE integrate_soc_atoms

END MODULE xas_tdp_atom
