!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Couple the PW augmentation density to the compensated Gaussian GAPW residual.
!>        The same cached tensor supplies the energy, both density adjoints and pair derivatives.
! **************************************************************************************************
MODULE qs_uspp_gapw
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE gapw_coulomb_integrals,          ONLY: gapw_coulomb_init,&
                                              gapw_periodic_augmentation_integrals
   USE gapw_coulomb_types,              ONLY: gapw_coulomb_cache_type,&
                                              gapw_coulomb_kind_type,&
                                              gapw_coulomb_term_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indso,&
                                              nsoset
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_poisson_types,                ONLY: pw_poisson_analytic,&
                                              pw_poisson_hockney,&
                                              pw_poisson_mt,&
                                              pw_poisson_multipole,&
                                              pw_poisson_periodic,&
                                              pw_poisson_type,&
                                              pw_poisson_wavelet
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_local_rho_types,              ONLY: local_rho_type
   USE qs_oce_methods,                  ONLY: prj_scatter
   USE qs_rho_atom_types,               ONLY: rho_atom_type
   USE qs_uspp_methods,                 ONLY: qs_uspp_atom_grid_apply_ks,&
                                              qs_uspp_atom_grid_prepare,&
                                              qs_uspp_atom_grid_type
   USE qs_uspp_types,                   ONLY: qs_uspp_type
   USE util,                            ONLY: get_limit
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_uspp_gapw'
   PUBLIC :: qs_uspp_gapw_coulomb

CONTAINS

! **************************************************************************************************
!> \brief Add J(Q, rho_h-rho_s+rho_z-rho0_h) to the GAPW Hartree functional.
!> \param qs_env Current geometry, local densities and force accumulators.
!> \param hamiltonian KS matrices, already initialized from the core Hamiltonian.
!> \param density AO density matrices by spin and image.
!> \param energy Existing one-center Hartree energy, incremented by the bilinear correction.
!> \param calculate_forces Add explicit center/strain and projector Pulay derivatives.
! **************************************************************************************************
   SUBROUTINE qs_uspp_gapw_coulomb(qs_env, hamiltonian, density, energy, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: hamiltonian(:, :)
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:, :)
      REAL(dp), INTENT(INOUT)                            :: energy
      LOGICAL, INTENT(IN)                                :: calculate_forces

      INTEGER                                            :: atom, axis, axis2, bo(2), center, &
                                                            handle, ia, ikind, j, k, m, n, nat, &
                                                            ng, spin, t
      INTEGER, POINTER                                   :: atom_list(:)
      LOGICAL                                            :: paw_atom, use_virial
      REAL(dp)                                           :: alpha_core, coefficient, correction, &
                                                            gradient(3), strain(3, 3), zeff
      REAL(dp), ALLOCATABLE                              :: coefficients(:), &
                                                            gaussian_potential(:, :), hard(:, :), &
                                                            occupations(:, :, :), projected(:), &
                                                            soft(:, :)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(local_rho_type), POINTER                      :: local_rho
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(qs_force_type), POINTER                       :: force(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_uspp_atom_grid_type), ALLOCATABLE          :: fields(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(rho_atom_type), POINTER                       :: rho_atom(:)
      TYPE(virial_type), POINTER                         :: virial

      CALL get_qs_env(qs_env, ks_env=ks_env)
      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL timeset('qs_uspp_gapw_coulomb', handle)
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=kinds, &
                      particle_set=particle_set, cell=cell, para_env=para_env, &
                      local_rho_set=local_rho, rho_atom_set=rho_atom, dft_control=dft_control, &
                      force=force, virial=virial, pw_env=pw_env)
      CALL pw_env_get(pw_env, poisson_env=poisson_env)
      SELECT CASE (poisson_env%method)
      CASE (pw_poisson_periodic, pw_poisson_analytic, pw_poisson_mt, pw_poisson_hockney, &
            pw_poisson_multipole, pw_poisson_wavelet)
         ! These solvers use the free-space Coulomb kernel and its physical images.
      CASE DEFAULT
         CPABORT("USPP/PAW GAPW Coulomb coupling requires a vacuum Poisson solver")
      END SELECT
      IF (poisson_env%parameters%has_dielectric) THEN
         CPABORT("USPP/PAW GAPW Coulomb coupling with a dielectric is not implemented")
      END IF
      IF (ANY(poisson_env%parameters%periodic /= cell%perd)) THEN
         CPABORT("USPP/PAW GAPW requires matching cell and Poisson periodicity")
      END IF
      ! Quantum nuclear densities require an additional product map.
      IF (ASSOCIATED(local_rho%rhoz_cneo_set)) THEN
         CPABORT("USPP GAPW Coulomb coupling with CNEO is not implemented")
      END IF
      IF (ALLOCATED(uspp%coulomb)) THEN
         IF (ANY(uspp%coulomb%hmat /= cell%hmat) .OR. ANY(uspp%coulomb%periodicity /= cell%perd) .OR. &
             uspp%coulomb%eta /= local_rho%rho0_mpole%zet0_h) DEALLOCATE (uspp%coulomb)
      END IF
      IF (ALLOCATED(uspp%coulomb)) THEN
         DO atom = 1, SIZE(particle_set)
            IF (ANY(uspp%coulomb%position(:, atom) /= particle_set(atom)%r)) THEN
               DEALLOCATE (uspp%coulomb)
               EXIT
            END IF
         END DO
      END IF
      IF (.NOT. ALLOCATED(uspp%coulomb)) THEN
         ALLOCATE (uspp%coulomb)
         CALL build_cache(uspp%coulomb, uspp, kinds, atomic_kind_set, particle_set, cell, &
                          para_env, local_rho, dft_control)
      END IF

      ALLOCATE (fields(SIZE(density, 1)))
      DO spin = 1, SIZE(fields)
         CALL qs_uspp_atom_grid_prepare(uspp, kinds, density(spin, :), fields(spin), calculate_forces)
      END DO
      ALLOCATE (occupations, SOURCE=fields(1)%density)
      DO spin = 2, SIZE(fields)
         occupations(:, :, :) = occupations + fields(spin)%density
      END DO
      correction = 0.0_dp
      use_virial = calculate_forces .AND. virial%pv_availability .AND. (.NOT. virial%pv_numer)
      DO ikind = 1, SIZE(kinds)
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=nat)
         CALL get_qs_kind(kinds(ikind), paw_atom=paw_atom, alpha_core_charge=alpha_core, zeff=zeff)
         IF (.NOT. paw_atom .AND. dft_control%qs_control%gapw_control%nopaw_as_gpw) CYCLE
         ASSOCIATE (map => uspp%coulomb%kind(ikind))
            ng = map%nprimitive
            ALLOCATE (hard(ng, ng), soft(ng, ng), gaussian_potential(ng, ng), &
                      coefficients(map%ncomponent), projected(map%ncomponent))
            bo = get_limit(nat, para_env%num_pe, para_env%mepos)
            DO ia = bo(1), bo(2)
               atom = atom_list(ia)
               hard = 0.0_dp
               IF (paw_atom) THEN
                  DO spin = 1, SIZE(fields)
                     soft = 0.0_dp
                     CALL prj_scatter(rho_atom(atom)%cpc_h(spin)%r_coef, soft, kinds(ikind))
                     hard(:, :) = hard + soft
                     soft = 0.0_dp
                     CALL prj_scatter(rho_atom(atom)%cpc_s(spin)%r_coef, soft, kinds(ikind))
                     hard(:, :) = hard - soft
                  END DO
               END IF
               coefficients = 0.0_dp
               DO t = 1, SIZE(map%term)
                  ASSOCIATE (term => map%term(t))
                     DO k = 1, SIZE(term%factor)
                        j = term%row(k)
                        n = term%col(k)
                        coefficient = hard(j, n)
                        IF (j /= n) coefficient = coefficient + hard(n, j)
                        m = term%offset + term%m(k)
                        coefficients(m) = coefficients(m) + coefficient*term%factor(k)
                     END DO
                  END ASSOCIATE
               END DO
               ! The last component is the ionic Gaussian minus its compensation monopole.
               coefficients(map%ncomponent) = -zeff*SQRT(fourpi)*(alpha_core/pi)**1.5_dp
               gaussian_potential = 0.0_dp
               DO center = 1, SIZE(uspp%center_atoms)
                  n = uspp%projector_sizes(center)
                  ASSOCIATE (tensor => uspp%coulomb%pair(atom, center)%integral)
                     DO m = 1, map%ncomponent
                        projected(m) = SUM(occupations(:n, :n, center)*tensor(:, :, m, 1))
                        ! Both spin Hamiltonians see the derivative with respect to the total density.
                        DO spin = 1, SIZE(fields)
                           fields(spin)%potential(:n, :n, center) = fields(spin)%potential(:n, :n, center) + &
                                                                    coefficients(m)*tensor(:, :, m, 1)
                        END DO
                     END DO
                     correction = correction + DOT_PRODUCT(coefficients, projected)
                     IF (calculate_forces) THEN
                        gradient = 0.0_dp
                        strain = 0.0_dp
                        DO m = 1, map%ncomponent
                           DO axis = 1, 3
                              gradient(axis) = gradient(axis) + coefficients(m)* &
                                               SUM(occupations(:n, :n, center)*tensor(:, :, m, 1 + axis))
                              IF (use_virial) THEN
                                 DO axis2 = 1, 3
                                    strain(axis, axis2) = strain(axis, axis2) + coefficients(m)* &
                                                          SUM(occupations(:n, :n, center)*tensor(:, :, m, 4 + axis + 3*(axis2 - 1)))
                                 END DO
                              END IF
                           END DO
                        END DO
                        ! Explicit derivatives are total-density terms and must be added only once.
                        fields(1)%center_gradient(:, center) = fields(1)%center_gradient(:, center) + gradient
                        force(ikind)%gth_ppnl(:, ia) = force(ikind)%gth_ppnl(:, ia) - gradient
                        fields(1)%center_strain = fields(1)%center_strain + strain
                     END IF
                  END ASSOCIATE
                  DO t = 1, SIZE(map%term)
                     ASSOCIATE (term => map%term(t))
                        DO k = 1, SIZE(term%factor)
                           j = term%row(k)
                           n = term%col(k)
                           m = term%offset + term%m(k)
                           coefficient = term%factor(k)*projected(m)
                           gaussian_potential(j, n) = gaussian_potential(j, n) + coefficient
                           IF (j /= n) gaussian_potential(n, j) = gaussian_potential(n, j) + coefficient
                        END DO
                     END ASSOCIATE
                  END DO
               END DO
               IF (paw_atom) THEN
                  DO spin = 1, SIZE(fields)
                     ! update_ks_atom applies hard minus soft, including both OCE Pulay terms.
                     rho_atom(atom)%ga_Vlocal_gb_h(spin)%r_coef = &
                        rho_atom(atom)%ga_Vlocal_gb_h(spin)%r_coef + gaussian_potential
                     rho_atom(atom)%ga_Vlocal_gb_s(spin)%r_coef = &
                        rho_atom(atom)%ga_Vlocal_gb_s(spin)%r_coef + gaussian_potential
                  END DO
               END IF
            END DO
            DEALLOCATE (hard, soft, gaussian_potential, coefficients, projected)
         END ASSOCIATE
      END DO
      CALL para_env%sum(correction)
      energy = energy + correction
      CALL qs_uspp_atom_grid_apply_ks(ks_env, fields, hamiltonian, density)
      CALL timestop(handle)
   END SUBROUTINE qs_uspp_gapw_coulomb

! **************************************************************************************************
!> \brief Construct sparse maps from padded GAPW primitive matrices to Gaussian multipoles.
!> \param kind Gaussian kind, with the same angular expansion used by calculate_rho_atom.
!> \param map Output product map; nuclear monopole occupies the final component.
! **************************************************************************************************
   SUBROUTINE build_product_map(kind, map)
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      TYPE(gapw_coulomb_kind_type), INTENT(OUT)          :: map

      INTEGER                                            :: g1, g2, i, ip1, ip2, is1, is2, iso1, &
                                                            iso2, l, m, maxso, nset, power, t, u
      INTEGER, POINTER                                   :: lmax(:), lmin(:), npgf(:)
      LOGICAL                                            :: found, paw_atom
      REAL(dp)                                           :: alpha, cg
      REAL(dp), POINTER                                  :: zet(:, :)
      TYPE(gapw_coulomb_term_type), ALLOCATABLE          :: terms(:)
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      ALLOCATE (map%term(0))
      map%ncomponent = 1
      CALL get_qs_kind(kind, paw_atom=paw_atom, harmonics=harmonics)
      IF (.NOT. paw_atom) RETURN
      CALL get_qs_kind(kind, basis_set=basis, basis_type='GAPW_1C')
      CALL get_gto_basis_set(basis, lmin=lmin, lmax=lmax, npgf=npgf, nset=nset, maxso=maxso, zet=zet)
      map%nprimitive = nset*maxso
      DO is1 = 1, nset
         DO ip1 = 1, npgf(is1)
            DO iso1 = nsoset(lmin(is1) - 1) + 1, nsoset(lmax(is1))
               g1 = (is1 - 1)*maxso + (ip1 - 1)*nsoset(lmax(is1)) + iso1
               DO is2 = 1, nset
                  DO ip2 = 1, npgf(is2)
                     alpha = zet(ip1, is1) + zet(ip2, is2)
                     DO iso2 = nsoset(lmin(is2) - 1) + 1, nsoset(lmax(is2))
                        g2 = (is2 - 1)*maxso + (ip2 - 1)*nsoset(lmax(is2)) + iso2
                        IF (g1 > g2) CYCLE
                        DO i = 1, harmonics%max_iso_not0
                           cg = harmonics%my_CG(iso1, iso2, i)
                           IF (cg == 0.0_dp) CYCLE
                           l = indso(1, i)
                           power = (indso(1, iso1) + indso(1, iso2) - l)/2
                           m = i - l*l
                           found = .FALSE.
                           DO u = 1, SIZE(map%term)
                              ASSOCIATE (term => map%term(u))
                                 IF (term%alpha /= alpha .OR. term%l /= l .OR. term%power /= power) CYCLE
                              END ASSOCIATE
                              t = u
                              found = .TRUE.
                              EXIT
                           END DO
                           IF (.NOT. found) THEN
                              t = SIZE(map%term) + 1
                              ALLOCATE (terms(t))
                              terms(:t - 1) = map%term
                              terms(t)%alpha = alpha
                              terms(t)%l = l
                              terms(t)%power = power
                              terms(t)%offset = map%ncomponent - 1
                              ALLOCATE (terms(t)%row(0), terms(t)%col(0), terms(t)%m(0), terms(t)%factor(0))
                              map%ncomponent = map%ncomponent + 2*l + 1
                              CALL MOVE_ALLOC(terms, map%term)
                           END IF
                           ASSOCIATE (term => map%term(t))
                              CALL append_product(term, g1, g2, m, cg)
                           END ASSOCIATE
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE build_product_map

! **************************************************************************************************
!> \brief Append one nonzero angular product to its sparse multipole map.
!> \param term Gaussian multipole map.
!> \param row First primitive index.
!> \param col Second primitive index.
!> \param m Multipole component.
!> \param factor Product Gaunt coefficient.
! **************************************************************************************************
   SUBROUTINE append_product(term, row, col, m, factor)
      TYPE(gapw_coulomb_term_type), INTENT(INOUT)        :: term
      INTEGER, INTENT(IN)                                :: row, col, m
      REAL(dp), INTENT(IN)                               :: factor

      INTEGER                                            :: n
      INTEGER, ALLOCATABLE                               :: cols(:), components(:), rows(:)
      REAL(dp), ALLOCATABLE                              :: factors(:)

      n = SIZE(term%factor)
      ALLOCATE (rows(n + 1), cols(n + 1), components(n + 1), factors(n + 1))
      rows(:n) = term%row
      cols(:n) = term%col
      components(:n) = term%m
      factors(:n) = term%factor
      rows(n + 1) = row
      cols(n + 1) = col
      components(n + 1) = m
      factors(n + 1) = factor
      CALL MOVE_ALLOC(rows, term%row)
      CALL MOVE_ALLOC(cols, term%col)
      CALL MOVE_ALLOC(components, term%m)
      CALL MOVE_ALLOC(factors, term%factor)
   END SUBROUTINE append_product

! **************************************************************************************************
!> \brief Cache periodic integrals on the same Gaussian-atom owners as the GAPW local matrices.
!> \param cache Output geometry cache.
!> \param uspp Original augmentation radial data and centers.
!> \param kinds Current Gaussian kinds.
!> \param atomic_kind_set Global per-kind atom lists.
!> \param particles Current positions.
!> \param cell Cell vectors and physical periodicity.
!> \param para_env Gaussian atom distribution.
!> \param local_rho Compensation exponents and represented multipoles.
!> \param dft_control GAPW settings.
! **************************************************************************************************
   SUBROUTINE build_cache(cache, uspp, kinds, atomic_kind_set, particles, cell, para_env, local_rho, dft_control)
      TYPE(gapw_coulomb_cache_type), INTENT(OUT)         :: cache
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
      TYPE(atomic_kind_type), INTENT(IN)                 :: atomic_kind_set(:)
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      TYPE(cell_type), INTENT(IN)                        :: cell
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      TYPE(local_rho_type), INTENT(IN)                   :: local_rho
      TYPE(dft_control_type), INTENT(IN)                 :: dft_control

      INTEGER                                            :: atom, bo(2), center, first, ia, ikind, &
                                                            last, maxl, maxpower, n, nat, t
      INTEGER, POINTER                                   :: atom_list(:)
      LOGICAL                                            :: paw_atom
      REAL(dp)                                           :: alpha_core, displacement(3), tolerance
      REAL(dp), ALLOCATABLE                              :: directions(:, :), weights(:)

      cache%eta = local_rho%rho0_mpole%zet0_h
      cache%hmat = cell%hmat
      cache%periodicity = cell%perd
      ALLOCATE (cache%kind(SIZE(kinds)), cache%pair(SIZE(particles), SIZE(uspp%center_atoms)), &
                cache%position(3, SIZE(particles)))
      DO atom = 1, SIZE(particles)
         cache%position(:, atom) = particles(atom)%r
      END DO
      maxl = 0
      maxpower = 0
      DO ikind = 1, SIZE(kinds)
         CALL build_product_map(kinds(ikind), cache%kind(ikind))
         DO t = 1, SIZE(cache%kind(ikind)%term)
            maxl = MAX(maxl, cache%kind(ikind)%term(t)%l)
            maxpower = MAX(maxpower, cache%kind(ikind)%term(t)%power)
         END DO
      END DO
      CALL gapw_coulomb_init(maxl, maxpower)
      tolerance = MIN(1.e-14_dp, dft_control%qs_control%gapw_control%eps_Vrho0)/ &
                  REAL(MAX(1, SIZE(particles)*SIZE(uspp%center_atoms)), dp)
      DO ikind = 1, SIZE(kinds)
         CALL get_qs_kind(kinds(ikind), paw_atom=paw_atom, alpha_core_charge=alpha_core)
         IF (.NOT. paw_atom .AND. dft_control%qs_control%gapw_control%nopaw_as_gpw) CYCLE
         IF (paw_atom) THEN
            IF (local_rho%rho0_mpole%lmax0_kind(ikind) < indso(1, kinds(ikind)%harmonics%max_iso_not0)) THEN
               CPABORT("Incomplete GAPW multipole compensation")
            END IF
         END IF
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=nat)
         bo = get_limit(nat, para_env%num_pe, para_env%mepos)
         DO ia = bo(1), bo(2)
            atom = atom_list(ia)
            DO center = 1, SIZE(uspp%center_atoms)
               n = uspp%projector_sizes(center)
               ALLOCATE (cache%pair(atom, center)%integral(n, n, cache%kind(ikind)%ncomponent, 13))
               displacement = uspp%position(:, center) - particles(atom)%r
               ! Use the augmentation kind's convergable angular rule, not the Gaussian target's rule.
               CALL angular_rule(kinds(uspp%center_kinds(center)), directions, weights)
               DO t = 1, SIZE(cache%kind(ikind)%term)
                  ASSOCIATE (term => cache%kind(ikind)%term(t))
                     first = term%offset + 1
                     last = term%offset + 2*term%l + 1
                     CALL gapw_periodic_augmentation_integrals(uspp%augmentation(uspp%center_kinds(center)), &
                                                   displacement, term%alpha, term%power, term%l, cache%eta, cell%hmat, cell%h_inv, &
                                           directions, weights, tolerance, cache%pair(atom, center)%integral(:, :, first:last, :), &
                                                               periodicity=cache%periodicity)
                  END ASSOCIATE
               END DO
               last = cache%kind(ikind)%ncomponent
               CALL gapw_periodic_augmentation_integrals(uspp%augmentation(uspp%center_kinds(center)), &
                                            displacement, alpha_core, 0, 0, cache%eta, cell%hmat, cell%h_inv, directions, weights, &
                                                         tolerance, cache%pair(atom, center)%integral(:, :, last:last, :), &
                                                         periodicity=cache%periodicity)
               DEALLOCATE (directions, weights)
            END DO
         END DO
      END DO
   END SUBROUTINE build_cache

! **************************************************************************************************
!> \brief Reuse the kind's Lebedev angular rule with weights normalized to 4*pi.
!> \param kind Initialized GAPW kind.
!> \param directions Cartesian unit vectors.
!> \param weights Solid-angle quadrature weights.
! **************************************************************************************************
   SUBROUTINE angular_rule(kind, directions, weights)
      TYPE(qs_kind_type), INTENT(IN)                     :: kind
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: directions(:, :), weights(:)

      INTEGER                                            :: n

      n = kind%grid_atom%ng_sphere
      ALLOCATE (directions(3, n), weights(n))
      directions(1, :) = kind%grid_atom%sin_pol*kind%grid_atom%cos_azi
      directions(2, :) = kind%grid_atom%sin_pol*kind%grid_atom%sin_azi
      directions(3, :) = kind%grid_atom%cos_pol
      weights(:) = kind%grid_atom%wa*(fourpi/SUM(kind%grid_atom%wa))
   END SUBROUTINE angular_rule

END MODULE qs_uspp_gapw
