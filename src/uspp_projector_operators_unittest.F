!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief MPI regression for AO/projector/grid coupling, with a dense Cartesian Gaussian oracle.
! **************************************************************************************************
PROGRAM uspp_projector_operators_unittest
   USE atom_upf,                        ONLY: atom_release_upf,&
                                              atom_upfpot_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_info, dbcsr_get_stored_coordinates, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, dbcsr_put_block, dbcsr_release, &
        dbcsr_set, dbcsr_type
   USE cp_subsys_types,                 ONLY: cp_subsys_release,&
                                              cp_subsys_type
   USE external_potential_types,        ONLY: allocate_potential,&
                                              deallocate_potential,&
                                              init_potential
   USE fft_tools,                       ONLY: finalize_fft,&
                                              init_fft
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: pi
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_comm_world
   USE orbital_pointers,                ONLY: indco,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: FULLSPACE,&
                                              HALFSPACE,&
                                              pw_grid_type
   USE pw_grids,                        ONLY: pw_grid_create,&
                                              pw_grid_release
   USE pw_methods,                      ONLY: pw_copy,&
                                              pw_derive,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              qs_ks_part_release,&
                                              qs_ks_release,&
                                              set_ks_env
   USE qs_subsys_types,                 ONLY: qs_subsys_set,&
                                              qs_subsys_type
   USE qs_uspp_methods,                 ONLY: &
        qs_uspp_atom_grid_apply, qs_uspp_atom_grid_collocate, qs_uspp_atom_grid_integrate, &
        qs_uspp_atom_grid_prepare, qs_uspp_atom_grid_type, qs_uspp_collocate, qs_uspp_core_matrix, &
        qs_uspp_init_ao_matrix, qs_uspp_integrate
   USE qs_uspp_projectors,              ONLY: qs_uspp_projector_gradient,&
                                              qs_uspp_projector_images,&
                                              qs_uspp_projector_matrix,&
                                              qs_uspp_projector_pair
   USE qs_uspp_types,                   ONLY: qs_uspp_init,&
                                              qs_uspp_type,&
                                              qs_uspp_update
   USE sap_kind_types,                  ONLY: release_sap_int,&
                                              sap_int_type,&
                                              sap_sort
   USE upf_local_integrals,             ONLY: upf_local_fourier,&
                                              upf_local_init,&
                                              upf_local_type
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_type
   USE uspp_augmentation,               ONLY: uspp_augmentation_init,&
                                              uspp_augmentation_type
   USE uspp_projector_operators,        ONLY: uspp_projector_add_operator,&
                                              uspp_projector_collocate,&
                                              uspp_projector_density,&
                                              uspp_projector_integrate,&
                                              uspp_projector_metric,&
                                              uspp_projector_pullback

   IMPLICIT NONE

   INTEGER, PARAMETER :: nao = 6, nproj = 6, ng = 7
   INTEGER, TARGET :: ao_sizes(3) = [2, 3, 1], pp_sizes(2) = [4, 2]
   INTEGER, TARGET :: ao_rows(3), ao_cols(3), pp_rows(2), pp_cols(2)
   INTEGER, ALLOCATABLE, TARGET :: grid(:, :)
   INTEGER, POINTER :: pgrid(:, :), rs(:), cs(:), rd(:), cd(:)
   INTEGER :: i, j, k, ip, iq, rank, nranks, nlocal, ig, case_id, nprows, npcols, idir, center
   INTEGER :: kind(2) = [1, 2], first(2) = [1, 5], last(2) = [4, 6]
   INTEGER, ALLOCATABLE :: local_index(:)
   REAL(KIND=dp), PARAMETER :: volume = 137.0_dp, step = 1.e-5_dp
   REAL(KIND=dp) :: b(nao, nproj), p(nao, nao), sg(nao, nao), metric(nao, nao), work(nao, nao)
   REAL(KIND=dp) :: h(nao, nao), q(nproj, nproj), d(nproj, nproj), a(nproj, nproj), projected(nproj, nproj)
   REAL(KIND=dp) :: db(nao, nproj), gradient(nao, nproj), positions(3, 2), gvecs(3, ng), weights(ng)
   REAL(KIND=dp) :: center_gradient(3, 2), displacement
   REAL(KIND=dp) :: charge, energy, expected, plus, minus, max_error, alpha(2) = [1.1_dp, 0.7_dp]
   REAL(KIND=dp), ALLOCATABLE :: local_g(:, :), local_weight(:)
   COMPLEX(KIND=dp) :: vg(ng), rho_ref(ng), qg(nproj, nproj), rho_all(ng)
   COMPLEX(KIND=dp), ALLOCATABLE :: local_rho(:), local_v(:)
   TYPE(atom_upfpot_type) :: pot(2)
   TYPE(upf_projector_type) :: beta(2)
   TYPE(uspp_augmentation_type), ALLOCATABLE :: aug(:)
   TYPE(dbcsr_type), TARGET :: bm
   TYPE(dbcsr_type) :: pm, sm, gm, hm, dm, am, gradm, assembled
   TYPE(dbcsr_distribution_type) :: dist_ao, dist_b, dist_pp
   TYPE(mp_comm_type) :: comm

   CALL cp2k_init()
   CALL comm%set_handle(mp_comm_world%get_handle())
   rank = comm%mepos
   nranks = comm%num_pe
   CALL init_fft("FFTSG", alltoall=.FALSE., pool_limit=4, wisdom_file="", plan_style=1)
   max_error = 0.0_dp
   ALLOCATE (aug(2))
   DO k = 1, 2
      CALL make_augmentation(k)
   END DO
   positions(:, 1) = [0.21_dp, -0.4_dp, 0.7_dp]
   positions(:, 2) = [1.3_dp, 0.9_dp, -0.8_dp]
   gvecs(:, 1) = 0.0_dp
   gvecs(:, 2) = [0.3_dp, -0.7_dp, 1.1_dp]
   gvecs(:, 3) = [2.1_dp, 0.2_dp, -0.4_dp]
   gvecs(:, 4) = [0.0_dp, 0.0_dp, 2.8_dp]
   gvecs(:, 5:7) = -gvecs(:, 2:4)
   vg(1) = CMPLX(0.37_dp, 0.0_dp, dp)
   vg(2:4) = [CMPLX(0.2_dp, -0.3_dp, dp), CMPLX(-0.5_dp, 0.8_dp, dp), CMPLX(1.1_dp, 0.4_dp, dp)]
   vg(5:7) = CONJG(vg(2:4))
   CALL add_ionic_potential()
   DO j = 1, nproj
      DO i = 1, nao
         b(i, j) = 0.13_dp*SIN(REAL(i*j, dp))
         db(i, j) = 0.17_dp*COS(REAL(2*i + j, dp))
      END DO
   END DO
   b(1:2, 5:6) = 0.0_dp
   b(6, 1:4) = 0.0_dp
   db(1:2, 5:6) = 0.0_dp
   db(6, 1:4) = 0.0_dp
   p = MATMUL(b, TRANSPOSE(b))
   sg = 0.0_dp
   DO i = 1, nao
      p(i, i) = p(i, i) + 0.7_dp
      sg(i, i) = 1.0_dp + 0.03_dp*i
   END DO
   q = 0.0_dp
   DO k = 1, 2
      q(first(k):last(k), first(k):last(k)) = aug(k)%q
   END DO

   ! Exercise both tall and wide MPI processor grids, including more ranks than centers.
   DO case_id = 1, 2
      nprows = nranks
      npcols = 1
      IF (case_id == 2) THEN
         nprows = 1
         npcols = nranks
      END IF
      ALLOCATE (grid(0:nprows - 1, 0:npcols - 1))
      grid(:, :) = RESHAPE([(i, i=0, nranks - 1)], [nprows, npcols])
      pgrid => grid
      DO i = 1, 3
         ao_rows(i) = MODULO(i - 1, nprows)
         ao_cols(i) = MODULO(i, npcols)
      END DO
      DO i = 1, 2
         pp_rows(i) = MODULO(i, nprows)
         pp_cols(i) = MODULO(i - 1, npcols)
      END DO
      rd => ao_rows; cd => ao_cols
      CALL dbcsr_distribution_new(dist_ao, group=comm%get_handle(), pgrid=pgrid, row_dist=rd, col_dist=cd)
      rd => ao_rows; cd => pp_cols
      CALL dbcsr_distribution_new(dist_b, group=comm%get_handle(), pgrid=pgrid, row_dist=rd, col_dist=cd)
      rd => pp_rows; cd => pp_cols
      CALL dbcsr_distribution_new(dist_pp, group=comm%get_handle(), pgrid=pgrid, row_dist=rd, col_dist=cd)
      rs => ao_sizes; cs => pp_sizes
      CALL make_matrix(bm, 'USPP B', dist_b, rs, cs, b, 'N')
      rs => ao_sizes; cs => ao_sizes
      CALL make_matrix(pm, 'USPP P', dist_ao, rs, cs, p, 'S')
      CALL make_matrix(sm, 'Gaussian S', dist_ao, rs, cs, sg, 'S')
      CALL check_sap_assembly()
      CALL check_image_convolution()
      CALL uspp_projector_metric(bm, aug, kind, sm, gm)
      CALL gather(gm, metric)
      work = sg + MATMUL(b, MATMUL(q, TRANSPOSE(b)))
      CALL check(MAXVAL(ABS(metric - work)), 'generalized metric including new AO blocks')
      IF (MAXVAL(ABS(metric(1:2, 3:5))) < 1.e-6_dp) ERROR STOP 'Missing projector-mediated AO coupling'
      CALL gather(sm, work)
      CALL check(MAXVAL(ABS(work - sg)), 'Gaussian overlap preserved')
      CALL uspp_projector_density(bm, pm, dm)
      CALL gather(dm, projected)
      d = MATMUL(TRANSPOSE(b), MATMUL(p, b))
      d(1:4, 5:6) = 0.0_dp
      d(5:6, 1:4) = 0.0_dp
      CALL check(MAXVAL(ABS(projected - d)), 'distributed on-site density, both triangles of P')
      CALL dbcsr_release(dm)

      ! Full and conjugate-pair G lists, then G=0 alone (some MPI ranks have no G).
      DO ip = 1, 3
         iq = ng
         weights = 1.0_dp
         IF (ip == 2) THEN
            iq = 4
            weights(2:4) = 2.0_dp
         ELSE IF (ip == 3) THEN
            iq = 1
         END IF
         nlocal = COUNT([(MODULO(ig - 1, nranks) == rank, ig=1, iq)])
         ALLOCATE (local_index(nlocal), local_g(3, nlocal), local_weight(nlocal), local_rho(nlocal), local_v(nlocal))
         j = 0
         DO ig = 1, iq
            IF (MODULO(ig - 1, nranks) /= rank) CYCLE
            j = j + 1
            local_index(j) = ig
            local_g(:, j) = gvecs(:, ig)
            local_v(j) = vg(ig)
            local_weight(j) = weights(ig)
         END DO
         rho_ref = CMPLX(0.0_dp, 0.0_dp, dp)
         a = 0.0_dp
         DO ig = 1, iq
            CALL fourier_reference(gvecs(:, ig), qg)
            rho_ref(ig) = SUM(d*qg)
            a = a + volume*weights(ig)*REAL(CONJG(qg)*vg(ig), dp)
         END DO
         local_rho = CMPLX(0.0_dp, 0.0_dp, dp)
         CALL uspp_projector_collocate(bm, pm, aug, kind, positions, local_g, volume, local_rho, charge)
         rho_all = CMPLX(0.0_dp, 0.0_dp, dp)
         rho_all(local_index) = local_rho
         CALL comm%sum(rho_all)
         CALL check(volume*MAXVAL(ABS(rho_all - rho_ref)), 'distributed density versus Cartesian Gaussian integrals')
         CALL check(ABS(charge - volume*REAL(rho_all(1), dp)), 'augmentation G=0 charge')
         CALL check(ABS(SUM(p*sg) + charge - SUM(p*metric)), 'total charge equals Tr(P S)')
         CALL dbcsr_copy(hm, sm)
         CALL uspp_projector_integrate(bm, aug, kind, positions, local_g, volume, local_v, local_weight, hm, &
                                       pm, gradm, center_gradient)
         CALL gather(gradm, gradient)
         CALL check(MAXVAL(ABS(gradient - 2.0_dp*MATMUL(p, MATMUL(b, a)))), 'screened projector pullback')
         CALL dbcsr_release(gradm)
         CALL comm%sum(center_gradient)
         DO center = 1, 2
            DO idir = 1, 3
               displacement = positions(idir, center)
               positions(idir, center) = displacement + step
               plus = 0.0_dp
               DO ig = 1, iq
                  CALL fourier_reference(gvecs(:, ig), qg)
                  plus = plus + volume*weights(ig)*REAL(SUM(d*CONJG(qg))*vg(ig), dp)
               END DO
               positions(idir, center) = displacement - step
               minus = 0.0_dp
               DO ig = 1, iq
                  CALL fourier_reference(gvecs(:, ig), qg)
                  minus = minus + volume*weights(ig)*REAL(SUM(d*CONJG(qg))*vg(ig), dp)
               END DO
               positions(idir, center) = displacement
               CALL check(ABS(center_gradient(idir, center) - (plus - minus)/(2.0_dp*step)), &
                          'augmentation-center derivative versus independent Gaussian energy', 2.e-9_dp)
            END DO
         END DO
         CALL gather(hm, h)
         CALL check(MAXVAL(ABS(h - sg - MATMUL(b, MATMUL(a, TRANSPOSE(b))))), 'screened Hamiltonian')
         energy = volume*SUM(local_weight*REAL(CONJG(local_v)*local_rho, dp))
         CALL comm%sum(energy)
         CALL check(ABS(energy - SUM(p*(h - sg))), 'AO/grid adjoint identity')
         IF (ip == 3) CALL check(MAXVAL(ABS(h - sg - REAL(vg(1), dp)*(metric - sg))), 'constant-potential metric shift')
         ! Repeat collocation to test accumulation and stable, independent spin use.
         CALL uspp_projector_collocate(bm, pm, aug, kind, positions, local_g, volume, local_rho, charge)
         IF (nlocal > 0) CALL check(volume*MAXVAL(ABS(local_rho - 2.0_dp*rho_ref(local_index))), 'density accumulation')
         CALL dbcsr_release(hm)
         DEALLOCATE (local_index, local_g, local_weight, local_rho, local_v)
      END DO

      ! Variational derivative with respect to P, through the sparse AO operator.
      rs => pp_sizes; cs => pp_sizes
      CALL make_matrix(am, 'USPP onsite', dist_pp, rs, cs, a, 'N')
      CALL dbcsr_copy(hm, sm)
      CALL uspp_projector_add_operator(bm, am, hm)
      CALL gather(hm, h)
      work = 0.5_dp*(db + TRANSPOSE(db))
      plus = SUM(MATMUL(TRANSPOSE(b), MATMUL(p + step*work, b))*a)
      minus = SUM(MATMUL(TRANSPOSE(b), MATMUL(p - step*work, b))*a)
      expected = SUM(work*(h - sg))
      CALL check(ABS((plus - minus)/(2.0_dp*step) - expected), 'density derivative of augmentation energy', 1.e-9_dp)
      ! B derivative and sign of the generalized-metric Pulay contraction.
      DO ip = 1, 2
         expected = REAL(3 - 2*ip, dp)
         CALL uspp_projector_pullback(bm, pm, am, expected, gradm)
         CALL gather(gradm, gradient)
         plus = expected*SUM(MATMUL(TRANSPOSE(b + step*db), MATMUL(p, b + step*db))*a)
         minus = expected*SUM(MATMUL(TRANSPOSE(b - step*db), MATMUL(p, b - step*db))*a)
         CALL check(ABS((plus - minus)/(2.0_dp*step) - SUM(gradient*db)), 'projector variational derivative', 1.e-9_dp)
         CALL dbcsr_release(gradm)
      END DO
      IF (case_id == 1) CALL check_pw_grid()
      CALL dbcsr_release(am)
      CALL dbcsr_release(hm)
      CALL dbcsr_release(bm)
      CALL dbcsr_release(pm)
      CALL dbcsr_release(sm)
      CALL dbcsr_release(gm)
      CALL dbcsr_release(assembled)
      CALL dbcsr_distribution_release(dist_ao)
      CALL dbcsr_distribution_release(dist_b)
      CALL dbcsr_distribution_release(dist_pp)
      DEALLOCATE (grid)
   END DO
   CALL comm%max(max_error)
   IF (rank == 0) WRITE (*, '(A,I0,A,ES23.15)') 'USPP distributed AO/grid checks: ranks=', nranks, ', maximum error=', max_error
   DO k = 1, 2
      CALL atom_release_upf(pot(k))
      BLOCK
         TYPE(upf_projector_type) :: empty
         beta(k) = empty
      END BLOCK
   END DO
   DEALLOCATE (aug)
   CALL finalize_fft(comm, "")
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Include local potentials on both augmentation centers and a third ordinary ion.
! **************************************************************************************************
   SUBROUTINE add_ionic_potential()
      COMPLEX(KIND=dp)                                   :: analytic(ng), values(ng)
      INTEGER                                            :: center, ig, ir
      REAL(KIND=dp)                                      :: amplitude, exponent, position(3), &
                                                            potential, radius
      TYPE(atom_upfpot_type)                             :: ion
      TYPE(upf_local_type)                               :: local

      ALLOCATE (ion%r(4001), ion%rab(4001), ion%vlocal(4001))
      ion%zion = 3.0_dp
      DO center = 1, 3
         amplitude = (-1.0_dp)**center*(0.5_dp + center)
         exponent = 0.35_dp + 0.6_dp*center
         IF (center <= 2) THEN
            position = positions(:, center)
         ELSE
            position = [-0.3_dp, 1.6_dp, 0.2_dp]
         END IF
         DO ir = 1, SIZE(ion%r)
            radius = REAL(ir - 1, dp)*0.003_dp
            ion%r(ir) = radius
            ion%rab(ir) = 0.003_dp
            potential = amplitude*EXP(-exponent*radius**2)
            IF (radius > 0.0_dp) THEN
               ion%vlocal(ir) = potential - ion%zion*ERF(SQRT(0.9_dp)*radius)/radius
            ELSE
               ion%vlocal(ir) = potential - 2.0_dp*ion%zion*SQRT(0.9_dp/pi)
            END IF
         END DO
         CALL upf_local_init(ion, 0.9_dp, 0, local)
         CALL upf_local_fourier(local, gvecs, position, volume, values)
         DO ig = 1, ng
            analytic(ig) = amplitude*(pi/exponent)**1.5_dp/volume* &
                           EXP(CMPLX(-SUM(gvecs(:, ig)**2)/(4.0_dp*exponent), -DOT_PRODUCT(gvecs(:, ig), position), dp))
         END DO
         CALL check(MAXVAL(ABS(values - analytic)), 'all-center local ionic Fourier potential')
         vg = vg + values
      END DO
      CALL atom_release_upf(ion)
   END SUBROUTINE add_ionic_potential

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************

   SUBROUTINE check_pw_grid()
      INTEGER                                            :: atom, idir, ix, iy, iz, jdir, k, nx, &
                                                            parity, span
      LOGICAL                                            :: odd
      REAL(KIND=dp) :: atom_gradient(3, 4), cell(3, 3), cell_gradient(3, 3), egrid, eminus, eplus, &
         expected_charge, hpw_dense(nao, nao), norm_g2, original, total_charge
      TYPE(atomic_kind_type), TARGET                     :: atomic_kinds(4)
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(dbcsr_type)                                   :: hpw
      TYPE(particle_type)                                :: particles(4)
      TYPE(pw_c1d_gs_type)                               :: derivative_g, reference_g, rho_g, v_g
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(pw_r3d_rs_type)                               :: derivative_r, reference_r, rho_r, v_r
      TYPE(qs_kind_type)                                 :: qs_kinds(4)
      TYPE(qs_ks_env_type)                               :: ks
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(qs_uspp_type), POINTER                        :: no_uspp, uspp
      TYPE(sap_int_type)                                 :: empty_sap(0)

      ! Exercise the production adapter, using borrowed B and owned per-kind radial data.
      ALLOCATE (subsys, cp_subsys)
      CALL qs_subsys_set(subsys, cp_subsys=cp_subsys)
      CALL cp_subsys_release(cp_subsys)
      CALL set_ks_env(ks, subsys=subsys)
      ALLOCATE (uspp)
      NULLIFY (no_uspp)
      uspp%overlaps => bm
      ALLOCATE (uspp%center_atoms, SOURCE=[1, 2])
      ALLOCATE (uspp%center_kinds, SOURCE=kind)
      ALLOCATE (uspp%projector_sizes, SOURCE=pp_sizes)
      ALLOCATE (uspp%position, SOURCE=positions)
      ALLOCATE (uspp%augmentation, SOURCE=aug)
      CALL set_ks_env(ks, uspp=uspp)
      DO k = 1, 4
         atomic_kinds(k)%kind_number = k
         particles(k)%atomic_kind => atomic_kinds(k)
         CALL allocate_potential(qs_kinds(k)%sgp_potential)
         atom = MIN(k, 2)
         ALLOCATE (qs_kinds(k)%sgp_potential%upf_data, SOURCE=pot(atom))
         qs_kinds(k)%sgp_potential%upf_projectors = beta(atom)
         qs_kinds(k)%sgp_potential%alpha_core_charge = 0.8_dp
      END DO
      qs_kinds(3)%sgp_potential%upf_data%is_ultrasoft = .FALSE.
      qs_kinds(4)%ghost = .TRUE.
      particles(1)%r = positions(:, 1)
      particles(2)%r = positions(:, 2)
      particles(3)%r = [-0.6_dp, 1.4_dp, 0.3_dp]
      particles(4)%r = [0.9_dp, -0.1_dp, 1.2_dp]
      expected_charge = SUM(MATMUL(TRANSPOSE(b), MATMUL(p, b))*q)
      NULLIFY (pw_grid)
      cell = 0.0_dp
      cell(1, 1) = 5.0_dp; cell(2, 2) = 6.0_dp; cell(3, 3) = 4.5_dp
      cell(1, 2) = 0.7_dp
      DO span = FULLSPACE, HALFSPACE, -1
         DO parity = 1, 2
            ! Quickstep uses HALFSPACE only with odd grids (pw_env_methods).
            IF (span == HALFSPACE .AND. parity == 2) CYCLE
            odd = parity == 1
            nx = 4 + parity
            CALL pw_grid_create(pw_grid, comm, cell, grid_span=span, npts=[nx, nx, nx], &
                                odd=odd, fft_usage=.TRUE., rs_dims=[nranks, 1])
            CALL rho_g%create(pw_grid)
            CALL v_g%create(pw_grid)
            CALL derivative_g%create(pw_grid)
            CALL reference_g%create(pw_grid)
            CALL rho_r%create(pw_grid)
            CALL v_r%create(pw_grid)
            CALL derivative_r%create(pw_grid)
            CALL reference_r%create(pw_grid)
            rho_r%array = 0.19_dp
            CALL pw_transfer(rho_r, rho_g)
            total_charge = -0.19_dp*pw_grid%vol
            CALL qs_uspp_collocate(ks, pm, rho_g, rho_r, total_charge)
            CALL check(ABS(total_charge + 0.19_dp*pw_grid%vol + expected_charge), 'signed Quickstep charge')
            egrid = SUM(rho_r%array)*pw_grid%dvol
            CALL comm%sum(egrid)
            CALL check(ABS(egrid + total_charge), 'CP2K FFT integrated augmentation charge')
            CALL pw_transfer(rho_r, reference_g)
            CALL check(MAXVAL(ABS(reference_g%array - rho_g%array)), 'consistent real and reciprocal density')
            DO idir = 1, 3
               derivative_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
               CALL qs_uspp_collocate(ks, pm, derivative_g, derivative_r, &
                                      derivative=[MERGE(1, 0, idir == 1), MERGE(1, 0, idir == 2), MERGE(1, 0, idir == 3)])
               CALL pw_copy(rho_g, reference_g)
               CALL pw_derive(reference_g, [MERGE(1, 0, idir == 1), MERGE(1, 0, idir == 2), MERGE(1, 0, idir == 3)])
               CALL pw_transfer(reference_g, reference_r)
               CALL check(MAXVAL(ABS(reference_r%array - derivative_r%array)), 'collocated augmentation gradient')
            END DO
            DO iz = LBOUND(v_r%array, 3), UBOUND(v_r%array, 3)
               DO iy = LBOUND(v_r%array, 2), UBOUND(v_r%array, 2)
                  DO ix = LBOUND(v_r%array, 1), UBOUND(v_r%array, 1)
                     v_r%array(ix, iy, iz) = 0.37_dp + SIN(0.17_dp*REAL(ix*iy, dp)) + COS(0.23_dp*REAL(iz + 2*ix, dp))
                  END DO
               END DO
            END DO
            CALL dbcsr_copy(hpw, sm)
            CALL qs_uspp_integrate(ks, v_r, hpw, weighted=.FALSE.)
            CALL gather(hpw, hpw_dense)
            egrid = SUM(v_r%array*(rho_r%array - 0.19_dp))*pw_grid%dvol
            CALL comm%sum(egrid)
            CALL check(ABS(egrid - SUM(p*(hpw_dense - sg))), 'Quickstep real-space FFT/DBCSR adjoint identity')
            CALL dbcsr_copy(hpw, sm)
            CALL pw_scale(v_r, pw_grid%dvol)
            CALL qs_uspp_integrate(ks, v_r, hpw, weighted=.TRUE.)
            CALL gather(hpw, hpw_dense)
            CALL check(ABS(egrid - SUM(p*(hpw_dense - sg))), 'Quickstep weighted potential normalization')
            ! The core contribution must include the third, ordinary ion and exclude the ghost.
            v_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
            DO atom = 1, 3
               k = MIN(atom, 2)
               DO ig = 1, SIZE(v_g%array)
                  norm_g2 = SUM(pw_grid%g(:, ig)**2)
                  v_g%array(ig) = v_g%array(ig) + REAL(k, dp)*(pi/(0.6_dp + 0.2_dp*k))**1.5_dp &
                                  *EXP(-norm_g2/(4.0_dp*(0.6_dp + 0.2_dp*k)))/pw_grid%vol &
                                  *EXP(CMPLX(0.0_dp, -DOT_PRODUCT(pw_grid%g(:, ig), particles(atom)%r), dp))
               END DO
            END DO
            CALL pw_transfer(v_g, v_r)
            CALL dbcsr_copy(hpw, sm)
            CALL qs_uspp_core_matrix(uspp, qs_kinds, particles, pw_grid, hpw)
            CALL gather(hpw, hpw_dense)
            egrid = SUM(v_r%array*(rho_r%array - 0.19_dp))*pw_grid%dvol
            CALL comm%sum(egrid)
            DO k = 1, 2
               egrid = egrid + SUM(MATMUL(TRANSPOSE(b(:, first(k):last(k))), &
                                          MATMUL(p, b(:, first(k):last(k))))*beta(k)%d)
            END DO
            CALL check(ABS(egrid - SUM(p*(hpw_dense - sg))), 'bare D0 plus all-ion augmentation core energy')
            atom_gradient = 0.0_dp
            cell_gradient = 0.0_dp
            CALL dbcsr_copy(hpw, sm)
            CALL qs_uspp_core_matrix(uspp, qs_kinds, particles, pw_grid, hpw, pm, pm, empty_sap, atom_gradient, cell_gradient)
            CALL comm%sum(atom_gradient)
            CALL comm%sum(cell_gradient)
            DO jdir = 1, 3
               DO idir = 1, 3
                  CALL core_strain_energy(uspp, qs_kinds, particles, cell, span, nx, idir, jdir, step, eplus)
                  CALL core_strain_energy(uspp, qs_kinds, particles, cell, span, nx, idir, jdir, -step, eminus)
                  CALL check(ABS(cell_gradient(idir, jdir) - (eplus - eminus)/(2.0_dp*step)), &
                             'FFT augmentation and all-ion strain derivative', 3.e-9_dp)
               END DO
            END DO
            DO atom = 1, 4
               DO idir = 1, 3
                  original = particles(atom)%r(idir)
                  particles(atom)%r(idir) = original + step
                  IF (atom <= 2) uspp%position(idir, atom) = original + step
                  CALL dbcsr_copy(hpw, sm)
                  CALL qs_uspp_core_matrix(uspp, qs_kinds, particles, pw_grid, hpw)
                  CALL gather(hpw, hpw_dense)
                  eplus = SUM(p*hpw_dense)
                  particles(atom)%r(idir) = original - step
                  IF (atom <= 2) uspp%position(idir, atom) = original - step
                  CALL dbcsr_copy(hpw, sm)
                  CALL qs_uspp_core_matrix(uspp, qs_kinds, particles, pw_grid, hpw)
                  CALL gather(hpw, hpw_dense)
                  eminus = SUM(p*hpw_dense)
                  particles(atom)%r(idir) = original
                  IF (atom <= 2) uspp%position(idir, atom) = original
                  CALL check(ABS(atom_gradient(idir, atom) - (eplus - eminus)/(2.0_dp*step)), &
                             'FFT core augmentation center and all-ion derivative', 3.e-9_dp)
               END DO
            END DO
            CALL dbcsr_copy(hpw, sm)
            CALL qs_uspp_core_matrix(uspp, qs_kinds, particles, pw_grid, hpw)
            CALL gather(hpw, hpw_dense)
            ! Ordinary calculations must leave their arrays and operators untouched.
            CALL set_ks_env(ks, uspp=no_uspp)
            CALL pw_copy(rho_g, reference_g)
            CALL qs_uspp_collocate(ks, pm, rho_g, rho_r, total_charge)
            CALL check(MAXVAL(ABS(reference_g%array - rho_g%array)), 'ordinary density unchanged')
            CALL qs_uspp_integrate(ks, v_r, hpw, weighted=.FALSE.)
            CALL gather(hpw, work)
            CALL check(MAXVAL(ABS(work - hpw_dense)), 'ordinary Hamiltonian unchanged')
            CALL set_ks_env(ks, uspp=uspp)
            CALL dbcsr_release(hpw)
            CALL rho_g%release()
            CALL v_g%release()
            CALL derivative_g%release()
            CALL reference_g%release()
            CALL rho_r%release()
            CALL v_r%release()
            CALL derivative_r%release()
            CALL reference_r%release()
            CALL pw_grid_release(pw_grid)
         END DO
      END DO
      CALL set_ks_env(ks, uspp=no_uspp)
      NULLIFY (uspp%overlaps)
      DEALLOCATE (uspp)
      CALL qs_ks_release(ks)
      DO k = 1, 4
         CALL deallocate_potential(qs_kinds(k)%sgp_potential)
      END DO
   END SUBROUTINE check_pw_grid

! **************************************************************************************************
!> \brief Independent core energy on a strained FFT grid with fixed fractional centers and B.
!> \param uspp Unstrained augmentation state; restored on return.
!> \param qs_kinds Original ionic potentials.
!> \param particles Unstrained particle coordinates.
!> \param cell Unstrained cell.
!> \param span Reciprocal representation.
!> \param nx Grid size.
!> \param idir Cartesian displacement direction.
!> \param jdir Cartesian coordinate supplying the displacement.
!> \param amount Strain component.
!> \param energy Strained core energy at fixed AO density and projector overlaps.
! **************************************************************************************************
   SUBROUTINE core_strain_energy(uspp, qs_kinds, particles, cell, span, nx, idir, jdir, amount, energy)
      TYPE(qs_uspp_type), INTENT(INOUT)                  :: uspp
      TYPE(qs_kind_type), INTENT(INOUT)                  :: qs_kinds(:)
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      REAL(KIND=dp), INTENT(IN)                          :: cell(3, 3)
      INTEGER, INTENT(IN)                                :: span, nx, idir, jdir
      REAL(KIND=dp), INTENT(IN)                          :: amount
      REAL(KIND=dp), INTENT(OUT)                         :: energy

      INTEGER                                            :: atom, axis
      REAL(KIND=dp)                                      :: deformation(3, 3), dense(nao, nao), &
                                                            strained_cell(3, 3)
      REAL(KIND=dp), ALLOCATABLE                         :: saved_position(:, :)
      TYPE(dbcsr_type)                                   :: hcore
      TYPE(particle_type)                                :: moved(SIZE(particles))
      TYPE(pw_grid_type), POINTER                        :: strained_grid

      deformation = 0.0_dp
      DO axis = 1, 3
         deformation(axis, axis) = 1.0_dp
      END DO
      deformation(idir, jdir) = deformation(idir, jdir) + amount
      strained_cell = MATMUL(deformation, cell)
      moved = particles
      DO atom = 1, SIZE(particles)
         moved(atom)%r = MATMUL(deformation, particles(atom)%r)
      END DO
      ALLOCATE (saved_position, SOURCE=uspp%position)
      uspp%position(:, :) = MATMUL(deformation, saved_position)
      NULLIFY (strained_grid)
      CALL pw_grid_create(strained_grid, comm, strained_cell, grid_span=span, npts=[nx, nx, nx], &
                          odd=MODULO(nx, 2) == 1, fft_usage=.TRUE., rs_dims=[nranks, 1])
      CALL dbcsr_copy(hcore, sm)
      CALL qs_uspp_core_matrix(uspp, qs_kinds, moved, strained_grid, hcore)
      CALL gather(hcore, dense)
      energy = SUM(p*(dense - sg))
      CALL dbcsr_release(hcore)
      CALL pw_grid_release(strained_grid)
      uspp%position(:, :) = saved_position
   END SUBROUTINE core_strain_energy

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE check_sap_assembly()
      INTEGER                                            :: atom, axis, center, ENTRY, idir, image, &
                                                            jdir, moved, mycol, myrow, r
      REAL(KIND=dp) :: assembled_b(nao, nproj), atom_gradient(3, 3), cell_gradient(3, 3), &
         distance(3), factor, fd, perturb(nao, nproj)
      TYPE(sap_int_type), POINTER                        :: sap(:)

      ALLOCATE (sap(1))
      ALLOCATE (sap(1)%alist(3))
      CALL dbcsr_get_info(sm, my_prow=myrow, my_pcol=mycol)
      r = 1
      DO atom = 1, 3
         ! Emulate ABBA's union of the local AO row and column atom lists.
         IF (ao_rows(atom) /= myrow .AND. ao_cols(atom) /= mycol) THEN
            r = r + ao_sizes(atom)
            CYCLE
         END IF
         sap(1)%alist(atom)%aatom = atom
         ALLOCATE (sap(1)%alist(atom)%clist(6))
         sap(1)%alist(atom)%clist(:)%catom = 0
         DO center = 1, 2
            IF (atom == 1 .AND. center == 2) CYCLE
            IF (atom == 3 .AND. center == 1) CYCLE
            DO image = 1, 2
               ENTRY = 2*(center - 1) + image
               sap(1)%alist(atom)%clist(ENTRY)%catom = center
               sap(1)%alist(atom)%clist(ENTRY)%cell = [image - 1, 0, 0]
               sap(1)%alist(atom)%clist(ENTRY)%rac = &
                  [0.4_dp*(center - atom) + 5.0_dp*(image - 1), 0.7_dp*(center - atom), -0.2_dp*(center - atom)]
               ALLOCATE (sap(1)%alist(atom)%clist(ENTRY)%acint(ao_sizes(atom), pp_sizes(center), 4))
               factor = 0.3_dp
               IF (image == 2) factor = 0.7_dp
               sap(1)%alist(atom)%clist(ENTRY)%acint(:, :, 1) = factor*b(r:r + ao_sizes(atom) - 1, first(center):last(center))
               DO idir = 1, 3
                  sap(1)%alist(atom)%clist(ENTRY)%acint(:, :, idir + 1) = &
                     factor*(0.2_dp*idir + 0.1_dp*image)*db(r:r + ao_sizes(atom) - 1, first(center):last(center))
               END DO
            END DO
         END DO
         ! One non-USPP center and one unused neighbor must not enter B.
         sap(1)%alist(atom)%clist(5)%catom = 3
         sap(1)%alist(atom)%clist(6)%catom = 0
         r = r + ao_sizes(atom)
      END DO
      CALL qs_uspp_projector_matrix(sm, sap, [1, 2], pp_sizes, assembled)
      CALL gather(assembled, assembled_b)
      CALL check(MAXVAL(ABS(assembled_b - b)), 'Quickstep ABBA ownership and periodic images')
      atom_gradient = 0.0_dp
      cell_gradient = 0.0_dp
      CALL qs_uspp_projector_gradient(assembled, sap, [1, 2], atom_gradient, cell_gradient)
      CALL comm%sum(atom_gradient)
      CALL comm%sum(cell_gradient)
      DO moved = 1, 3
         DO axis = 1, 3
            perturb = 0.0_dp
            r = 1
            DO atom = 1, 3
               DO center = 1, 2
                  IF (atom == 1 .AND. center == 2) CYCLE
                  IF (atom == 3 .AND. center == 1) CYCLE
                  ! Sum 0.3 and 0.7 weighted image derivatives independently.
                  factor = (0.2_dp*axis + 0.17_dp)*(MERGE(1, 0, moved == atom) - MERGE(1, 0, moved == center))
                  perturb(r:r + ao_sizes(atom) - 1, first(center):last(center)) = &
                     factor*db(r:r + ao_sizes(atom) - 1, first(center):last(center))
               END DO
               r = r + ao_sizes(atom)
            END DO
            fd = (SUM((b + step*perturb)**2) - SUM((b - step*perturb)**2))/(4.0_dp*step)
            CALL check(ABS(atom_gradient(axis, moved) - fd), 'image-resolved projector atomic derivative', 2.e-10_dp)
         END DO
      END DO
      CALL check(MAXVAL(ABS(SUM(atom_gradient, DIM=2))), 'projector translation sum rule')
      DO jdir = 1, 3
         DO axis = 1, 3
            perturb = 0.0_dp
            r = 1
            DO atom = 1, 3
               DO center = 1, 2
                  IF (atom == 1 .AND. center == 2) CYCLE
                  IF (atom == 3 .AND. center == 1) CYCLE
                  DO image = 1, 2
                     distance = [0.4_dp*(atom - center) - 5.0_dp*(image - 1), &
                                 0.7_dp*(atom - center), -0.2_dp*(atom - center)]
                     factor = MERGE(0.3_dp, 0.7_dp, image == 1)*(0.2_dp*axis + 0.1_dp*image)*distance(jdir)
                     perturb(r:r + ao_sizes(atom) - 1, first(center):last(center)) = &
                        perturb(r:r + ao_sizes(atom) - 1, first(center):last(center)) &
                        + factor*db(r:r + ao_sizes(atom) - 1, first(center):last(center))
                  END DO
               END DO
               r = r + ao_sizes(atom)
            END DO
            fd = (SUM((b + step*perturb)**2) - SUM((b - step*perturb)**2))/(4.0_dp*step)
            CALL check(ABS(cell_gradient(axis, jdir) - fd), 'image-resolved projector strain derivative', 2.e-10_dp)
         END DO
      END DO
      CALL check_state(sap)
      CALL release_sap_int(sap)
      ! Feed the assembled matrix to all subsequent AO/grid and variational tests.
      CALL dbcsr_release(bm)
      CALL dbcsr_copy(bm, assembled)
   END SUBROUTINE check_sap_assembly

! **************************************************************************************************
!> \brief Compare image-pair contractions with independent complex Bloch matrix products.
!>        Three non-collinear images per AO/center, mixed centers, an omitted center, unequal
!>        AO blocks and both natural atom orders exercise the real-space storage conventions.
! **************************************************************************************************
   SUBROUTINE check_image_convolution()
      COMPLEX(KIND=dp) :: bgrad(nao, nproj), bk(nao, nproj), coeff(nao, nao), &
         datom(nao, nproj, 3, 3), dcell(nao, nproj, 3, 3), hk(nao, nao), href(nao, nao), phase, &
         pkmat(nao, nao)
      INTEGER :: atom, atom_a, atom_b, axis, center, ENTRY, ia, ib, image, ix, iy, iz, ja, jaxis, &
         jb, kp, moved, mycol, myrow, nc, owner, ra, rb, rcell(3), translations(3, 3, 2, 3)
      REAL(KIND=dp) :: atom_grad(3, 3), cell_grad(3, 3), distance(3, 3, 2, 3), &
         dlocal(nproj, nproj), dref(nproj, nproj), factor, image_error, &
         image_values(nao, nproj, 3, 4), kpoints(3, 4), kvec(3), lattice(3, 3), &
         onsite(nproj, nproj), pair_energy, previous_error, reference_energy, sites(3, 3), value
      REAL(KIND=dp), ALLOCATABLE                         :: pair_h(:, :), pair_p(:, :), &
                                                            pair_reverse(:, :)
      TYPE(sap_int_type), POINTER                        :: sap(:)

      previous_error = max_error
      max_error = 0.0_dp
      sites(:, 1:2) = positions
      sites(:, 3) = [-0.6_dp, 0.8_dp, 1.7_dp]
      lattice(:, 1) = [4.1_dp, 0.2_dp, 0.3_dp]
      lattice(:, 2) = [0.0_dp, 4.7_dp, -0.4_dp]
      lattice(:, 3) = [0.6_dp, 0.0_dp, 5.3_dp]
      kpoints(:, 1) = 0.0_dp
      kpoints(:, 2) = [0.173_dp, -0.219_dp, 0.137_dp]
      kpoints(:, 3) = -kpoints(:, 2)
      kpoints(:, 4) = [0.5_dp, 0.0_dp, 0.0_dp]
      onsite = 0.0_dp
      DO center = 1, 2
         DO ib = first(center), last(center)
            DO ia = first(center), last(center)
               onsite(ia, ib) = 0.09_dp*SIN(REAL(ia + ib + center, dp))
            END DO
            onsite(ib, ib) = onsite(ib, ib) + 0.3_dp
         END DO
      END DO
      image_values = 0.0_dp
      ALLOCATE (sap(1))
      ALLOCATE (sap(1)%alist(3))
      CALL dbcsr_get_info(sm, my_prow=myrow, my_pcol=mycol)
      ra = 1
      DO atom = 1, 3
         ia = ra + ao_sizes(atom) - 1
         sap(1)%alist(atom)%aatom = atom
         IF (ao_rows(atom) == myrow .OR. ao_cols(atom) == mycol) THEN
            ALLOCATE (sap(1)%alist(atom)%clist(8))
            sap(1)%alist(atom)%clist(:)%catom = 0
            sap(1)%alist(atom)%clist(7)%catom = 3
         END IF
         DO center = 1, 2
            IF (atom == 1 .AND. center == 2) CYCLE
            IF (atom == 3 .AND. center == 1) CYCLE
            ja = first(center)
            jb = last(center)
            nc = pp_sizes(center)
            translations(:, atom, center, 1) = 0
            translations(:, atom, center, 2) = [1, -1, MODULO(atom + center, 2)]
            translations(:, atom, center, 3) = [-1, MODULO(atom, 2), -1]
            DO image = 1, 3
               image_values(ra:ia, ja:jb, image, 1) = &
                  0.2_dp*image*b(ra:ia, ja:jb) + 0.02_dp*(image - 1)*db(ra:ia, ja:jb)
               DO axis = 1, 3
                  image_values(ra:ia, ja:jb, image, axis + 1) = &
                     (0.2_dp*axis + 0.1_dp*image)*db(ra:ia, ja:jb)
               END DO
               distance(:, atom, center, image) = &
                  sites(:, center) - sites(:, atom) + MATMUL(lattice, REAL(translations(:, atom, center, image), dp))
               IF (.NOT. ASSOCIATED(sap(1)%alist(atom)%clist)) CYCLE
               ENTRY = 3*(center - 1) + image
               sap(1)%alist(atom)%clist(ENTRY)%catom = center
               sap(1)%alist(atom)%clist(ENTRY)%cell = translations(:, atom, center, image)
               sap(1)%alist(atom)%clist(ENTRY)%rac = distance(:, atom, center, image)
               ALLOCATE (sap(1)%alist(atom)%clist(ENTRY)%acint(ao_sizes(atom), nc, 4))
               sap(1)%alist(atom)%clist(ENTRY)%acint = image_values(ra:ia, ja:jb, image, :)
            END DO
         END DO
         ra = ia + 1
      END DO
      CALL sap_sort(sap)
      DO kp = 1, SIZE(kpoints, 2)
         kvec = kpoints(:, kp)
         bk = CMPLX(0.0_dp, 0.0_dp, dp)
         datom = CMPLX(0.0_dp, 0.0_dp, dp)
         dcell = CMPLX(0.0_dp, 0.0_dp, dp)
         ra = 1
         DO atom = 1, 3
            ia = ra + ao_sizes(atom) - 1
            DO center = 1, 2
               IF (atom == 1 .AND. center == 2) CYCLE
               IF (atom == 3 .AND. center == 1) CYCLE
               ja = first(center)
               jb = last(center)
               DO image = 1, 3
                  phase = EXP(CMPLX(0.0_dp, 2.0_dp*pi*DOT_PRODUCT(kvec, REAL(translations(:, atom, center, image), dp)), dp))
                  bk(ra:ia, ja:jb) = bk(ra:ia, ja:jb) + phase*image_values(ra:ia, ja:jb, image, 1)
                  DO axis = 1, 3
                     datom(ra:ia, ja:jb, axis, atom) = datom(ra:ia, ja:jb, axis, atom) &
                                                       + phase*image_values(ra:ia, ja:jb, image, axis + 1)
                     datom(ra:ia, ja:jb, axis, center) = datom(ra:ia, ja:jb, axis, center) &
                                                         - phase*image_values(ra:ia, ja:jb, image, axis + 1)
                     DO jaxis = 1, 3
                        dcell(ra:ia, ja:jb, axis, jaxis) = dcell(ra:ia, ja:jb, axis, jaxis) &
                                                           - phase*image_values(ra:ia, ja:jb, image, axis + 1) &
                                                           *distance(jaxis, atom, center, image)
                     END DO
                  END DO
               END DO
            END DO
            ra = ia + 1
         END DO
         DO ib = 1, nao
            DO ia = 1, nao
               coeff(ia, ib) = CMPLX(0.2_dp*SIN(0.3_dp*ia + 0.2_dp*ib), 0.1_dp*COS(0.1_dp*ia - 0.4_dp*ib), dp)
            END DO
         END DO
         IF (kp == 1) coeff = CMPLX(REAL(coeff, dp), 0.0_dp, dp)
         IF (kp == 3) coeff = CONJG(coeff)
         pkmat = MATMUL(coeff, CONJG(TRANSPOSE(coeff)))
         href = MATMUL(bk, MATMUL(CMPLX(onsite, 0.0_dp, dp), CONJG(TRANSPOSE(bk))))
         dref = REAL(MATMUL(CONJG(TRANSPOSE(bk)), MATMUL(pkmat, bk)), dp)
         dref(1:4, 5:6) = 0.0_dp
         dref(5:6, 1:4) = 0.0_dp
         bgrad = 2.0_dp*MATMUL(pkmat, MATMUL(bk, CMPLX(onsite, 0.0_dp, dp)))
         hk = CMPLX(0.0_dp, 0.0_dp, dp)
         dlocal = 0.0_dp
         atom_grad = 0.0_dp
         cell_grad = 0.0_dp
         pair_energy = 0.0_dp
         DO atom_a = 1, 3
            DO atom_b = atom_a, 3
               CALL dbcsr_get_stored_coordinates(sm, atom_a, atom_b, owner)
               IF (owner /= rank) CYCLE
               ia = atom_a
               ib = atom_b
               ! The second processor layout uses the reversed natural atom order.
               IF (case_id == 2) THEN
                  ia = atom_b
                  ib = atom_a
               END IF
               ra = 1 + SUM(ao_sizes(:ia - 1))
               rb = 1 + SUM(ao_sizes(:ib - 1))
               ja = ra + ao_sizes(ia) - 1
               jb = rb + ao_sizes(ib) - 1
               ALLOCATE (pair_h(ao_sizes(ia), ao_sizes(ib)), pair_p(ao_sizes(ia), ao_sizes(ib)), &
                         pair_reverse(ao_sizes(ib), ao_sizes(ia)))
               factor = MERGE(1.0_dp, 2.0_dp, ia == ib)
               DO iz = -2, 2
                  DO iy = -2, 2
                     DO ix = -2, 2
                        rcell = [ix, iy, iz]
                        phase = EXP(CMPLX(0.0_dp, 2.0_dp*pi*DOT_PRODUCT(kvec, REAL(rcell, dp)), dp))
                        pair_p(:, :) = REAL(pkmat(ra:ja, rb:jb)*CONJG(phase), dp)
                        pair_h = 0.0_dp
                        pair_reverse = 0.0_dp
                        DO center = 1, 2
                           CALL qs_uspp_projector_pair(sap(1)%alist(ia), sap(1)%alist(ib), center, rcell, &
                                                       onsite=onsite(first(center):last(center), first(center):last(center)), &
                                                       ao_operator=pair_h, density=pair_p, &
                                                       projected=dlocal(first(center):last(center), first(center):last(center)), &
                                                       weight=factor, atom_gradient=atom_grad, cell_gradient=cell_grad)
                           CALL qs_uspp_projector_pair(sap(1)%alist(ib), sap(1)%alist(ia), center, -rcell, &
                                                       onsite=onsite(first(center):last(center), first(center):last(center)), &
                                                       ao_operator=pair_reverse)
                        END DO
                        CALL check(MAXVAL(ABS(pair_h - TRANSPOSE(pair_reverse))), 'image-pair transpose and reversed cell')
                        hk(ra:ja, rb:jb) = hk(ra:ja, rb:jb) + phase*pair_h
                        IF (ia /= ib) hk(rb:jb, ra:ja) = hk(rb:jb, ra:ja) + CONJG(phase)*TRANSPOSE(pair_h)
                        pair_energy = pair_energy + factor*SUM(pair_p*pair_h)
                     END DO
                  END DO
               END DO
               ! A cell outside every image difference has no operator contribution.
               pair_h = 0.0_dp
               CALL qs_uspp_projector_pair(sap(1)%alist(ia), sap(1)%alist(ib), 1, [3, 0, 0], &
                                           onsite=onsite(1:4, 1:4), ao_operator=pair_h)
               CALL check(MAXVAL(ABS(pair_h)), 'empty projector image pair')
               DEALLOCATE (pair_h, pair_p, pair_reverse)
            END DO
         END DO
         CALL comm%sum(hk)
         CALL comm%sum(dlocal)
         CALL comm%sum(atom_grad)
         CALL comm%sum(cell_grad)
         CALL comm%sum(pair_energy)
         CALL check(MAXVAL(ABS(hk - href)), 'image operator versus complex Bloch product')
         CALL check(MAXVAL(ABS(dlocal - dref)), 'image density versus complex Bloch projection')
         reference_energy = REAL(SUM(CONJG(pkmat)*href), dp)
         CALL check(ABS(pair_energy - reference_energy), 'image energy versus complex Bloch trace')
         CALL check(ABS(SUM(onsite*dlocal) - reference_energy), 'image operator-density adjoint')
         DO moved = 1, 3
            DO axis = 1, 3
               value = REAL(SUM(CONJG(bgrad)*datom(:, :, axis, moved)), dp)
               CALL check(ABS(atom_grad(axis, moved) - value), 'image atomic derivative versus complex Bloch pullback')
               value = REAL(SUM(CONJG(bgrad)*dcell(:, :, axis, moved)), dp)
               CALL check(ABS(cell_grad(axis, moved) - value), 'image strain derivative versus complex Bloch pullback')
            END DO
         END DO
         CALL check(MAXVAL(ABS(SUM(atom_grad, DIM=2))), 'image convolution translation sum rule')
         CALL check_image_driver(sap, onsite, pkmat, kvec, href, dref, atom_grad, cell_grad)
      END DO
      CALL release_sap_int(sap)
      image_error = max_error
      CALL comm%max(image_error)
      IF (rank == 0) WRITE (*, '(A,I0,A,ES23.15)') 'USPP image/Bloch checks: layout=', case_id, ', maximum error=', image_error
      max_error = MAX(previous_error, image_error)
   END SUBROUTINE check_image_convolution

! **************************************************************************************************
!> \brief Check distributed image storage and reductions against the complex Bloch oracle.
!> \param sap Image-resolved projector overlaps on the ABBA distribution.
!> \param onsite Dense on-site operator.
!> \param pkmat Complex Bloch density.
!> \param kvec Reduced Bloch vector.
!> \param href Independent Bloch operator.
!> \param dref Independent on-site projected density.
!> \param atom_ref Verified atomic energy gradients.
!> \param cell_ref Verified strain energy gradients.
! **************************************************************************************************
   SUBROUTINE check_image_driver(sap, onsite, pkmat, kvec, href, dref, atom_ref, cell_ref)
      TYPE(sap_int_type), INTENT(IN)                     :: sap(:)
      REAL(KIND=dp), INTENT(IN)                          :: onsite(nproj, nproj)
      COMPLEX(KIND=dp), INTENT(IN)                       :: pkmat(nao, nao)
      REAL(KIND=dp), INTENT(IN)                          :: kvec(3)
      COMPLEX(KIND=dp), INTENT(IN)                       :: href(nao, nao)
      REAL(KIND=dp), INTENT(IN)                          :: dref(nproj, nproj), atom_ref(3, 3), &
                                                            cell_ref(3, 3)

      REAL(KIND=dp), PARAMETER                           :: scale = -0.7_dp

      COMPLEX(KIND=dp)                                   :: hk(nao, nao), phase, stored_phase
      INTEGER                                            :: cells(3, 125), ia, ib, image, ix, iy, &
                                                            iz, ja, jb, pairs(8, 750), ra, rb, slot
      REAL(KIND=dp)                                      :: atom_grad(3, 3), cell_grad(3, 3), &
                                                            dense(nao, nao), &
                                                            projected(nproj, nproj)
      TYPE(dbcsr_p_type)                                 :: h_images(125), p_images(125)
      TYPE(dbcsr_type)                                   :: density, OPERATOR

      rs => pp_sizes; cs => pp_sizes
      CALL make_matrix(OPERATOR, 'image onsite', dist_pp, rs, cs, onsite, 'N')
      rs => ao_sizes; cs => ao_sizes
      slot = 0
      image = 0
      DO iz = -2, 2
         DO iy = -2, 2
            DO ix = -2, 2
               image = image + 1
               cells(:, image) = [ix, iy, iz]
               phase = EXP(CMPLX(0.0_dp, 2.0_dp*pi*DOT_PRODUCT(kvec, REAL(cells(:, image), dp)), dp))
               dense = REAL(pkmat*CONJG(phase), dp)
               DO ia = 1, 3
                  ra = 1 + SUM(ao_sizes(:ia - 1)); ja = ra + ao_sizes(ia) - 1
                  DO ib = ia, 3
                     rb = 1 + SUM(ao_sizes(:ib - 1)); jb = rb + ao_sizes(ib) - 1
                     slot = slot + 1
                     ! Replicate metadata: only the stored AO block owner may contribute.
                     pairs(:, slot) = [ia, ib, 1, 1, image, cells(:, image)]
                     IF (case_id == 2 .AND. ia /= ib) THEN
                        pairs(1:2, slot) = [ib, ia]
                        dense(ra:ja, rb:jb) = REAL(pkmat(ra:ja, rb:jb)*phase, dp)
                     END IF
                  END DO
               END DO
               ALLOCATE (h_images(image)%matrix, p_images(image)%matrix)
               CALL make_matrix(p_images(image)%matrix, 'image density', dist_ao, rs, cs, dense, 'S')
               dense = 1.0_dp
               CALL make_matrix(h_images(image)%matrix, 'image operator', dist_ao, rs, cs, dense, 'S')
               CALL dbcsr_set(h_images(image)%matrix, 0.0_dp)
            END DO
         END DO
      END DO
      CALL check_image_atom_grid(sap, pairs, cells, h_images, p_images, kvec, pkmat, dref)
      atom_grad = 0.0_dp
      cell_grad = 0.0_dp
      CALL qs_uspp_projector_images(bm, sap, pairs, 1, [1, 2], [1, 1], OPERATOR, h_images, p_images, density, &
                                    atom_grad, cell_grad, alpha=scale)
      CALL gather(density, projected)
      CALL comm%sum(atom_grad)
      CALL comm%sum(cell_grad)
      CALL check(MAXVAL(ABS(projected - dref)), 'distributed image density versus Bloch projection')
      CALL check(MAXVAL(ABS(atom_grad - scale*atom_ref)), 'distributed image atomic gradients')
      CALL check(MAXVAL(ABS(cell_grad - scale*cell_ref)), 'distributed image strain gradients')
      hk = CMPLX(0.0_dp, 0.0_dp, dp)
      DO image = 1, SIZE(h_images)
         CALL gather(h_images(image)%matrix, dense)
         phase = EXP(CMPLX(0.0_dp, 2.0_dp*pi*DOT_PRODUCT(kvec, REAL(cells(:, image), dp)), dp))
         DO ia = 1, 3
            ra = 1 + SUM(ao_sizes(:ia - 1)); ja = ra + ao_sizes(ia) - 1
            DO ib = ia, 3
               rb = 1 + SUM(ao_sizes(:ib - 1)); jb = rb + ao_sizes(ib) - 1
               stored_phase = phase
               IF (case_id == 2 .AND. ia /= ib) stored_phase = CONJG(phase)
               hk(ra:ja, rb:jb) = hk(ra:ja, rb:jb) + stored_phase*dense(ra:ja, rb:jb)
               IF (ia /= ib) hk(rb:jb, ra:ja) = hk(rb:jb, ra:ja) + CONJG(stored_phase)*TRANSPOSE(dense(ra:ja, rb:jb))
            END DO
         END DO
         CALL dbcsr_release(h_images(image)%matrix)
         CALL dbcsr_release(p_images(image)%matrix)
         DEALLOCATE (h_images(image)%matrix, p_images(image)%matrix)
      END DO
      CALL check(MAXVAL(ABS(hk - scale*href)), 'distributed image operator versus Bloch product')
      CALL dbcsr_release(OPERATOR)
      CALL dbcsr_release(density)
   END SUBROUTINE check_image_driver

! **************************************************************************************************
!> \brief Exercise the atom-grid adapter with complex Bloch densities and image Hamiltonians.
!> \param sap Projector image lists, borrowed without changing their ownership.
!> \param pairs Natural-order AO image-pair metadata.
!> \param cells Lattice cells of the stored density and Hamiltonian matrices.
!> \param h_images Zero AO operators, reset to zero on return.
!> \param p_images AO image densities.
!> \param kvec Reduced Bloch vector.
!> \param pkmat Independent complex Bloch density.
!> \param dref Independent on-site density from complex projector multiplication.
! **************************************************************************************************
   SUBROUTINE check_image_atom_grid(sap, pairs, cells, h_images, p_images, kvec, pkmat, dref)
      TYPE(sap_int_type), INTENT(IN)                     :: sap(:)
      INTEGER, INTENT(IN)                                :: pairs(:, :), cells(:, :)
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: h_images(:)
      TYPE(dbcsr_p_type), INTENT(IN)                     :: p_images(:)
      REAL(KIND=dp), INTENT(IN)                          :: kvec(3)
      COMPLEX(KIND=dp), INTENT(IN)                       :: pkmat(nao, nao)
      REAL(KIND=dp), INTENT(IN)                          :: dref(nproj, nproj)

      COMPLEX(KIND=dp)                                   :: hk(nao, nao), phase, stored_phase
      INTEGER                                            :: center, ia, ib, image, info, ja, jb, n, &
                                                            np, ra, rb
      REAL(KIND=dp)                                      :: dense(nao, nao), energy, grad(3, 1), &
                                                            points(3, 1), rho(1), vg(3, 1)
      TYPE(cell_type)                                    :: cell
      TYPE(qs_kind_type)                                 :: qs_kinds(2)
      TYPE(qs_uspp_atom_grid_type)                       :: fields
      TYPE(qs_uspp_type)                                 :: state
      TYPE(sap_int_type), TARGET                         :: image_sap(4)

      state%overlaps => bm
      ALLOCATE (state%center_atoms, SOURCE=[1, 2])
      ALLOCATE (state%center_kinds, SOURCE=[1, 2])
      ALLOCATE (state%projector_sizes, SOURCE=pp_sizes)
      ALLOCATE (state%position, SOURCE=positions)
      ALLOCATE (state%augmentation, SOURCE=aug)
      ALLOCATE (state%image_pairs, SOURCE=pairs)
      ! All AO kinds are one in this synthetic fixture; the two projector kinds differ.
      image_sap(1) = sap(1)
      image_sap(3) = sap(1)
      state%image_integrals => image_sap
      DO center = 1, 2
         CALL allocate_potential(qs_kinds(center)%sgp_potential)
         ALLOCATE (qs_kinds(center)%sgp_potential%upf_data, SOURCE=pot(center))
      END DO
      CALL qs_uspp_atom_grid_prepare(state, qs_kinds, p_images, fields, .FALSE.)
      DO center = 1, 2
         n = pp_sizes(center)
         CALL check(MAXVAL(ABS(fields%density(:n, :n, center) - &
                               dref(first(center):last(center), first(center):last(center)))), &
                    'atom-grid image projection versus complex Bloch density')
      END DO
      cell%hmat(:, 1) = [3.1_dp, 0.0_dp, 0.2_dp]
      cell%hmat(:, 2) = [0.4_dp, 3.7_dp, 0.0_dp]
      cell%hmat(:, 3) = [0.0_dp, 0.3_dp, 4.2_dp]
      cell%h_inv = cell%hmat
      CALL invmat(cell%h_inv, info)
      IF (info /= 0) ERROR STOP 'Invalid image atom-grid cell'
      np = MERGE(1, 0, rank == 0)
      points(:, 1) = [0.3_dp, 0.1_dp, 0.7_dp]
      rho = 0.0_dp
      grad = 0.0_dp
      vg(:, 1) = [-0.3_dp, 0.2_dp, 0.7_dp]
      CALL qs_uspp_atom_grid_collocate(state, fields, cell, points(:, :np), rho(:np), grad(:, :np))
      energy = SUM(0.31_dp*(0.2_dp*rho(:np) + SUM(vg(:, :np)*grad(:, :np), DIM=1)))
      CALL comm%sum(energy)
      CALL qs_uspp_atom_grid_integrate(state, fields, cell, points(:, :np), &
                                       SPREAD(0.31_dp, 1, np), SPREAD(0.2_dp, 1, np), vg(:, :np))
      CALL qs_uspp_atom_grid_apply(state, fields, h_images)
      hk = CMPLX(0.0_dp, 0.0_dp, dp)
      DO image = 1, SIZE(h_images)
         CALL gather(h_images(image)%matrix, dense)
         phase = EXP(CMPLX(0.0_dp, 2.0_dp*pi*DOT_PRODUCT(kvec, REAL(cells(:, image), dp)), dp))
         DO ia = 1, 3
            ra = 1 + SUM(ao_sizes(:ia - 1)); ja = ra + ao_sizes(ia) - 1
            DO ib = ia, 3
               rb = 1 + SUM(ao_sizes(:ib - 1)); jb = rb + ao_sizes(ib) - 1
               stored_phase = phase
               IF (case_id == 2 .AND. ia /= ib) stored_phase = CONJG(phase)
               hk(ra:ja, rb:jb) = hk(ra:ja, rb:jb) + stored_phase*dense(ra:ja, rb:jb)
               IF (ia /= ib) hk(rb:jb, ra:ja) = hk(rb:jb, ra:ja) + CONJG(stored_phase)*TRANSPOSE(dense(ra:ja, rb:jb))
            END DO
         END DO
         CALL dbcsr_set(h_images(image)%matrix, 0.0_dp)
      END DO
      CALL check(ABS(energy - REAL(SUM(CONJG(pkmat)*hk), dp)), 'atom-grid complex Bloch energy pairing')
      DO center = 1, 2
         CALL deallocate_potential(qs_kinds(center)%sgp_potential)
      END DO
   END SUBROUTINE check_image_atom_grid

! **************************************************************************************************
!> \brief Check all-center atom quadrature against a Cartesian Gaussian oracle and energy differences.
!> \param uspp Geometry state with two different projector kinds.
!> \param qs_kinds Original UPF kinds.
!> \param sap Synthetic image-resolved projector values and derivatives.
! **************************************************************************************************
   SUBROUTINE check_atom_grid(uspp, qs_kinds, sap)
      TYPE(qs_uspp_type), INTENT(INOUT)                  :: uspp
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kinds(:)
      TYPE(sap_int_type), INTENT(IN)                     :: sap(:)

      INTEGER                                            :: atom, axis, center, image, info, jdir, &
                                                            n, r, sample, sign
      REAL(KIND=dp) :: atom_grad(3, 3), cell_grad(3, 3), dense(nao, nao), distance(3), &
         energies(2), factor, grad(3, 3), grad_ref(3, 3), local_energy, old_positions(3, 2), &
         perturb(nao, nproj), point_change(3, 3), point_grad(3, 3), points(3, 3), &
         projected(nproj, nproj), rho(3), rho_ref(3), shifted(3, 3), vgrad(3, 3), vrho(3), w(3)
      TYPE(cell_type)                                    :: cell, varied_cell
      TYPE(dbcsr_p_type)                                 :: dmats(1), hmats(1)
      TYPE(qs_uspp_atom_grid_type), ALLOCATABLE          :: fields, trial

      ALLOCATE (fields, trial)
      old_positions = uspp%position
      cell%hmat(:, 1) = [3.1_dp, 0.0_dp, 0.2_dp]
      cell%hmat(:, 2) = [0.4_dp, 3.7_dp, 0.0_dp]
      cell%hmat(:, 3) = [0.0_dp, 0.3_dp, 4.2_dp]
      cell%h_inv = cell%hmat
      CALL invmat(cell%h_inv, info)
      IF (info /= 0) ERROR STOP 'Invalid atom-grid test cell'
      points(:, 1) = positions(:, 1)
      points(:, 2) = [0.3_dp, 1.1_dp, -0.7_dp]
      points(:, 3) = [1.7_dp, 0.0_dp, 0.4_dp]
      w = [0.3_dp, -0.5_dp, 0.8_dp]
      DO sample = 1, 3
         IF (MODULO(sample - 1, nranks) /= rank) w(sample) = 0.0_dp
      END DO
      ALLOCATE (dmats(1)%matrix, hmats(1)%matrix)
      CALL dbcsr_copy(dmats(1)%matrix, pm)
      CALL dbcsr_copy(hmats(1)%matrix, sm)
      CALL dbcsr_set(hmats(1)%matrix, 0.0_dp)
      CALL qs_uspp_atom_grid_prepare(uspp, qs_kinds, dmats, fields, .TRUE.)
      projected = MATMUL(TRANSPOSE(b), MATMUL(p, b))
      DO center = 1, 2
         n = pp_sizes(center)
         CALL check(MAXVAL(ABS(fields%density(:n, :n, center) - &
                               projected(first(center):last(center), first(center):last(center)))), &
                    'atom-grid replicated projected density')
      END DO
      rho = 0.4_dp
      grad = 0.07_dp
      CALL qs_uspp_atom_grid_collocate(uspp, fields, cell, points, rho, grad)
      CALL atom_grid_gaussian_reference(cell, points, projected, rho_ref, grad_ref)
      CALL check(MAXVAL(ABS(rho - rho_ref - 0.4_dp)), 'periodic atom-grid Gaussian density oracle', 2.e-10_dp)
      CALL check(MAXVAL(ABS(grad - grad_ref - 0.07_dp)), 'periodic atom-grid Gaussian gradient oracle', 2.e-8_dp)
      shifted = points + SPREAD(MATMUL(cell%hmat, [2.0_dp, -1.0_dp, 1.0_dp]), 2, 3)
      rho_ref = 0.4_dp
      grad_ref = 0.07_dp
      CALL qs_uspp_atom_grid_collocate(uspp, fields, cell, shifted, rho_ref, grad_ref)
      CALL check(MAXVAL(ABS(rho - rho_ref)) + MAXVAL(ABS(grad - grad_ref)), 'atom-grid cell-representative invariance')
      ! Signed quadrature weights and a nonlinear GGA-like energy exercise cross-center terms.
      vrho = rho**2
      vgrad = grad
      point_grad = 0.0_dp
      ! Rank-local point batches are unequal and empty when rank >= three.
      DO sample = rank + 1, 3, nranks
         CALL qs_uspp_atom_grid_integrate(uspp, fields, cell, points(:, sample:sample), w(sample:sample), &
                                          vrho(sample:sample), vgrad(:, sample:sample), point_grad(:, sample:sample))
      END DO
      CALL qs_uspp_atom_grid_collocate(uspp, fields, cell, points(:, :0), rho_ref(:0), grad_ref(:, :0))
      CALL qs_uspp_atom_grid_integrate(uspp, fields, cell, points(:, :0), w(:0), vrho(:0), vgrad(:, :0), &
                                       point_grad(:, :0))
      atom_grad = 0.0_dp
      cell_grad = 0.0_dp
      CALL qs_uspp_atom_grid_apply(uspp, fields, hmats, dmats, sap, atom_grad, cell_grad)
      CALL comm%sum(atom_grad)
      CALL comm%sum(cell_grad)
      CALL comm%sum(point_grad)
      CALL gather(hmats(1)%matrix, dense)
      local_energy = SUM(w*(vrho*(rho - 0.4_dp) + SUM(vgrad*(grad - 0.07_dp), DIM=1)))
      CALL comm%sum(local_energy)
      CALL check(ABS(SUM(p*dense) - local_energy), 'all-center atom-grid AO/field adjoint')
      CALL check(MAXVAL(ABS(SUM(atom_grad, DIM=2) + SUM(point_grad, DIM=2))), 'moving-grid translation sum rule')
      ! A full symmetric AO density direction, including cross-atom blocks.
      DO sign = 1, 2
         CALL dbcsr_release(dmats(1)%matrix)
         rs => ao_sizes; cs => ao_sizes
         CALL make_matrix(dmats(1)%matrix, 'atom-grid trial density', dist_ao, rs, cs, &
                          p + REAL(3 - 2*sign, dp)*step*sg, 'S')
         CALL qs_uspp_atom_grid_prepare(uspp, qs_kinds, dmats, trial, .FALSE.)
         energies(sign) = atom_grid_energy(uspp, trial, cell, points, w)
      END DO
      CALL check(ABS((energies(1) - energies(2))/(2.0_dp*step) - SUM(dense*sg)), &
                 'atom-grid nonlinear AO density derivative', 2.e-9_dp)
      DO atom = 1, 3
         DO axis = 1, 3
            perturb = 0.0_dp
            r = 1
            DO center = 1, 2
               r = 1
               DO sample = 1, 3
                  factor = (0.2_dp*axis + 0.17_dp)*(MERGE(1, 0, atom == sample) - MERGE(1, 0, atom == center))
                  perturb(r:r + ao_sizes(sample) - 1, first(center):last(center)) = &
                     factor*db(r:r + ao_sizes(sample) - 1, first(center):last(center))
                  r = r + ao_sizes(sample)
               END DO
            END DO
            DO sign = 1, 2
               factor = REAL(3 - 2*sign, dp)*step
               trial%density(:, :, :) = fields%density
               projected = MATMUL(TRANSPOSE(b + factor*perturb), MATMUL(p, b + factor*perturb))
               DO center = 1, 2
                  n = pp_sizes(center)
                  trial%density(:n, :n, center) = projected(first(center):last(center), first(center):last(center))
               END DO
               IF (atom <= 2) uspp%position(axis, atom) = old_positions(axis, atom) + factor
               energies(sign) = atom_grid_energy(uspp, trial, cell, points, w)
            END DO
            uspp%position(:, :) = old_positions
            CALL check(ABS((energies(1) - energies(2))/(2.0_dp*step) - atom_grad(axis, atom)), &
                       'atom-grid center plus projector energy derivative', 2.e-9_dp)
            point_change = 0.0_dp
            point_change(axis, atom) = step
            energies(1) = atom_grid_energy(uspp, fields, cell, points + point_change, w)
            energies(2) = atom_grid_energy(uspp, fields, cell, points - point_change, w)
            CALL check(ABS((energies(1) - energies(2))/(2.0_dp*step) - point_grad(axis, atom)), &
                       'atom-grid moving-point energy derivative', 2.e-9_dp)
         END DO
      END DO
      DO jdir = 1, 3
         DO axis = 1, 3
            perturb = 0.0_dp
            r = 1
            DO atom = 1, 3
               DO center = 1, 2
                  DO image = 1, 2
                     distance = [0.4_dp*(atom - center) - 5.0_dp*(image - 1), &
                                 0.7_dp*(atom - center), -0.2_dp*(atom - center)]
                     factor = MERGE(0.3_dp, 0.7_dp, image == 1)*(0.2_dp*axis + 0.1_dp*image)*distance(jdir)
                     perturb(r:r + ao_sizes(atom) - 1, first(center):last(center)) = &
                        perturb(r:r + ao_sizes(atom) - 1, first(center):last(center)) &
                        + factor*db(r:r + ao_sizes(atom) - 1, first(center):last(center))
                  END DO
               END DO
               r = r + ao_sizes(atom)
            END DO
            DO sign = 1, 2
               factor = REAL(3 - 2*sign, dp)*step
               trial%density(:, :, :) = fields%density
               projected = MATMUL(TRANSPOSE(b + factor*perturb), MATMUL(p, b + factor*perturb))
               DO center = 1, 2
                  n = pp_sizes(center)
                  trial%density(:n, :n, center) = projected(first(center):last(center), first(center):last(center))
               END DO
               uspp%position(axis, :) = old_positions(axis, :) + factor*old_positions(jdir, :)
               varied_cell = cell
               varied_cell%hmat(axis, :) = cell%hmat(axis, :) + factor*cell%hmat(jdir, :)
               varied_cell%h_inv = varied_cell%hmat
               CALL invmat(varied_cell%h_inv, info)
               IF (info /= 0) ERROR STOP 'Invalid strained test cell'
               energies(sign) = atom_grid_energy(uspp, trial, varied_cell, points, w)
            END DO
            uspp%position(:, :) = old_positions
            CALL check(ABS((energies(1) - energies(2))/(2.0_dp*step) - cell_grad(axis, jdir)), &
                       'atom-grid source-center and projector strain derivative', 2.e-9_dp)
         END DO
      END DO
      CALL dbcsr_release(dmats(1)%matrix)
      CALL dbcsr_release(hmats(1)%matrix)
      DEALLOCATE (dmats(1)%matrix, hmats(1)%matrix)
   END SUBROUTINE check_atom_grid

! **************************************************************************************************
!> \brief Independent energy used for finite differences of the complete atom-grid map.
!> \param uspp Geometry state.
!> \param fields Projected density.
!> \param cell Auxiliary cell.
!> \param points Fixed Cartesian quadrature points.
!> \param weights Rank-local signed weights.
!> \return Globally reduced nonlinear energy with a fixed density/gradient environment.
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION atom_grid_energy(uspp, fields, cell, points, weights) RESULT(energy)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(qs_uspp_atom_grid_type), INTENT(IN)           :: fields
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(KIND=dp), INTENT(IN)                          :: points(:, :), weights(:)

      REAL(KIND=dp)                                      :: rho(SIZE(weights)), &
                                                            grad(3, SIZE(weights))

      rho = 0.4_dp
      grad = 0.07_dp
      CALL qs_uspp_atom_grid_collocate(uspp, fields, cell, points, rho, grad)
      energy = SUM(weights*(rho**3/3.0_dp + SUM(grad**2, DIM=1)/2.0_dp))
      CALL comm%sum(energy)
   END FUNCTION atom_grid_energy

! **************************************************************************************************
!> \brief Cartesian Gaussian oracle with a fixed, independent image enumeration.
!> \param cell Triclinic test cell.
!> \param points Test points in the central cell.
!> \param density Dense on-site density.
!> \param rho Augmentation density.
!> \param gradient Cartesian density gradient.
! **************************************************************************************************
   SUBROUTINE atom_grid_gaussian_reference(cell, points, density, rho, gradient)
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(KIND=dp), INTENT(IN)                          :: points(:, :), density(nproj, nproj)
      REAL(KIND=dp), INTENT(OUT)                         :: rho(:), gradient(:, :)

      INTEGER                                            :: axis, c, center, i, ix, iy, iz, j, l, m, &
                                                            p, sample, u, v
      REAL(KIND=dp)                                      :: angular(4), angular_d(3, 4), coeff, &
                                                            radial, x(3)

      rho = 0.0_dp
      gradient = 0.0_dp
      DO center = 1, 2
         DO iz = -5, 5
            DO iy = -5, 5
               DO ix = -5, 5
                  DO sample = 1, SIZE(rho)
                     x = points(:, sample) - positions(:, center) - MATMUL(cell%hmat, REAL([ix, iy, iz], dp))
                     radial = EXP(-alpha(center)*SUM(x*x))
                     angular = 0.0_dp
                     angular_d = 0.0_dp
                     DO p = 1, pp_sizes(center)
                        i = beta(center)%radial_index(p)
                        l = pot(center)%lbeta(i)
                        m = beta(center)%angular_index(p) - l*l
                        DO c = 1, nco(l)
                           IF (l == 0) THEN
                              angular(p) = orbtramat(l)%slm(m, c)
                           ELSE
                              DO axis = 1, 3
                                 IF (indco(axis, ncoset(l - 1) + c) == 1) THEN
                                    angular(p) = angular(p) + orbtramat(l)%slm(m, c)*x(axis)
                                    angular_d(axis, p) = orbtramat(l)%slm(m, c)
                                 END IF
                              END DO
                           END IF
                        END DO
                     END DO
                     DO v = 1, pp_sizes(center)
                        j = beta(center)%radial_index(v)
                        DO u = 1, pp_sizes(center)
                           i = beta(center)%radial_index(u)
                           coeff = coefficient(i, j, center)*radial*density(first(center) + u - 1, first(center) + v - 1)
                           rho(sample) = rho(sample) + coeff*angular(u)*angular(v)
                           gradient(:, sample) = gradient(:, sample) + coeff* &
                                                 (angular_d(:, u)*angular(v) + angular(u)*angular_d(:, v) &
                                                  - 2.0_dp*alpha(center)*x*angular(u)*angular(v))
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE atom_grid_gaussian_reference

! **************************************************************************************************
!> \brief Exercise mixed kinds and KS state ownership.
!> \param sap Image-resolved projector integrals.
! **************************************************************************************************
   SUBROUTINE check_state(sap)
      TYPE(sap_int_type), INTENT(INOUT)                  :: sap(:)

      INTEGER                                            :: atom, ENTRY, iteration, k
      REAL(KIND=dp)                                      :: dense(nao, nao), expected(nao, nao)
      TYPE(atomic_kind_type), TARGET                     :: atomic_kinds(4)
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(dbcsr_p_type), POINTER                        :: gaussian(:), images(:, :), metric(:)
      TYPE(dbcsr_type)                                   :: density_template
      TYPE(particle_type)                                :: particles(3)
      TYPE(qs_kind_type)                                 :: qs_kinds(4)
      TYPE(qs_ks_env_type)                               :: ks
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(qs_uspp_type), POINTER                        :: uspp

      NULLIFY (uspp)
      DO k = 1, 4
         atomic_kinds(k)%kind_number = k
      END DO
      ! Kind numbers differ from projector-center numbers; kind 2 is ordinary, 4 a ghost.
      CALL allocate_potential(qs_kinds(3)%sgp_potential)
      ALLOCATE (qs_kinds(3)%sgp_potential%upf_data, SOURCE=pot(1))
      CALL allocate_potential(qs_kinds(1)%sgp_potential)
      ALLOCATE (qs_kinds(1)%sgp_potential%upf_data, SOURCE=pot(2))
      CALL allocate_potential(qs_kinds(4)%sgp_potential)
      ALLOCATE (qs_kinds(4)%sgp_potential%upf_data, SOURCE=pot(1))
      qs_kinds(4)%ghost = .TRUE.
      particles(1)%atomic_kind => atomic_kinds(3)
      particles(2)%atomic_kind => atomic_kinds(1)
      particles(3)%atomic_kind => atomic_kinds(4)
      particles(1)%r = positions(:, 1)
      particles(2)%r = positions(:, 2)
      CALL qs_uspp_init(uspp, qs_kinds, particles(3:))
      IF (ASSOCIATED(uspp)) ERROR STOP 'Ghost atom acquired USPP augmentation'
      particles(3)%atomic_kind => atomic_kinds(2)
      CALL qs_uspp_init(uspp, qs_kinds, particles(3:))
      IF (ASSOCIATED(uspp)) ERROR STOP 'NC atom acquired USPP augmentation'
      CALL init_potential(qs_kinds(3)%sgp_potential)
      IF (qs_kinds(3)%sgp_potential%n_nonlocal /= 0) ERROR STOP 'Unexpected Gaussian nonlocal fit'
      IF (qs_kinds(3)%sgp_potential%nppnl /= beta(1)%nproj) ERROR STOP 'Native projector count missing'
      IF (qs_kinds(3)%sgp_potential%lmax /= MAXVAL(pot(1)%lbeta)) ERROR STOP 'Native projector lmax missing'
      IF (.NOT. ALL(qs_kinds(3)%sgp_potential%is_nonlocal(pot(1)%lbeta))) ERROR STOP 'Native neighbor-list metadata missing'
      CALL qs_uspp_init(uspp, qs_kinds, particles)
      IF (.NOT. ASSOCIATED(uspp)) ERROR STOP 'Missing mixed-kind USPP state'
      IF (ANY(uspp%center_atoms /= [1, 2])) ERROR STOP 'Incorrect augmentation-center mapping'
      IF (ANY(uspp%center_kinds /= [3, 1])) ERROR STOP 'Incorrect augmentation-kind mapping'
      IF (ANY(uspp%projector_sizes /= pp_sizes)) ERROR STOP 'Incorrect magnetic-projector count'
      CALL qs_uspp_update(uspp, sm, sap, particles)
      CALL check_atom_grid(uspp, qs_kinds, sap)

      ALLOCATE (subsys, cp_subsys)
      CALL qs_subsys_set(subsys, cp_subsys=cp_subsys)
      CALL cp_subsys_release(cp_subsys)
      ALLOCATE (gaussian(1))
      ALLOCATE (gaussian(1)%matrix)
      CALL dbcsr_copy(gaussian(1)%matrix, sm)
      CALL set_ks_env(ks, subsys=subsys, matrix_s=gaussian)
      CALL get_ks_env(ks, matrix_metric=metric, matrix_metric_kp=images)
      IF (.NOT. ASSOCIATED(metric, gaussian)) ERROR STOP 'NC metric must alias Gaussian overlap'
      IF (.NOT. ASSOCIATED(images(1, 1)%matrix, gaussian(1)%matrix)) ERROR STOP 'NC image metric mismatch'
      CALL set_ks_env(ks, uspp=uspp)
      CALL get_ks_env(ks, matrix_metric=metric, matrix_metric_kp=images)
      IF (ASSOCIATED(metric(1)%matrix, gaussian(1)%matrix)) ERROR STOP 'USPP metric aliases Gaussian overlap'
      IF (.NOT. ASSOCIATED(metric(1)%matrix, images(1, 1)%matrix)) ERROR STOP 'USPP image metric mismatch'
      CALL gather(metric(1)%matrix, dense)
      expected = sg + MATMUL(b, MATMUL(q, TRANSPOSE(b)))
      CALL check(MAXVAL(ABS(dense - expected)), 'KS generalized metric access')

      ! Higher augmentation multipoles can couple these blocks even when every q vanishes.
      uspp%augmentation(1)%q = 0.0_dp
      uspp%augmentation(3)%q = 0.0_dp
      CALL qs_uspp_update(uspp, sm, sap, particles)
      CALL get_ks_env(ks, matrix_metric=metric)
      CALL gather(metric(1)%matrix, dense)
      CALL check(MAXVAL(ABS(dense - sg)), 'zero integrated augmentation metric')
      CALL gather(uspp%ao_pattern, dense)
      expected = 1.0_dp
      expected(1:2, 6) = 0.0_dp
      expected(6, 1:2) = 0.0_dp
      CALL check(MAXVAL(ABS(dense - expected)), 'projector connectivity independent of integrated q')
      CALL dbcsr_copy(density_template, sm)
      CALL qs_uspp_init_ao_matrix(ks, density_template)
      CALL gather(density_template, dense)
      CALL check(MAXVAL(ABS(dense)), 'expanded density initialized to zero')
      CALL dbcsr_set(density_template, 1.0_dp)
      CALL gather(density_template, dense)
      CALL check(MAXVAL(ABS(dense - expected)), 'density retains all projector-mediated blocks')
      CALL dbcsr_release(density_template)
      uspp%augmentation(1)%q(:, :) = aug(2)%q
      uspp%augmentation(3)%q(:, :) = aug(1)%q

      ! A changed geometry changes B and positions. Repeating that update must not accumulate q.
      DO atom = 1, SIZE(sap(1)%alist)
         IF (.NOT. ASSOCIATED(sap(1)%alist(atom)%clist)) CYCLE
         DO ENTRY = 1, 4
            IF (.NOT. ASSOCIATED(sap(1)%alist(atom)%clist(ENTRY)%acint)) CYCLE
            sap(1)%alist(atom)%clist(ENTRY)%acint = 0.6_dp*sap(1)%alist(atom)%clist(ENTRY)%acint
         END DO
      END DO
      particles(1)%r = particles(1)%r + [0.1_dp, -0.3_dp, 0.2_dp]
      expected = sg + 0.36_dp*MATMUL(b, MATMUL(q, TRANSPOSE(b)))
      DO iteration = 1, 2
         CALL qs_uspp_update(uspp, sm, sap, particles)
         CALL get_ks_env(ks, matrix_metric=metric)
         CALL gather(metric(1)%matrix, dense)
         CALL check(MAXVAL(ABS(dense - expected)), 'rebuilt metric without stale or doubled augmentation')
         CALL check(MAXVAL(ABS(uspp%position(:, 1) - particles(1)%r)), 'updated augmentation-center position')
      END DO
      CALL gather(gaussian(1)%matrix, dense)
      CALL check(MAXVAL(ABS(dense - sg)), 'KS Gaussian overlap preserved across USPP updates')
      CALL qs_ks_part_release(ks)
      CALL get_ks_env(ks, uspp=uspp, matrix_metric=metric)
      IF (ASSOCIATED(uspp)) ERROR STOP 'Partial KS release retained stale USPP state'
      IF (.NOT. ASSOCIATED(metric, gaussian)) ERROR STOP 'Metric fallback after partial release'
      CALL qs_uspp_init(uspp, qs_kinds, particles)
      CALL qs_uspp_update(uspp, sm, sap, particles)
      CALL set_ks_env(ks, uspp=uspp)
      CALL qs_ks_release(ks)
      NULLIFY (uspp, gaussian, metric, images, subsys)
      DO k = 1, 4
         IF (ASSOCIATED(qs_kinds(k)%sgp_potential)) CALL deallocate_potential(qs_kinds(k)%sgp_potential)
      END DO
   END SUBROUTINE check_state

! **************************************************************************************************
!> \brief ...
!> \param value ...
!> \param label ...
!> \param tolerance ...
! **************************************************************************************************
   SUBROUTINE check(value, label, tolerance)
      REAL(KIND=dp), INTENT(IN)                          :: value
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: tolerance

      REAL(KIND=dp)                                      :: tol

      tol = 2.e-12_dp
      IF (PRESENT(tolerance)) tol = tolerance
      max_error = MAX(max_error, value)
      IF (value > tol) THEN
         WRITE (*, *) rank, TRIM(label), value, tol
         ERROR STOP 'USPP distributed operator mismatch'
      END IF
   END SUBROUTINE check

! **************************************************************************************************
!> \brief ...
!> \param k ...
! **************************************************************************************************
   SUBROUTINE make_augmentation(k)
      INTEGER, INTENT(IN)                                :: k

      INTEGER, PARAMETER                                 :: nr = 1601

      INTEGER                                            :: i, j, l, pair

      pot(k)%number_of_proj = 2
      pot(k)%is_ultrasoft = .TRUE.
      ALLOCATE (pot(k)%lbeta(2), pot(k)%dion(2, 2), pot(k)%qij(2, 2))
      pot(k)%lbeta(:) = [0, 2 - k]
      pot(k)%dion = 0.0_dp
      pot(k)%dion(1, 1) = 0.7_dp
      pot(k)%dion(2, 2) = -0.3_dp
      IF (k == 2) THEN
         pot(k)%dion(1, 2) = 0.11_dp
         pot(k)%dion(2, 1) = 0.11_dp
      END IF
      pot(k)%zion = REAL(k, dp)
      pot(k)%qij = 42.0_dp
      ALLOCATE (pot(k)%r(nr), pot(k)%rab(nr), pot(k)%beta(nr, 2), pot(k)%qijl(nr, 3, 0:2*(2 - k)))
      pot(k)%r(:) = [(REAL(i - 1, dp)*0.00625_dp, i=1, nr)]
      pot(k)%rab = 0.00625_dp
      ALLOCATE (pot(k)%vlocal(nr))
      pot(k)%vlocal(1) = -2.0_dp*pot(k)%zion*SQRT(0.8_dp/pi)
      pot(k)%vlocal(2:) = -pot(k)%zion*ERF(SQRT(0.8_dp)*pot(k)%r(2:))/pot(k)%r(2:)
      pot(k)%vlocal(:) = pot(k)%vlocal + REAL(k, dp)*EXP(-(0.6_dp + 0.2_dp*k)*pot(k)%r**2)
      pot(k)%qijl = 0.0_dp
      DO j = 1, 2
         pot(k)%beta(:, j) = pot(k)%r**pot(k)%lbeta(j)*EXP(-alpha(k)*pot(k)%r**2)
         DO i = 1, j
            pair = j*(j - 1)/2 + i
            DO l = ABS(pot(k)%lbeta(i) - pot(k)%lbeta(j)), pot(k)%lbeta(i) + pot(k)%lbeta(j), 2
               pot(k)%qijl(:, pair, l) = coefficient(i, j, k)*pot(k)%r**(pot(k)%lbeta(i) + pot(k)%lbeta(j) + 2) &
                                         *EXP(-alpha(k)*pot(k)%r**2)
            END DO
         END DO
      END DO
      CALL upf_projector_init(pot(k), 1, beta(k))
      CALL uspp_augmentation_init(pot(k), beta(k), aug(k))
   END SUBROUTINE make_augmentation

! **************************************************************************************************
!> \brief ...
!> \param i ...
!> \param j ...
!> \param k ...
!> \return ...
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION coefficient(i, j, k)
      INTEGER, INTENT(IN)                                :: i, j, k

      coefficient = 0.04_dp*(i + j + k)
      IF (i /= j) coefficient = -coefficient
   END FUNCTION coefficient

! **************************************************************************************************
!> \brief ...
!> \param g ...
!> \param RESULT ...
! **************************************************************************************************
   SUBROUTINE fourier_reference(g, RESULT)
      REAL(KIND=dp), INTENT(IN)                          :: g(3)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: result(nproj, nproj)

      COMPLEX(KIND=dp)                                   :: moments(0:2, 3), phase, term
      INTEGER                                            :: a, axis, c, i, j, k, li, lj, mi, mj, &
                                                            powers(3), u, v

      RESULT = CMPLX(0.0_dp, 0.0_dp, dp)
      DO k = 1, 2
         DO axis = 1, 3
            moments(0, axis) = SQRT(pi/alpha(k))*EXP(-g(axis)**2/(4.0_dp*alpha(k)))
            moments(1, axis) = CMPLX(0.0_dp, -g(axis)/(2.0_dp*alpha(k)), dp)*moments(0, axis)
            moments(2, axis) = moments(0, axis)/(2.0_dp*alpha(k)) &
                               + CMPLX(0.0_dp, -g(axis)/(2.0_dp*alpha(k)), dp)*moments(1, axis)
         END DO
         phase = EXP(CMPLX(0.0_dp, -DOT_PRODUCT(g, positions(:, k)), dp))/volume
         DO v = 1, aug(k)%nproj
            j = beta(k)%radial_index(v); lj = pot(k)%lbeta(j); mj = beta(k)%angular_index(v) - lj*lj
            DO u = 1, aug(k)%nproj
               i = beta(k)%radial_index(u); li = pot(k)%lbeta(i); mi = beta(k)%angular_index(u) - li*li
               DO c = 1, nco(lj)
                  DO a = 1, nco(li)
                     powers = indco(:, ncoset(li - 1) + a) + indco(:, ncoset(lj - 1) + c)
                     term = PRODUCT([moments(powers(1), 1), moments(powers(2), 2), moments(powers(3), 3)])
                     RESULT(first(k) + u - 1, first(k) + v - 1) = RESULT(first(k) + u - 1, first(k) + v - 1) &
                                                 + orbtramat(li)%slm(mi, a)*orbtramat(lj)%slm(mj, c)*term*phase*coefficient(i, j, k)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE fourier_reference

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param name ...
!> \param dist ...
!> \param row_size ...
!> \param col_size ...
!> \param dense ...
!> \param sym ...
! **************************************************************************************************
   SUBROUTINE make_matrix(matrix, name, dist, row_size, col_size, dense, sym)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      TYPE(dbcsr_distribution_type), INTENT(IN)          :: dist
      INTEGER, INTENT(INOUT), POINTER                    :: row_size(:), col_size(:)
      REAL(KIND=dp), INTENT(IN)                          :: dense(:, :)
      CHARACTER, INTENT(IN)                              :: sym

      INTEGER                                            :: c, i, j, owner, r

      CALL dbcsr_create(matrix, name, dist, sym, row_size, col_size)
      r = 1
      DO i = 1, SIZE(row_size)
         c = 1
         DO j = 1, SIZE(col_size)
            IF (sym == 'S' .AND. j < i) THEN
               c = c + col_size(j)
               CYCLE
            END IF
            CALL dbcsr_get_stored_coordinates(matrix, i, j, owner)
            IF (owner == rank .AND. ANY(dense(r:r + row_size(i) - 1, c:c + col_size(j) - 1) /= 0.0_dp)) THEN
               CALL dbcsr_put_block(matrix, i, j, dense(r:r + row_size(i) - 1, c:c + col_size(j) - 1))
            END IF
            c = c + col_size(j)
         END DO
         r = r + row_size(i)
      END DO
      CALL dbcsr_finalize(matrix)
   END SUBROUTINE make_matrix

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param dense ...
! **************************************************************************************************
   SUBROUTINE gather(matrix, dense)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), INTENT(OUT)                         :: dense(:, :)

      CHARACTER                                          :: sym
      INTEGER                                            :: c, i, j, nc, nr, r
      INTEGER, POINTER                                   :: cols(:), rows(:)
      LOGICAL                                            :: transposed
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_get_info(matrix, row_blk_offset=rows, col_blk_offset=cols, matrix_type=sym)
      dense = 0.0_dp
      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, i, j, block, transposed=transposed)
         r = rows(i); c = cols(j)
         nr = SIZE(block, 1); nc = SIZE(block, 2)
         IF (transposed) THEN
            dense(r:r + nc - 1, c:c + nr - 1) = TRANSPOSE(block)
            IF (sym == 'S' .AND. i /= j) dense(c:c + nr - 1, r:r + nc - 1) = block
         ELSE
            dense(r:r + nr - 1, c:c + nc - 1) = block
            IF (sym == 'S' .AND. i /= j) dense(c:c + nc - 1, r:r + nr - 1) = TRANSPOSE(block)
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL comm%sum(dense)
   END SUBROUTINE gather
END PROGRAM uspp_projector_operators_unittest
