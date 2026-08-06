!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief GAPW reciprocal-space reconstruction and its discrete adjoint.
!> \author Matthias Krack (original density reconstruction)
! **************************************************************************************************
MODULE gapw_gspace_reconstruction
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: indco,&
                                              nco,&
                                              ncoset,&
                                              nso,&
                                              nsoset
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_types,                  ONLY: particle_type
   USE paw_basis_types,                 ONLY: get_paw_basis_info
   USE pw_grid_types,                   ONLY: HALFSPACE
   USE pw_grids,                        ONLY: get_pw_grid_info
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                              rho_atom_coeff,&
                                              rho_atom_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gapw_gspace_reconstruction'

   PUBLIC :: calculate_rhotot_elec_gspace, &
             integrate_gapw_composite_vxc_gspace

CONTAINS

! **************************************************************************************************
!> \brief  The total electronic density in reciprocal space (g-space) is
!>         calculated.
!> \param qs_env ...
!> \param auxbas_pw_pool ...
!> \param rhotot_elec_gspace ...
!> \param q_max ...
!> \param rho_hard ...
!> \param rho_soft ...
!> \param fsign ...
!> \param compute_tau reconstruct the kinetic-energy density instead of rho
!> \param rho_source density structure providing the smooth rho or tau field
!> \param allow_nonorthorhombic permit a general cell for non-XRD reconstruction
!> \date   14.03.2008 (splitted from the routine xray_diffraction_spectrum)
!> \author Matthias Krack
!> \note   This code assumes that the g-vectors are ordered (in gsq and %cc)
! **************************************************************************************************
   SUBROUTINE calculate_rhotot_elec_gspace(qs_env, auxbas_pw_pool, &
                                           rhotot_elec_gspace, q_max, rho_hard, &
                                           rho_soft, fsign, compute_tau, rho_source, &
                                           allow_nonorthorhombic)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rhotot_elec_gspace
      REAL(KIND=dp), INTENT(IN)                          :: q_max
      REAL(KIND=dp), INTENT(OUT)                         :: rho_hard, rho_soft
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: fsign
      LOGICAL, INTENT(IN), OPTIONAL                      :: compute_tau
      TYPE(qs_rho_type), OPTIONAL, POINTER               :: rho_source
      LOGICAL, INTENT(IN), OPTIONAL                      :: allow_nonorthorhombic

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_rhotot_elec_gspace'

      INTEGER :: atom, handle, iatom, ico, ico1_pgf, ico1_set, ikind, ipgf, iset, iso, iso1_pgf, &
         iso1_set, ison, ispin, jco, jco1_pgf, jco1_set, jpgf, jset, jso, jso1_pgf, jso1_set, &
         json, la, lb, maxco, maxso, na, natom, nb, ncoa, ncob, ncotot, nkind, nsatbas, nset, &
         nsoa, nsob, nsotot, nspin
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, lmax, lmin, npgf, o2nindex
      LOGICAL                                            :: my_compute_tau, orthorhombic, paw_atom, &
                                                            permit_nonorthorhombic
      REAL(KIND=dp)                                      :: alpha, eps_rho_gspace, rho_total, scale, &
                                                            volume
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: delta_cpc, pab, work, zet
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_c1d_gs_type)                               :: rho_elec_gspace
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r, source_r, tau_r
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: cpc_h, cpc_s
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      TYPE(rho_atom_type), POINTER                       :: rho_atom

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(auxbas_pw_pool))

      CALL timeset(routineN, handle)

      NULLIFY (atom_list)
      NULLIFY (atomic_kind_set)
      NULLIFY (qs_kind_set)
      NULLIFY (cell)
      NULLIFY (cpc_h)
      NULLIFY (cpc_s)
      NULLIFY (delta_cpc)
      NULLIFY (dft_control)
      NULLIFY (lmax)
      NULLIFY (lmin)
      NULLIFY (npgf)
      NULLIFY (basis_1c_set)
      NULLIFY (pab)
      NULLIFY (particle_set)
      NULLIFY (rho, rho_r, source_r, tau_r)
      NULLIFY (rho_atom)
      NULLIFY (rho_atom_set)
      NULLIFY (work)
      NULLIFY (zet)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      rho_atom_set=rho_atom_set)

      IF (PRESENT(rho_source)) THEN
         rho => rho_source
      ELSE
         CALL get_qs_env(qs_env=qs_env, rho=rho)
      END IF

      CALL qs_rho_get(rho, rho_r=rho_r, tau_r=tau_r)
      my_compute_tau = .FALSE.
      IF (PRESENT(compute_tau)) my_compute_tau = compute_tau
      IF (my_compute_tau) THEN
         CPASSERT(ASSOCIATED(tau_r))
         source_r => tau_r
      ELSE
         source_r => rho_r
      END IF
      eps_rho_gspace = dft_control%qs_control%eps_rho_gspace
      nkind = SIZE(atomic_kind_set)
      nspin = dft_control%nspins

      ! Load the soft contribution of the electronic density

      CALL auxbas_pw_pool%create_pw(pw=rho_elec_gspace)

      CALL pw_zero(rhotot_elec_gspace)

      DO ispin = 1, nspin
         CALL pw_zero(rho_elec_gspace)
         CALL pw_transfer(source_r(ispin), rho_elec_gspace)
         IF (PRESENT(fsign) .AND. (ispin == 2)) THEN
            alpha = fsign
         ELSE
            alpha = 1.0_dp
         END IF
         CALL pw_axpy(rho_elec_gspace, rhotot_elec_gspace, alpha=alpha)
      END DO

      ! Release the auxiliary PW grid for the calculation of the soft
      ! contribution

      CALL auxbas_pw_pool%give_back_pw(rho_elec_gspace)

      rho_soft = pw_integrate_function(rhotot_elec_gspace, isign=-1)

      CALL get_pw_grid_info(pw_grid=rhotot_elec_gspace%pw_grid, vol=volume, &
                            orthorhombic=orthorhombic)
      permit_nonorthorhombic = .FALSE.
      IF (PRESENT(allow_nonorthorhombic)) permit_nonorthorhombic = allow_nonorthorhombic
      IF (.NOT. orthorhombic .AND. .NOT. permit_nonorthorhombic) THEN
         CALL cp_abort(__LOCATION__, &
                       "The calculation of XRD spectra for non-orthorhombic cells is not implemented")
      END IF

      CALL pw_scale(rhotot_elec_gspace, volume)

      ! Add the hard contribution of the electronic density

      ! Each process has to loop over all PAW atoms, since the g-space grid
      ! is already distributed over all processes

      DO ikind = 1, nkind

         CALL get_atomic_kind(atomic_kind_set(ikind), &
                              atom_list=atom_list, &
                              natom=natom)

         CALL get_qs_kind(qs_kind_set(ikind), &
                          basis_set=basis_1c_set, &
                          basis_type="GAPW_1C", &
                          paw_atom=paw_atom)

         IF (.NOT. paw_atom) CYCLE ! no PAW atom: nothing to do

         CALL get_paw_basis_info(basis_1c_set, o2nindex=o2nindex, nsatbas=nsatbas)

         CALL get_gto_basis_set(gto_basis_set=basis_1c_set, &
                                lmax=lmax, &
                                lmin=lmin, &
                                maxco=maxco, &
                                maxso=maxso, &
                                npgf=npgf, &
                                nset=nset, &
                                zet=zet)

         ncotot = maxco*nset
         nsotot = maxso*nset
         CALL reallocate(delta_cpc, 1, nsatbas, 1, nsatbas)
         CALL reallocate(pab, 1, ncotot, 1, ncotot)
         CALL reallocate(work, 1, maxso, 1, maxco)

         DO iatom = 1, natom

            atom = atom_list(iatom)
            rho_atom => rho_atom_set(atom)

            CALL get_rho_atom(rho_atom=rho_atom, &
                              cpc_h=cpc_h, &
                              cpc_s=cpc_s)

            ra(:) = pbc(particle_set(atom)%r, cell)

            delta_cpc = 0.0_dp

            DO ispin = 1, nspin
               IF (PRESENT(fsign) .AND. (ispin == 2)) THEN
                  alpha = fsign
               ELSE
                  alpha = 1.0_dp
               END IF
               delta_cpc = delta_cpc + alpha*(cpc_h(ispin)%r_coef - cpc_s(ispin)%r_coef)
            END DO

            scale = 1.0_dp

            DO iset = 1, nset
               ico1_set = (iset - 1)*maxco + 1
               iso1_set = (iset - 1)*maxso + 1
               ncoa = ncoset(lmax(iset))
               nsoa = nsoset(lmax(iset))
               DO jset = 1, nset
                  jco1_set = (jset - 1)*maxco + 1
                  jso1_set = (jset - 1)*maxso + 1
                  ncob = ncoset(lmax(jset))
                  nsob = nsoset(lmax(jset))
                  DO ipgf = 1, npgf(iset)
                     ico1_pgf = ico1_set + (ipgf - 1)*ncoa
                     iso1_pgf = iso1_set + (ipgf - 1)*nsoa
                     DO jpgf = 1, npgf(jset)
                        jco1_pgf = jco1_set + (jpgf - 1)*ncob
                        jso1_pgf = jso1_set + (jpgf - 1)*nsob
                        ico = ico1_pgf + ncoset(lmin(iset) - 1)
                        iso = iso1_pgf + nsoset(lmin(iset) - 1)

                        ! Transformation spherical to Cartesian

                        DO la = lmin(iset), lmax(iset)
                           jco = jco1_pgf + ncoset(lmin(jset) - 1)
                           jso = jso1_pgf + nsoset(lmin(jset) - 1)
                           DO lb = lmin(jset), lmax(jset)
                              ison = o2nindex(iso)
                              json = o2nindex(jso)
                              CALL dgemm("N", "N", nso(la), nco(lb), nso(lb), 1.0_dp, &
                                         delta_cpc(ison:ison + nso(la) - 1, json), SIZE(delta_cpc, 1), &
                                         orbtramat(lb)%slm, nso(lb), 0.0_dp, work, &
                                         maxso)
                              CALL dgemm("T", "N", nco(la), nco(lb), nso(la), 1.0_dp, &
                                         orbtramat(la)%slm, nso(la), work, maxso, &
                                         0.0_dp, pab(ico:ico + nco(la) - 1, jco), SIZE(pab, 1))
                              jco = jco + nco(lb)
                              jso = jso + nso(lb)
                           END DO ! next lb
                           ico = ico + nco(la)
                           iso = iso + nso(la)
                        END DO ! next la

                        ! Collocate current product of primitive Cartesian functions

                        na = ico1_pgf - 1
                        nb = jco1_pgf - 1

                        CALL collocate_pgf_product_gspace( &
                           la_max=lmax(iset), &
                           zeta=zet(ipgf, iset), &
                           la_min=lmin(iset), &
                           lb_max=lmax(jset), &
                           zetb=zet(jpgf, jset), &
                           lb_min=lmin(jset), &
                           ra=ra, &
                           rab=[0.0_dp, 0.0_dp, 0.0_dp], &
                           rab2=0.0_dp, &
                           scale=scale, &
                           pab=pab, &
                           na=na, &
                           nb=nb, &
                           eps_rho_gspace=eps_rho_gspace, &
                           gsq_max=q_max*q_max, &
                           pw=rhotot_elec_gspace, &
                           compute_tau=my_compute_tau)

                     END DO ! next primitive Gaussian function "jpgf"
                  END DO ! next primitive Gaussian function "ipgf"
               END DO ! next shell set "jset"
            END DO ! next shell set "iset"
         END DO ! next atom "iatom" of atomic kind "ikind"
         DEALLOCATE (o2nindex)
      END DO ! next atomic kind "ikind"

      rho_total = pw_integrate_function(rhotot_elec_gspace, isign=-1)/volume

      rho_hard = rho_total - rho_soft

      ! Release work storage

      IF (ASSOCIATED(delta_cpc)) THEN
         DEALLOCATE (delta_cpc)
      END IF

      IF (ASSOCIATED(work)) THEN
         DEALLOCATE (work)
      END IF

      IF (ASSOCIATED(pab)) THEN
         DEALLOCATE (pab)
      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_rhotot_elec_gspace

! **************************************************************************************************
!> \brief Apply the discrete adjoint of the common-grid GAPW hard-minus-soft reconstruction.
!> \param qs_env ...
!> \param auxbas_pw_pool ...
!> \param vxc_rho derivative of the XC energy with respect to rho
!> \param vxc_tau derivative of the XC energy with respect to tau
!> \param atom_force explicit force from translating the one-center reconstruction
!> \param strain_virial derivative of the one-center reconstruction with respect to strain
!> \param one_center_contraction contraction of VXC with the hard-minus-soft fields
!> \param one_center_rho_contraction rho contribution to one_center_contraction
!> \param one_center_tau_contraction tau contribution to one_center_contraction
! **************************************************************************************************
   SUBROUTINE integrate_gapw_composite_vxc_gspace(qs_env, auxbas_pw_pool, &
                                                  vxc_rho, vxc_tau, atom_force, strain_virial, &
                                                  one_center_contraction, one_center_rho_contraction, &
                                                  one_center_tau_contraction)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), INTENT(IN)     :: vxc_rho, vxc_tau
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL                                        :: atom_force, strain_virial
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: one_center_contraction, &
                                                            one_center_rho_contraction, &
                                                            one_center_tau_contraction

      CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_gapw_composite_vxc_gspace'

      INTEGER :: atom, handle, iatom, ico, ico1_pgf, ico1_set, ikind, ipgf, iset, iso, iso1_pgf, &
         iso1_set, ison, ispin, jco, jco1_pgf, jco1_set, jpgf, jset, jso, jso1_pgf, jso1_set, &
         json, la, lb, maxco, maxso, na, natom, nb, ncoa, ncob, ncotot, nkind, nsatbas, nset, &
         nsoa, nsob, nsotot, nspin
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, lmax, lmin, npgf, o2nindex
      LOGICAL                                            :: adjoint_consistent, paw_atom
      REAL(KIND=dp)                                      :: adjoint_cart, adjoint_rho_cart, &
                                                            adjoint_spherical, adjoint_tau_cart, &
                                                            eps_rho_gspace, q_max
      REAL(KIND=dp), DIMENSION(3)                        :: center_force, ra
      REAL(KIND=dp), DIMENSION(3, 3)                     :: atom_strain_virial
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: delta_cpc, hab, hab_rho, hab_tau, pab, &
                                                            v_spherical, work, zet
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_c1d_gs_type), ALLOCATABLE, DIMENSION(:)    :: vrho_g, vtau_g
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: cpc_h, cpc_s, int_h, int_s
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      TYPE(rho_atom_type), POINTER                       :: rho_atom

      CALL timeset(routineN, handle)

      NULLIFY (atom_list, atomic_kind_set, basis_1c_set, cell, cpc_h, cpc_s, delta_cpc, &
               dft_control, hab, hab_rho, hab_tau, int_h, int_s, lmax, lmin, npgf, o2nindex, &
               pab, para_env, &
               particle_set, qs_kind_set, rho_atom, rho_atom_set, work, zet)
      NULLIFY (v_spherical)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(auxbas_pw_pool))
      CPASSERT(SIZE(vxc_rho) == SIZE(vxc_tau))
      IF (PRESENT(one_center_contraction)) one_center_contraction = 0.0_dp
      IF (PRESENT(one_center_rho_contraction)) one_center_rho_contraction = 0.0_dp
      IF (PRESENT(one_center_tau_contraction)) one_center_tau_contraction = 0.0_dp
      IF (PRESENT(atom_force)) THEN
         CPASSERT(SIZE(atom_force, 1) == 3)
      END IF
      IF (PRESENT(strain_virial)) THEN
         CPASSERT(SIZE(strain_virial, 1) == 3)
         CPASSERT(SIZE(strain_virial, 2) == 3)
      END IF

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      particle_set=particle_set, &
                      rho_atom_set=rho_atom_set)

      eps_rho_gspace = dft_control%qs_control%eps_rho_gspace
      nkind = SIZE(atomic_kind_set)
      nspin = SIZE(vxc_rho)

      ALLOCATE (vrho_g(nspin), vtau_g(nspin))
      DO ispin = 1, nspin
         CALL auxbas_pw_pool%create_pw(vrho_g(ispin))
         CALL auxbas_pw_pool%create_pw(vtau_g(ispin))
         CALL pw_transfer(vxc_rho(ispin), vrho_g(ispin))
         CALL pw_transfer(vxc_tau(ispin), vtau_g(ispin))
      END DO
      q_max = SQRT(MAXVAL(vrho_g(1)%pw_grid%gsq))

      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=basis_1c_set, &
                          basis_type="GAPW_1C", paw_atom=paw_atom)
         IF (.NOT. paw_atom) CYCLE

         CALL get_paw_basis_info(basis_1c_set, o2nindex=o2nindex, nsatbas=nsatbas)
         CALL get_gto_basis_set(gto_basis_set=basis_1c_set, &
                                lmax=lmax, lmin=lmin, maxco=maxco, maxso=maxso, &
                                npgf=npgf, nset=nset, zet=zet)
         ncotot = maxco*nset
         nsotot = maxso*nset
         CALL reallocate(hab, 1, ncotot, 1, ncotot)
         CALL reallocate(hab_rho, 1, ncotot, 1, ncotot)
         CALL reallocate(hab_tau, 1, ncotot, 1, ncotot)
         CALL reallocate(pab, 1, ncotot, 1, ncotot)
         CALL reallocate(delta_cpc, 1, nsatbas, 1, nsatbas)
         CALL reallocate(v_spherical, 1, nsotot, 1, nsotot)
         CALL reallocate(work, 1, maxso, 1, maxco)
         pab = 0.0_dp

         DO iatom = 1, natom
            atom = atom_list(iatom)
            rho_atom => rho_atom_set(atom)
            CALL get_rho_atom(rho_atom=rho_atom, cpc_h=cpc_h, cpc_s=cpc_s, &
                              ga_Vlocal_gb_h=int_h, ga_Vlocal_gb_s=int_s)
            ra(:) = pbc(particle_set(atom)%r, cell)

            DO ispin = 1, nspin
               hab = 0.0_dp
               hab_rho = 0.0_dp
               hab_tau = 0.0_dp
               pab = 0.0_dp
               center_force = 0.0_dp
               atom_strain_virial = 0.0_dp
               delta_cpc = cpc_h(ispin)%r_coef - cpc_s(ispin)%r_coef
               DO iset = 1, nset
                  ico1_set = (iset - 1)*maxco + 1
                  iso1_set = (iset - 1)*maxso + 1
                  ncoa = ncoset(lmax(iset))
                  nsoa = nsoset(lmax(iset))
                  DO jset = 1, nset
                     jco1_set = (jset - 1)*maxco + 1
                     jso1_set = (jset - 1)*maxso + 1
                     ncob = ncoset(lmax(jset))
                     nsob = nsoset(lmax(jset))
                     DO ipgf = 1, npgf(iset)
                        ico1_pgf = ico1_set + (ipgf - 1)*ncoa
                        iso1_pgf = iso1_set + (ipgf - 1)*nsoa
                        DO jpgf = 1, npgf(jset)
                           jco1_pgf = jco1_set + (jpgf - 1)*ncob
                           jso1_pgf = jso1_set + (jpgf - 1)*nsob
                           ico = ico1_pgf + ncoset(lmin(iset) - 1)
                           iso = iso1_pgf + nsoset(lmin(iset) - 1)
                           DO la = lmin(iset), lmax(iset)
                              jco = jco1_pgf + ncoset(lmin(jset) - 1)
                              jso = jso1_pgf + nsoset(lmin(jset) - 1)
                              DO lb = lmin(jset), lmax(jset)
                                 ison = o2nindex(iso)
                                 json = o2nindex(jso)
                                 CALL dgemm("N", "N", nso(la), nco(lb), nso(lb), 1.0_dp, &
                                            delta_cpc(ison:ison + nso(la) - 1, json), &
                                            SIZE(delta_cpc, 1), orbtramat(lb)%slm, nso(lb), &
                                            0.0_dp, work, maxso)
                                 CALL dgemm("T", "N", nco(la), nco(lb), nso(la), 1.0_dp, &
                                            orbtramat(la)%slm, nso(la), work, maxso, 0.0_dp, &
                                            pab(ico:ico + nco(la) - 1, jco), SIZE(pab, 1))
                                 jco = jco + nco(lb)
                                 jso = jso + nso(lb)
                              END DO
                              ico = ico + nco(la)
                              iso = iso + nso(la)
                           END DO
                           na = ico1_pgf - 1
                           nb = jco1_pgf - 1
                           IF (PRESENT(strain_virial)) THEN
                              CALL collocate_pgf_product_gspace( &
                                 la_max=lmax(iset), zeta=zet(ipgf, iset), la_min=lmin(iset), &
                                 lb_max=lmax(jset), zetb=zet(jpgf, jset), lb_min=lmin(jset), &
                                 ra=ra, rab=[0.0_dp, 0.0_dp, 0.0_dp], rab2=0.0_dp, scale=1.0_dp, &
                                 pab=pab, na=na, nb=nb, eps_rho_gspace=eps_rho_gspace, &
                                 gsq_max=q_max*q_max, pw=vrho_g(ispin), &
                                 potential_gspace=vrho_g(ispin), hab=hab_rho, &
                                 center_force=center_force, strain_virial=atom_strain_virial)
                              CALL collocate_pgf_product_gspace( &
                                 la_max=lmax(iset), zeta=zet(ipgf, iset), la_min=lmin(iset), &
                                 lb_max=lmax(jset), zetb=zet(jpgf, jset), lb_min=lmin(jset), &
                                 ra=ra, rab=[0.0_dp, 0.0_dp, 0.0_dp], rab2=0.0_dp, scale=1.0_dp, &
                                 pab=pab, na=na, nb=nb, eps_rho_gspace=eps_rho_gspace, &
                                 gsq_max=q_max*q_max, pw=vtau_g(ispin), compute_tau=.TRUE., &
                                 potential_gspace=vtau_g(ispin), hab=hab_tau, &
                                 center_force=center_force, strain_virial=atom_strain_virial)
                           ELSE IF (PRESENT(atom_force)) THEN
                              CALL collocate_pgf_product_gspace( &
                                 la_max=lmax(iset), zeta=zet(ipgf, iset), la_min=lmin(iset), &
                                 lb_max=lmax(jset), zetb=zet(jpgf, jset), lb_min=lmin(jset), &
                                 ra=ra, rab=[0.0_dp, 0.0_dp, 0.0_dp], rab2=0.0_dp, scale=1.0_dp, &
                                 pab=pab, na=na, nb=nb, eps_rho_gspace=eps_rho_gspace, &
                                 gsq_max=q_max*q_max, pw=vrho_g(ispin), &
                                 potential_gspace=vrho_g(ispin), hab=hab_rho, &
                                 center_force=center_force)
                              CALL collocate_pgf_product_gspace( &
                                 la_max=lmax(iset), zeta=zet(ipgf, iset), la_min=lmin(iset), &
                                 lb_max=lmax(jset), zetb=zet(jpgf, jset), lb_min=lmin(jset), &
                                 ra=ra, rab=[0.0_dp, 0.0_dp, 0.0_dp], rab2=0.0_dp, scale=1.0_dp, &
                                 pab=pab, na=na, nb=nb, eps_rho_gspace=eps_rho_gspace, &
                                 gsq_max=q_max*q_max, pw=vtau_g(ispin), compute_tau=.TRUE., &
                                 potential_gspace=vtau_g(ispin), hab=hab_tau, &
                                 center_force=center_force)
                           ELSE
                              CALL collocate_pgf_product_gspace( &
                                 la_max=lmax(iset), zeta=zet(ipgf, iset), la_min=lmin(iset), &
                                 lb_max=lmax(jset), zetb=zet(jpgf, jset), lb_min=lmin(jset), &
                                 ra=ra, rab=[0.0_dp, 0.0_dp, 0.0_dp], rab2=0.0_dp, scale=1.0_dp, &
                                 pab=pab, na=na, nb=nb, eps_rho_gspace=eps_rho_gspace, &
                                 gsq_max=q_max*q_max, pw=vrho_g(ispin), &
                                 potential_gspace=vrho_g(ispin), hab=hab_rho)
                              CALL collocate_pgf_product_gspace( &
                                 la_max=lmax(iset), zeta=zet(ipgf, iset), la_min=lmin(iset), &
                                 lb_max=lmax(jset), zetb=zet(jpgf, jset), lb_min=lmin(jset), &
                                 ra=ra, rab=[0.0_dp, 0.0_dp, 0.0_dp], rab2=0.0_dp, scale=1.0_dp, &
                                 pab=pab, na=na, nb=nb, eps_rho_gspace=eps_rho_gspace, &
                                 gsq_max=q_max*q_max, pw=vtau_g(ispin), compute_tau=.TRUE., &
                                 potential_gspace=vtau_g(ispin), hab=hab_tau)
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
               CALL para_env%sum(hab_rho)
               CALL para_env%sum(hab_tau)
               ! Atom forces and virials remain rank-local here.  Quickstep performs the
               ! corresponding global reductions after all force components are assembled.
               hab = hab_rho + hab_tau
               adjoint_rho_cart = SUM(pab*hab_rho)
               adjoint_tau_cart = SUM(pab*hab_tau)
               adjoint_cart = adjoint_rho_cart + adjoint_tau_cart
               IF (PRESENT(one_center_contraction)) THEN
                  one_center_contraction = one_center_contraction + adjoint_cart
               END IF
               IF (PRESENT(one_center_rho_contraction)) THEN
                  one_center_rho_contraction = one_center_rho_contraction + adjoint_rho_cart
               END IF
               IF (PRESENT(one_center_tau_contraction)) THEN
                  one_center_tau_contraction = one_center_tau_contraction + adjoint_tau_cart
               END IF

               IF (ASSOCIATED(int_h) .AND. ASSOCIATED(int_s)) THEN
                  IF (ASSOCIATED(int_h(ispin)%r_coef) .AND. ASSOCIATED(int_s(ispin)%r_coef)) THEN
                     v_spherical = 0.0_dp
                     DO iset = 1, nset
                        ico1_set = (iset - 1)*maxco + 1
                        iso1_set = (iset - 1)*maxso + 1
                        ncoa = ncoset(lmax(iset))
                        nsoa = nsoset(lmax(iset))
                        DO jset = 1, nset
                           jco1_set = (jset - 1)*maxco + 1
                           jso1_set = (jset - 1)*maxso + 1
                           ncob = ncoset(lmax(jset))
                           nsob = nsoset(lmax(jset))
                           DO ipgf = 1, npgf(iset)
                              ico1_pgf = ico1_set + (ipgf - 1)*ncoa
                              iso1_pgf = iso1_set + (ipgf - 1)*nsoa
                              DO jpgf = 1, npgf(jset)
                                 jco1_pgf = jco1_set + (jpgf - 1)*ncob
                                 jso1_pgf = jso1_set + (jpgf - 1)*nsob
                                 ico = ico1_pgf + ncoset(lmin(iset) - 1)
                                 iso = iso1_pgf + nsoset(lmin(iset) - 1)
                                 DO la = lmin(iset), lmax(iset)
                                    jco = jco1_pgf + ncoset(lmin(jset) - 1)
                                    jso = jso1_pgf + nsoset(lmin(jset) - 1)
                                    DO lb = lmin(jset), lmax(jset)
                                       CALL dgemm("N", "N", nso(la), nco(lb), nco(la), 1.0_dp, &
                                                  orbtramat(la)%slm, nso(la), &
                                                  hab(ico:ico + nco(la) - 1, jco), &
                                                  SIZE(hab, 1), 0.0_dp, work, maxso)
                                       CALL dgemm("N", "T", nso(la), nso(lb), nco(lb), 1.0_dp, &
                                                  work, maxso, orbtramat(lb)%slm, nso(lb), 1.0_dp, &
                                                  v_spherical(iso:iso + nso(la) - 1, jso), nsotot)
                                       jco = jco + nco(lb)
                                       jso = jso + nso(lb)
                                    END DO
                                    ico = ico + nco(la)
                                    iso = iso + nso(la)
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                     int_h(ispin)%r_coef = int_h(ispin)%r_coef + v_spherical
                     int_s(ispin)%r_coef = int_s(ispin)%r_coef + v_spherical
                     adjoint_spherical = 0.0_dp
                     DO jso = 1, nsotot
                        json = o2nindex(jso)
                        IF (json == 0) CYCLE
                        DO iso = 1, nsotot
                           ison = o2nindex(iso)
                           IF (ison == 0) CYCLE
                           adjoint_spherical = adjoint_spherical + &
                                               delta_cpc(ison, json)*v_spherical(iso, jso)
                        END DO
                     END DO
                     adjoint_consistent = ABS(adjoint_cart - adjoint_spherical) <= &
                        1.0E-9_dp*MAX(1.0_dp, ABS(adjoint_cart))
                     IF (.NOT. adjoint_consistent) THEN
                        WRITE (*, '(A,3I6,3ES24.15)') &
                           " SKALA_GPW| Composite adjoint mismatch", ikind, atom, ispin, &
                           adjoint_cart, adjoint_spherical, adjoint_cart - adjoint_spherical
                     END IF
                     CPASSERT(adjoint_consistent)
                  END IF
               END IF
               IF (PRESENT(atom_force)) atom_force(:, atom) = atom_force(:, atom) + center_force
               IF (PRESENT(strain_virial)) strain_virial = strain_virial + atom_strain_virial
            END DO
         END DO

         DEALLOCATE (delta_cpc, hab, hab_rho, hab_tau, pab, v_spherical, work, o2nindex)
      END DO

      DO ispin = 1, nspin
         CALL auxbas_pw_pool%give_back_pw(vrho_g(ispin))
         CALL auxbas_pw_pool%give_back_pw(vtau_g(ispin))
      END DO
      DEALLOCATE (vrho_g, vtau_g)

      CALL timestop(handle)

   END SUBROUTINE integrate_gapw_composite_vxc_gspace

! **************************************************************************************************
!> \brief low level collocation of primitive gaussian functions in g-space
!> \param la_max ...
!> \param zeta ...
!> \param la_min ...
!> \param lb_max ...
!> \param zetb ...
!> \param lb_min ...
!> \param ra ...
!> \param rab ...
!> \param rab2 ...
!> \param scale ...
!> \param pab ...
!> \param na ...
!> \param nb ...
!> \param eps_rho_gspace ...
!> \param gsq_max ...
!> \param pw ...
!> \param compute_tau apply the symmetric kinetic-density DADB transformation
!> \param potential_gspace integrate each primitive product against this potential
!> \param hab Cartesian primitive-product potential integrals
!> \param center_force force from translating the collocated product density
!> \param strain_virial derivative of the product density with respect to strain
! **************************************************************************************************
   SUBROUTINE collocate_pgf_product_gspace(la_max, zeta, la_min, &
                                           lb_max, zetb, lb_min, &
                                           ra, rab, rab2, scale, pab, na, nb, &
                                           eps_rho_gspace, gsq_max, pw, compute_tau, &
                                           potential_gspace, hab, center_force, strain_virial)

      INTEGER, INTENT(IN)                                :: la_max
      REAL(dp), INTENT(IN)                               :: zeta
      INTEGER, INTENT(IN)                                :: la_min, lb_max
      REAL(dp), INTENT(IN)                               :: zetb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rab
      REAL(dp), INTENT(IN)                               :: rab2, scale
      REAL(dp), DIMENSION(:, :), POINTER                 :: pab
      INTEGER, INTENT(IN)                                :: na, nb
      REAL(dp), INTENT(IN)                               :: eps_rho_gspace, gsq_max
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: pw
      LOGICAL, INTENT(IN), OPTIONAL                      :: compute_tau
      TYPE(pw_c1d_gs_type), INTENT(IN), OPTIONAL         :: potential_gspace
      REAL(dp), DIMENSION(:, :), INTENT(INOUT), OPTIONAL :: hab
      REAL(dp), DIMENSION(3), INTENT(INOUT), OPTIONAL    :: center_force
      REAL(dp), DIMENSION(3, 3), INTENT(INOUT), OPTIONAL :: strain_virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'collocate_pgf_product_gspace'

      COMPLEX(dp)                                        :: plane_product, rag, rbg, shifted_product, &
                                                            tau_product
      COMPLEX(dp), DIMENSION(3)                          :: plane_gradient
      COMPLEX(dp), DIMENSION(:, :, :), ALLOCATABLE       :: cubeaxis
      INTEGER                                            :: ax, ay, az, bx, by, bz, handle, i, ico, &
                                                            idir, jco, jdir, la, la_rec_max, lb, &
                                                            lb_grid, lb_rec_max, ub_grid
      LOGICAL                                            :: integrate_potential, integrate_strain, &
                                                            my_compute_tau
      REAL(dp)                                           :: f, fa, fb, gcoord, pij, potential_term, &
                                                            prefactor, reciprocal_weight, rzetp, &
                                                            twozetp, zetp
      REAL(dp), DIMENSION(3)                             :: fap, fbp, force_product, gvec, rap, rbp, rp

      CALL timeset(routineN, handle)

      my_compute_tau = .FALSE.
      IF (PRESENT(compute_tau)) my_compute_tau = compute_tau
      integrate_potential = PRESENT(potential_gspace) .OR. PRESENT(hab)
      integrate_strain = PRESENT(strain_virial)
      CPASSERT(PRESENT(potential_gspace) .EQV. PRESENT(hab))
      CPASSERT(.NOT. integrate_strain .OR. integrate_potential)
      IF (integrate_potential) THEN
         CPASSERT(ASSOCIATED(pw%pw_grid, potential_gspace%pw_grid))
      END IF

      la_rec_max = la_max
      lb_rec_max = lb_max
      IF (my_compute_tau) THEN
         la_rec_max = la_rec_max + 1
         lb_rec_max = lb_rec_max + 1
      END IF
      IF (integrate_strain) la_rec_max = la_rec_max + 1

      zetp = zeta + zetb
      rzetp = 1.0_dp/zetp
      f = zetb*rzetp
      rap(:) = f*rab(:)
      rbp(:) = rap(:) - rab(:)
      rp(:) = ra(:) + rap(:)
      twozetp = 2.0_dp*zetp
      fap(:) = twozetp*rap(:)
      fbp(:) = twozetp*rbp(:)
      prefactor = scale*SQRT((pi*rzetp)**3)*EXP(-zeta*f*rab2)

      ALLOCATE (cubeaxis(3, 0:la_rec_max, 0:lb_rec_max))
      lb_grid = LBOUND(pw%array, 1)
      ub_grid = UBOUND(pw%array, 1)

      DO i = lb_grid, ub_grid
         IF (pw%pw_grid%gsq(i) > gsq_max) CYCLE
         gvec = pw%pw_grid%g(:, i)
         CALL build_cubeaxis(gvec)
         reciprocal_weight = 1.0_dp
         IF (pw%pw_grid%grid_span == HALFSPACE .AND. pw%pw_grid%gsq(i) > 0.0_dp) THEN
            reciprocal_weight = 2.0_dp
         END IF

         DO ico = ncoset(la_min - 1) + 1, ncoset(la_max)
            ax = indco(1, ico)
            ay = indco(2, ico)
            az = indco(3, ico)

            DO jco = ncoset(lb_min - 1) + 1, ncoset(lb_max)
               bx = indco(1, jco)
               by = indco(2, jco)
               bz = indco(3, jco)
               pij = prefactor*pab(na + ico, nb + jco)
               IF (.NOT. integrate_potential .AND. ABS(pij) < eps_rho_gspace) CYCLE

               IF (.NOT. my_compute_tau) THEN
                  plane_product = cubeaxis(1, ax, bx)*cubeaxis(2, ay, by)*cubeaxis(3, az, bz)
               ELSE
                  ! Equivalent to GRID_FUNC_DADB: 0.5 * grad(phi_a).grad(phi_b).
                  tau_product = 2.0_dp*zeta*zetb* &
                                cubeaxis(1, ax + 1, bx + 1)* &
                                cubeaxis(2, ay, by)*cubeaxis(3, az, bz)
                  IF (ax > 0) THEN
                     tau_product = tau_product - REAL(ax, dp)*zetb* &
                                   cubeaxis(1, ax - 1, bx + 1)* &
                                   cubeaxis(2, ay, by)*cubeaxis(3, az, bz)
                  END IF
                  IF (bx > 0) THEN
                     tau_product = tau_product - zeta*REAL(bx, dp)* &
                                   cubeaxis(1, ax + 1, bx - 1)* &
                                   cubeaxis(2, ay, by)*cubeaxis(3, az, bz)
                  END IF
                  IF (ax > 0 .AND. bx > 0) THEN
                     tau_product = tau_product + 0.5_dp*REAL(ax*bx, dp)* &
                                   cubeaxis(1, ax - 1, bx - 1)* &
                                   cubeaxis(2, ay, by)*cubeaxis(3, az, bz)
                  END IF

                  tau_product = tau_product + 2.0_dp*zeta*zetb* &
                                cubeaxis(1, ax, bx)* &
                                cubeaxis(2, ay + 1, by + 1)*cubeaxis(3, az, bz)
                  IF (ay > 0) THEN
                     tau_product = tau_product - REAL(ay, dp)*zetb* &
                                   cubeaxis(1, ax, bx)* &
                                   cubeaxis(2, ay - 1, by + 1)*cubeaxis(3, az, bz)
                  END IF
                  IF (by > 0) THEN
                     tau_product = tau_product - zeta*REAL(by, dp)* &
                                   cubeaxis(1, ax, bx)* &
                                   cubeaxis(2, ay + 1, by - 1)*cubeaxis(3, az, bz)
                  END IF
                  IF (ay > 0 .AND. by > 0) THEN
                     tau_product = tau_product + 0.5_dp*REAL(ay*by, dp)* &
                                   cubeaxis(1, ax, bx)* &
                                   cubeaxis(2, ay - 1, by - 1)*cubeaxis(3, az, bz)
                  END IF

                  tau_product = tau_product + 2.0_dp*zeta*zetb* &
                                cubeaxis(1, ax, bx)*cubeaxis(2, ay, by)* &
                                cubeaxis(3, az + 1, bz + 1)
                  IF (az > 0) THEN
                     tau_product = tau_product - REAL(az, dp)*zetb* &
                                   cubeaxis(1, ax, bx)*cubeaxis(2, ay, by)* &
                                   cubeaxis(3, az - 1, bz + 1)
                  END IF
                  IF (bz > 0) THEN
                     tau_product = tau_product - zeta*REAL(bz, dp)* &
                                   cubeaxis(1, ax, bx)*cubeaxis(2, ay, by)* &
                                   cubeaxis(3, az + 1, bz - 1)
                  END IF
                  IF (az > 0 .AND. bz > 0) THEN
                     tau_product = tau_product + 0.5_dp*REAL(az*bz, dp)* &
                                   cubeaxis(1, ax, bx)*cubeaxis(2, ay, by)* &
                                   cubeaxis(3, az - 1, bz - 1)
                  END IF
                  plane_product = tau_product
               END IF

               IF (integrate_strain) THEN
                  IF (.NOT. my_compute_tau) THEN
                     shifted_product = cubeaxis(1, ax + 1, bx)* &
                                       cubeaxis(2, ay, by)*cubeaxis(3, az, bz)
                     plane_gradient(1) = CMPLX(0.0_dp, -1.0_dp, KIND=dp)* &
                                         (shifted_product + ra(1)*plane_product)
                     shifted_product = cubeaxis(1, ax, bx)* &
                                       cubeaxis(2, ay + 1, by)*cubeaxis(3, az, bz)
                     plane_gradient(2) = CMPLX(0.0_dp, -1.0_dp, KIND=dp)* &
                                         (shifted_product + ra(2)*plane_product)
                     shifted_product = cubeaxis(1, ax, bx)* &
                                       cubeaxis(2, ay, by)*cubeaxis(3, az + 1, bz)
                     plane_gradient(3) = CMPLX(0.0_dp, -1.0_dp, KIND=dp)* &
                                         (shifted_product + ra(3)*plane_product)
                  ELSE
                     DO idir = 1, 3
                        shifted_product = tau_plane_product_shifted(idir)
                        plane_gradient(idir) = CMPLX(0.0_dp, -1.0_dp, KIND=dp)* &
                                               (shifted_product + ra(idir)*plane_product)
                     END DO
                  END IF
               END IF

               IF (integrate_potential) THEN
                  potential_term = reciprocal_weight*prefactor* &
                                   REAL(CONJG(potential_gspace%array(i))*plane_product, KIND=dp)
                  hab(na + ico, nb + jco) = hab(na + ico, nb + jco) + potential_term
                  IF (PRESENT(center_force) .OR. integrate_strain) THEN
                     force_product = reciprocal_weight*pij* &
                                     REAL(CONJG(potential_gspace%array(i))* &
                                          CMPLX(0.0_dp, -gvec, KIND=dp)*plane_product, KIND=dp)
                     IF (PRESENT(center_force)) center_force = center_force + force_product
                  END IF
                  IF (integrate_strain) THEN
                     DO idir = 1, 3
                        DO jdir = 1, 3
                           strain_virial(idir, jdir) = strain_virial(idir, jdir) - &
                                                       reciprocal_weight*pij*gvec(idir)* &
                                                       REAL(CONJG(potential_gspace%array(i))*plane_gradient(jdir), KIND=dp) + &
                                                       force_product(idir)*ra(jdir)
                        END DO
                        strain_virial(idir, idir) = strain_virial(idir, idir) - &
                                                    pab(na + ico, nb + jco)*potential_term
                     END DO
                  END IF
               ELSE
                  pw%array(i) = pw%array(i) + pij*plane_product
               END IF
            END DO
         END DO
      END DO

      DEALLOCATE (cubeaxis)
      CALL timestop(handle)

   CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param gvector ...
! **************************************************************************************************
      SUBROUTINE build_cubeaxis(gvector)
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: gvector

         cubeaxis = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO idir = 1, 3
            gcoord = gvector(idir)
            cubeaxis(idir, 0, 0) = &
               EXP(CMPLX(-0.25_dp*rzetp*gcoord*gcoord, -rp(idir)*gcoord, KIND=dp))

            IF (la_rec_max > 0) THEN
               rag = CMPLX(fap(idir), -gcoord, KIND=dp)
               cubeaxis(idir, 1, 0) = rag*cubeaxis(idir, 0, 0)
               DO la = 2, la_rec_max
                  fa = REAL(la - 1, dp)*twozetp
                  cubeaxis(idir, la, 0) = rag*cubeaxis(idir, la - 1, 0) + &
                                          fa*cubeaxis(idir, la - 2, 0)
               END DO
               IF (lb_rec_max > 0) THEN
                  rbg = CMPLX(fbp(idir), -gcoord, KIND=dp)
                  fa = twozetp
                  cubeaxis(idir, 0, 1) = rbg*cubeaxis(idir, 0, 0)
                  cubeaxis(idir, 1, 1) = rbg*cubeaxis(idir, 1, 0) + &
                                         fa*cubeaxis(idir, 0, 0)
                  DO lb = 2, lb_rec_max
                     fb = REAL(lb - 1, dp)*twozetp
                     cubeaxis(idir, 0, lb) = rbg*cubeaxis(idir, 0, lb - 1) + &
                                             fb*cubeaxis(idir, 0, lb - 2)
                     cubeaxis(idir, 1, lb) = rbg*cubeaxis(idir, 1, lb - 1) + &
                                             fb*cubeaxis(idir, 1, lb - 2) + &
                                             fa*cubeaxis(idir, 0, lb - 1)
                  END DO
                  DO la = 2, la_rec_max
                     fa = REAL(la, dp)*twozetp
                     cubeaxis(idir, la, 1) = rbg*cubeaxis(idir, la, 0) + &
                                             fa*cubeaxis(idir, la - 1, 0)
                     DO lb = 2, lb_rec_max
                        fb = REAL(lb - 1, dp)*twozetp
                        cubeaxis(idir, la, lb) = rbg*cubeaxis(idir, la, lb - 1) + &
                                                 fb*cubeaxis(idir, la, lb - 2) + &
                                                 fa*cubeaxis(idir, la - 1, lb - 1)
                     END DO
                  END DO
               END IF
            ELSE IF (lb_rec_max > 0) THEN
               rbg = CMPLX(fbp(idir), -gcoord, KIND=dp)
               cubeaxis(idir, 0, 1) = rbg*cubeaxis(idir, 0, 0)
               DO lb = 2, lb_rec_max
                  fb = REAL(lb - 1, dp)*twozetp
                  cubeaxis(idir, 0, lb) = rbg*cubeaxis(idir, 0, lb - 1) + &
                                          fb*cubeaxis(idir, 0, lb - 2)
               END DO
            END IF

            DO la = 0, la_rec_max
               DO lb = 0, lb_rec_max
                  IF (la + lb == 0) CYCLE
                  fa = (1.0_dp/twozetp)**(la + lb)
                  cubeaxis(idir, la, lb) = fa*cubeaxis(idir, la, lb)
               END DO
            END DO
         END DO
      END SUBROUTINE build_cubeaxis

! **************************************************************************************************
!> \brief ...
!> \param shift_axis ...
!> \return ...
! **************************************************************************************************
      FUNCTION tau_plane_product_shifted(shift_axis) RESULT(value)
      INTEGER, INTENT(IN)                                :: shift_axis
      COMPLEX(dp)                                        :: value

      INTEGER                                            :: sx, sy, sz

         sx = 0
         sy = 0
         sz = 0
         IF (shift_axis == 1) sx = 1
         IF (shift_axis == 2) sy = 1
         IF (shift_axis == 3) sz = 1

         value = 2.0_dp*zeta*zetb* &
                 cubeaxis(1, ax + 1 + sx, bx + 1)* &
                 cubeaxis(2, ay + sy, by)*cubeaxis(3, az + sz, bz)
         IF (ax > 0) THEN
            value = value - REAL(ax, dp)*zetb* &
                    cubeaxis(1, ax - 1 + sx, bx + 1)* &
                    cubeaxis(2, ay + sy, by)*cubeaxis(3, az + sz, bz)
         END IF
         IF (bx > 0) THEN
            value = value - zeta*REAL(bx, dp)* &
                    cubeaxis(1, ax + 1 + sx, bx - 1)* &
                    cubeaxis(2, ay + sy, by)*cubeaxis(3, az + sz, bz)
         END IF
         IF (ax > 0 .AND. bx > 0) THEN
            value = value + 0.5_dp*REAL(ax*bx, dp)* &
                    cubeaxis(1, ax - 1 + sx, bx - 1)* &
                    cubeaxis(2, ay + sy, by)*cubeaxis(3, az + sz, bz)
         END IF

         value = value + 2.0_dp*zeta*zetb* &
                 cubeaxis(1, ax + sx, bx)* &
                 cubeaxis(2, ay + 1 + sy, by + 1)*cubeaxis(3, az + sz, bz)
         IF (ay > 0) THEN
            value = value - REAL(ay, dp)*zetb* &
                    cubeaxis(1, ax + sx, bx)* &
                    cubeaxis(2, ay - 1 + sy, by + 1)*cubeaxis(3, az + sz, bz)
         END IF
         IF (by > 0) THEN
            value = value - zeta*REAL(by, dp)* &
                    cubeaxis(1, ax + sx, bx)* &
                    cubeaxis(2, ay + 1 + sy, by - 1)*cubeaxis(3, az + sz, bz)
         END IF
         IF (ay > 0 .AND. by > 0) THEN
            value = value + 0.5_dp*REAL(ay*by, dp)* &
                    cubeaxis(1, ax + sx, bx)* &
                    cubeaxis(2, ay - 1 + sy, by - 1)*cubeaxis(3, az + sz, bz)
         END IF

         value = value + 2.0_dp*zeta*zetb* &
                 cubeaxis(1, ax + sx, bx)*cubeaxis(2, ay + sy, by)* &
                 cubeaxis(3, az + 1 + sz, bz + 1)
         IF (az > 0) THEN
            value = value - REAL(az, dp)*zetb* &
                    cubeaxis(1, ax + sx, bx)*cubeaxis(2, ay + sy, by)* &
                    cubeaxis(3, az - 1 + sz, bz + 1)
         END IF
         IF (bz > 0) THEN
            value = value - zeta*REAL(bz, dp)* &
                    cubeaxis(1, ax + sx, bx)*cubeaxis(2, ay + sy, by)* &
                    cubeaxis(3, az + 1 + sz, bz - 1)
         END IF
         IF (az > 0 .AND. bz > 0) THEN
            value = value + 0.5_dp*REAL(az*bz, dp)* &
                    cubeaxis(1, ax + sx, bx)*cubeaxis(2, ay + sy, by)* &
                    cubeaxis(3, az - 1 + sz, bz - 1)
         END IF
      END FUNCTION tau_plane_product_shifted

   END SUBROUTINE collocate_pgf_product_gspace

! **************************************************************************************************
!> \brief Buffered reference retained while validating the streaming implementation.
!> \param la_max ...
!> \param zeta ...
!> \param la_min ...
!> \param lb_max ...
!> \param zetb ...
!> \param lb_min ...
!> \param ra ...
!> \param rab ...
!> \param rab2 ...
!> \param scale ...
!> \param pab ...
!> \param na ...
!> \param nb ...
!> \param eps_rho_gspace ...
!> \param gsq_max ...
!> \param pw ...
!> \param compute_tau ...
!> \param potential_gspace ...
!> \param hab ...
!> \param center_force ...
!> \param strain_virial ...
! **************************************************************************************************
   SUBROUTINE collocate_pgf_product_gspace_buffered(la_max, zeta, la_min, &
                                                    lb_max, zetb, lb_min, &
                                                    ra, rab, rab2, scale, pab, na, nb, &
                                                    eps_rho_gspace, gsq_max, pw, compute_tau, &
                                                    potential_gspace, hab, center_force, strain_virial)

      ! NOTE: this routine is much slower than the real-space version of collocate_pgf_product

      INTEGER, INTENT(IN)                                :: la_max
      REAL(dp), INTENT(IN)                               :: zeta
      INTEGER, INTENT(IN)                                :: la_min, lb_max
      REAL(dp), INTENT(IN)                               :: zetb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rab
      REAL(dp), INTENT(IN)                               :: rab2, scale
      REAL(dp), DIMENSION(:, :), POINTER                 :: pab
      INTEGER, INTENT(IN)                                :: na, nb
      REAL(dp), INTENT(IN)                               :: eps_rho_gspace, gsq_max
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: pw
      LOGICAL, INTENT(IN), OPTIONAL                      :: compute_tau
      TYPE(pw_c1d_gs_type), INTENT(IN), OPTIONAL         :: potential_gspace
      REAL(dp), DIMENSION(:, :), INTENT(INOUT), OPTIONAL :: hab
      REAL(dp), DIMENSION(3), INTENT(INOUT), OPTIONAL    :: center_force
      REAL(dp), DIMENSION(3, 3), INTENT(INOUT), OPTIONAL :: strain_virial

      CHARACTER(LEN=*), PARAMETER :: routineN = 'collocate_pgf_product_gspace'

      COMPLEX(dp)                                        :: plane_product, rag, rbg, shifted_product, &
                                                            tau_product
      COMPLEX(dp), DIMENSION(3)                          :: plane_gradient
      COMPLEX(dp), DIMENSION(:, :, :, :), POINTER        :: cubeaxis
      INTEGER                                            :: ax, ay, az, bx, by, bz, handle, i, ico, &
                                                            idir, jco, jdir, la, la_rec_max, lb, &
                                                            lb_grid, lb_rec_max, ub_grid
      LOGICAL                                            :: integrate_potential, integrate_strain, my_compute_tau
      REAL(dp)                                           :: f, fa, fb, g0_contribution, gcoord, pij, &
                                                            potential_product, prefactor, rzetp, &
                                                            twozetp, zetp
      REAL(dp), DIMENSION(3)                             :: force_product, gvec
      REAL(dp), DIMENSION(3)                             :: fap, fbp, rap, rbp, rp
      REAL(dp), DIMENSION(3, 3)                          :: strain_product

      CALL timeset(routineN, handle)

      my_compute_tau = .FALSE.
      IF (PRESENT(compute_tau)) my_compute_tau = compute_tau
      integrate_potential = PRESENT(potential_gspace) .OR. PRESENT(hab)
      integrate_strain = PRESENT(strain_virial)
      CPASSERT(PRESENT(potential_gspace) .EQV. PRESENT(hab))
      CPASSERT(.NOT. integrate_strain .OR. integrate_potential)
      IF (integrate_potential) THEN
         CPASSERT(ASSOCIATED(pw%pw_grid, potential_gspace%pw_grid))
      END IF
      la_rec_max = la_max
      lb_rec_max = lb_max
      IF (my_compute_tau) THEN
         la_rec_max = la_rec_max + 1
         lb_rec_max = lb_rec_max + 1
      END IF
      IF (integrate_strain) la_rec_max = la_rec_max + 1

      zetp = zeta + zetb
      rzetp = 1.0_dp/zetp
      f = zetb*rzetp
      rap(:) = f*rab(:)
      rbp(:) = rap(:) - rab(:)
      rp(:) = ra(:) + rap(:)
      twozetp = 2.0_dp*zetp
      fap(:) = twozetp*rap(:)
      fbp(:) = twozetp*rbp(:)

      prefactor = scale*SQRT((pi*rzetp)**3)*EXP(-zeta*f*rab2)
      lb_grid = LBOUND(pw%array, 1)
      ub_grid = UBOUND(pw%array, 1)
      NULLIFY (cubeaxis)
      CALL reallocate(cubeaxis, lb_grid, ub_grid, 1, 3, 0, la_rec_max, 0, lb_rec_max)

      DO idir = 1, 3
         DO i = lb_grid, ub_grid
            gcoord = pw%pw_grid%g(idir, i)
            cubeaxis(i, idir, 0, 0) = &
               EXP(CMPLX(-0.25_dp*rzetp*gcoord*gcoord, -rp(idir)*gcoord, KIND=dp))
         END DO

         IF (la_rec_max > 0) THEN
            DO i = lb_grid, ub_grid
               gcoord = pw%pw_grid%g(idir, i)
               rag = CMPLX(fap(idir), -gcoord, KIND=dp)
               cubeaxis(i, idir, 1, 0) = rag*cubeaxis(i, idir, 0, 0)
            END DO
            DO la = 2, la_rec_max
               fa = REAL(la - 1, dp)*twozetp
               DO i = lb_grid, ub_grid
                  gcoord = pw%pw_grid%g(idir, i)
                  rag = CMPLX(fap(idir), -gcoord, KIND=dp)
                  cubeaxis(i, idir, la, 0) = rag*cubeaxis(i, idir, la - 1, 0) + &
                                             fa*cubeaxis(i, idir, la - 2, 0)
               END DO
            END DO
            IF (lb_rec_max > 0) THEN
               fa = twozetp
               DO i = lb_grid, ub_grid
                  gcoord = pw%pw_grid%g(idir, i)
                  rbg = CMPLX(fbp(idir), -gcoord, KIND=dp)
                  cubeaxis(i, idir, 0, 1) = rbg*cubeaxis(i, idir, 0, 0)
                  cubeaxis(i, idir, 1, 1) = rbg*cubeaxis(i, idir, 1, 0) + &
                                            fa*cubeaxis(i, idir, 0, 0)
               END DO
               DO lb = 2, lb_rec_max
                  fb = REAL(lb - 1, dp)*twozetp
                  DO i = lb_grid, ub_grid
                     gcoord = pw%pw_grid%g(idir, i)
                     rbg = CMPLX(fbp(idir), -gcoord, KIND=dp)
                     cubeaxis(i, idir, 0, lb) = rbg*cubeaxis(i, idir, 0, lb - 1) + &
                                                fb*cubeaxis(i, idir, 0, lb - 2)
                     cubeaxis(i, idir, 1, lb) = rbg*cubeaxis(i, idir, 1, lb - 1) + &
                                                fb*cubeaxis(i, idir, 1, lb - 2) + &
                                                fa*cubeaxis(i, idir, 0, lb - 1)
                  END DO
               END DO
               DO la = 2, la_rec_max
                  fa = REAL(la, dp)*twozetp
                  DO i = lb_grid, ub_grid
                     gcoord = pw%pw_grid%g(idir, i)
                     rbg = CMPLX(fbp(idir), -gcoord, KIND=dp)
                     cubeaxis(i, idir, la, 1) = rbg*cubeaxis(i, idir, la, 0) + &
                                                fa*cubeaxis(i, idir, la - 1, 0)
                  END DO
                  DO lb = 2, lb_rec_max
                     fb = REAL(lb - 1, dp)*twozetp
                     DO i = lb_grid, ub_grid
                        gcoord = pw%pw_grid%g(idir, i)
                        rbg = CMPLX(fbp(idir), -gcoord, KIND=dp)
                        cubeaxis(i, idir, la, lb) = rbg*cubeaxis(i, idir, la, lb - 1) + &
                                                    fb*cubeaxis(i, idir, la, lb - 2) + &
                                                    fa*cubeaxis(i, idir, la - 1, lb - 1)
                     END DO
                  END DO
               END DO
            END IF
         ELSE
            IF (lb_rec_max > 0) THEN
               DO i = lb_grid, ub_grid
                  gcoord = pw%pw_grid%g(idir, i)
                  rbg = CMPLX(fbp(idir), -gcoord, KIND=dp)
                  cubeaxis(i, idir, 0, 1) = rbg*cubeaxis(i, idir, 0, 0)
               END DO
               DO lb = 2, lb_rec_max
                  fb = REAL(lb - 1, dp)*twozetp
                  DO i = lb_grid, ub_grid
                     gcoord = pw%pw_grid%g(idir, i)
                     rbg = CMPLX(fbp(idir), -gcoord, KIND=dp)
                     cubeaxis(i, idir, 0, lb) = rbg*cubeaxis(i, idir, 0, lb - 1) + &
                                                fb*cubeaxis(i, idir, 0, lb - 2)
                  END DO
               END DO
            END IF
         END IF
      END DO

      DO la = 0, la_rec_max
         DO lb = 0, lb_rec_max
            IF (la + lb == 0) CYCLE
            fa = (1.0_dp/twozetp)**(la + lb)
            DO idir = 1, 3
               DO i = lb_grid, ub_grid
                  cubeaxis(i, idir, la, lb) = fa*cubeaxis(i, idir, la, lb)
               END DO
            END DO
         END DO
      END DO

      ! Add the current primitive Gaussian function product to grid

      DO ico = ncoset(la_min - 1) + 1, ncoset(la_max)

         ax = indco(1, ico)
         ay = indco(2, ico)
         az = indco(3, ico)

         DO jco = ncoset(lb_min - 1) + 1, ncoset(lb_max)

            IF (integrate_potential) THEN
               g0_contribution = 0.0_dp
               potential_product = 0.0_dp
               force_product = 0.0_dp
               strain_product = 0.0_dp
               IF (PRESENT(center_force) .OR. integrate_strain) THEN
                  pij = prefactor*pab(na + ico, nb + jco)
               END IF
            ELSE
               pij = prefactor*pab(na + ico, nb + jco)
               IF (ABS(pij) < eps_rho_gspace) CYCLE
            END IF

            bx = indco(1, jco)
            by = indco(2, jco)
            bz = indco(3, jco)

            DO i = lb_grid, ub_grid
               IF (pw%pw_grid%gsq(i) > gsq_max) CYCLE
               IF (.NOT. my_compute_tau) THEN
                  plane_product = cubeaxis(i, 1, ax, bx)* &
                                  cubeaxis(i, 2, ay, by)* &
                                  cubeaxis(i, 3, az, bz)
               ELSE
                  ! Equivalent to GRID_FUNC_DADB: 0.5 * grad(phi_a).grad(phi_b).
                  tau_product = 2.0_dp*zeta*zetb* &
                                cubeaxis(i, 1, ax + 1, bx + 1)* &
                                cubeaxis(i, 2, ay, by)*cubeaxis(i, 3, az, bz)
                  IF (ax > 0) THEN
                     tau_product = tau_product - REAL(ax, dp)*zetb* &
                                   cubeaxis(i, 1, ax - 1, bx + 1)* &
                                   cubeaxis(i, 2, ay, by)*cubeaxis(i, 3, az, bz)
                  END IF
                  IF (bx > 0) THEN
                     tau_product = tau_product - zeta*REAL(bx, dp)* &
                                   cubeaxis(i, 1, ax + 1, bx - 1)* &
                                   cubeaxis(i, 2, ay, by)*cubeaxis(i, 3, az, bz)
                  END IF
                  IF (ax > 0 .AND. bx > 0) THEN
                     tau_product = tau_product + 0.5_dp*REAL(ax*bx, dp)* &
                                   cubeaxis(i, 1, ax - 1, bx - 1)* &
                                   cubeaxis(i, 2, ay, by)*cubeaxis(i, 3, az, bz)
                  END IF

                  tau_product = tau_product + 2.0_dp*zeta*zetb* &
                                cubeaxis(i, 1, ax, bx)* &
                                cubeaxis(i, 2, ay + 1, by + 1)*cubeaxis(i, 3, az, bz)
                  IF (ay > 0) THEN
                     tau_product = tau_product - REAL(ay, dp)*zetb* &
                                   cubeaxis(i, 1, ax, bx)* &
                                   cubeaxis(i, 2, ay - 1, by + 1)*cubeaxis(i, 3, az, bz)
                  END IF
                  IF (by > 0) THEN
                     tau_product = tau_product - zeta*REAL(by, dp)* &
                                   cubeaxis(i, 1, ax, bx)* &
                                   cubeaxis(i, 2, ay + 1, by - 1)*cubeaxis(i, 3, az, bz)
                  END IF
                  IF (ay > 0 .AND. by > 0) THEN
                     tau_product = tau_product + 0.5_dp*REAL(ay*by, dp)* &
                                   cubeaxis(i, 1, ax, bx)* &
                                   cubeaxis(i, 2, ay - 1, by - 1)*cubeaxis(i, 3, az, bz)
                  END IF

                  tau_product = tau_product + 2.0_dp*zeta*zetb* &
                                cubeaxis(i, 1, ax, bx)*cubeaxis(i, 2, ay, by)* &
                                cubeaxis(i, 3, az + 1, bz + 1)
                  IF (az > 0) THEN
                     tau_product = tau_product - REAL(az, dp)*zetb* &
                                   cubeaxis(i, 1, ax, bx)*cubeaxis(i, 2, ay, by)* &
                                   cubeaxis(i, 3, az - 1, bz + 1)
                  END IF
                  IF (bz > 0) THEN
                     tau_product = tau_product - zeta*REAL(bz, dp)* &
                                   cubeaxis(i, 1, ax, bx)*cubeaxis(i, 2, ay, by)* &
                                   cubeaxis(i, 3, az + 1, bz - 1)
                  END IF
                  IF (az > 0 .AND. bz > 0) THEN
                     tau_product = tau_product + 0.5_dp*REAL(az*bz, dp)* &
                                   cubeaxis(i, 1, ax, bx)*cubeaxis(i, 2, ay, by)* &
                                   cubeaxis(i, 3, az - 1, bz - 1)
                  END IF
                  plane_product = tau_product
               END IF
               IF (integrate_strain) THEN
                  IF (.NOT. my_compute_tau) THEN
                     shifted_product = cubeaxis(i, 1, ax + 1, bx)* &
                                       cubeaxis(i, 2, ay, by)*cubeaxis(i, 3, az, bz)
                     plane_gradient(1) = CMPLX(0.0_dp, -1.0_dp, KIND=dp)* &
                                         (shifted_product + ra(1)*plane_product)
                     shifted_product = cubeaxis(i, 1, ax, bx)* &
                                       cubeaxis(i, 2, ay + 1, by)*cubeaxis(i, 3, az, bz)
                     plane_gradient(2) = CMPLX(0.0_dp, -1.0_dp, KIND=dp)* &
                                         (shifted_product + ra(2)*plane_product)
                     shifted_product = cubeaxis(i, 1, ax, bx)* &
                                       cubeaxis(i, 2, ay, by)*cubeaxis(i, 3, az + 1, bz)
                     plane_gradient(3) = CMPLX(0.0_dp, -1.0_dp, KIND=dp)* &
                                         (shifted_product + ra(3)*plane_product)
                  ELSE
                     DO idir = 1, 3
                        shifted_product = tau_plane_product_shifted(idir)
                        plane_gradient(idir) = CMPLX(0.0_dp, -1.0_dp, KIND=dp)* &
                                               (shifted_product + ra(idir)*plane_product)
                     END DO
                  END IF
               END IF
               IF (integrate_potential) THEN
                  potential_product = potential_product + prefactor* &
                                      REAL(CONJG(potential_gspace%array(i))*plane_product, KIND=dp)
                  IF (PRESENT(center_force) .OR. integrate_strain) THEN
                     gvec = pw%pw_grid%g(:, i)
                     force_product = force_product + pij* &
                                     REAL(CONJG(potential_gspace%array(i))* &
                                          CMPLX(0.0_dp, -gvec, KIND=dp)*plane_product, KIND=dp)
                  END IF
                  IF (integrate_strain) THEN
                     DO idir = 1, 3
                        DO jdir = 1, 3
                           strain_product(idir, jdir) = strain_product(idir, jdir) - &
                                                        pij*gvec(idir)*REAL(CONJG(potential_gspace%array(i))* &
                                                                            plane_gradient(jdir), KIND=dp)
                        END DO
                     END DO
                  END IF
                  IF (pw%pw_grid%have_g0 .AND. i == 1) THEN
                     g0_contribution = prefactor* &
                                       REAL(CONJG(potential_gspace%array(i))*plane_product, KIND=dp)
                  END IF
               ELSE
                  pw%array(i) = pw%array(i) + pij*plane_product
               END IF
            END DO

            IF (integrate_potential) THEN
               IF (pw%pw_grid%grid_span == HALFSPACE) THEN
                  potential_product = 2.0_dp*potential_product - g0_contribution
                  force_product = 2.0_dp*force_product
                  strain_product = 2.0_dp*strain_product
               END IF
               hab(na + ico, nb + jco) = hab(na + ico, nb + jco) + potential_product
               IF (PRESENT(center_force)) center_force = center_force + force_product
               IF (integrate_strain) THEN
                  DO idir = 1, 3
                     DO jdir = 1, 3
                        strain_product(idir, jdir) = strain_product(idir, jdir) + &
                                                     force_product(idir)*ra(jdir)
                     END DO
                     strain_product(idir, idir) = strain_product(idir, idir) - &
                                                  pab(na + ico, nb + jco)*potential_product
                  END DO
                  strain_virial = strain_virial + strain_product
               END IF
            END IF

         END DO

      END DO

      DEALLOCATE (cubeaxis)

      CALL timestop(handle)

   CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param shift_axis ...
!> \return ...
! **************************************************************************************************
      FUNCTION tau_plane_product_shifted(shift_axis) RESULT(value)
      INTEGER, INTENT(IN)                                :: shift_axis
      COMPLEX(dp)                                        :: value

      INTEGER                                            :: sx, sy, sz

         sx = 0
         sy = 0
         sz = 0
         IF (shift_axis == 1) sx = 1
         IF (shift_axis == 2) sy = 1
         IF (shift_axis == 3) sz = 1

         value = 2.0_dp*zeta*zetb* &
                 cubeaxis(i, 1, ax + 1 + sx, bx + 1)* &
                 cubeaxis(i, 2, ay + sy, by)*cubeaxis(i, 3, az + sz, bz)
         IF (ax > 0) THEN
            value = value - REAL(ax, dp)*zetb* &
                    cubeaxis(i, 1, ax - 1 + sx, bx + 1)* &
                    cubeaxis(i, 2, ay + sy, by)*cubeaxis(i, 3, az + sz, bz)
         END IF
         IF (bx > 0) THEN
            value = value - zeta*REAL(bx, dp)* &
                    cubeaxis(i, 1, ax + 1 + sx, bx - 1)* &
                    cubeaxis(i, 2, ay + sy, by)*cubeaxis(i, 3, az + sz, bz)
         END IF
         IF (ax > 0 .AND. bx > 0) THEN
            value = value + 0.5_dp*REAL(ax*bx, dp)* &
                    cubeaxis(i, 1, ax - 1 + sx, bx - 1)* &
                    cubeaxis(i, 2, ay + sy, by)*cubeaxis(i, 3, az + sz, bz)
         END IF

         value = value + 2.0_dp*zeta*zetb* &
                 cubeaxis(i, 1, ax + sx, bx)* &
                 cubeaxis(i, 2, ay + 1 + sy, by + 1)*cubeaxis(i, 3, az + sz, bz)
         IF (ay > 0) THEN
            value = value - REAL(ay, dp)*zetb* &
                    cubeaxis(i, 1, ax + sx, bx)* &
                    cubeaxis(i, 2, ay - 1 + sy, by + 1)*cubeaxis(i, 3, az + sz, bz)
         END IF
         IF (by > 0) THEN
            value = value - zeta*REAL(by, dp)* &
                    cubeaxis(i, 1, ax + sx, bx)* &
                    cubeaxis(i, 2, ay + 1 + sy, by - 1)*cubeaxis(i, 3, az + sz, bz)
         END IF
         IF (ay > 0 .AND. by > 0) THEN
            value = value + 0.5_dp*REAL(ay*by, dp)* &
                    cubeaxis(i, 1, ax + sx, bx)* &
                    cubeaxis(i, 2, ay - 1 + sy, by - 1)*cubeaxis(i, 3, az + sz, bz)
         END IF

         value = value + 2.0_dp*zeta*zetb* &
                 cubeaxis(i, 1, ax + sx, bx)*cubeaxis(i, 2, ay + sy, by)* &
                 cubeaxis(i, 3, az + 1 + sz, bz + 1)
         IF (az > 0) THEN
            value = value - REAL(az, dp)*zetb* &
                    cubeaxis(i, 1, ax + sx, bx)*cubeaxis(i, 2, ay + sy, by)* &
                    cubeaxis(i, 3, az - 1 + sz, bz + 1)
         END IF
         IF (bz > 0) THEN
            value = value - zeta*REAL(bz, dp)* &
                    cubeaxis(i, 1, ax + sx, bx)*cubeaxis(i, 2, ay + sy, by)* &
                    cubeaxis(i, 3, az + 1 + sz, bz - 1)
         END IF
         IF (az > 0 .AND. bz > 0) THEN
            value = value + 0.5_dp*REAL(az*bz, dp)* &
                    cubeaxis(i, 1, ax + sx, bx)*cubeaxis(i, 2, ay + sy, by)* &
                    cubeaxis(i, 3, az - 1 + sz, bz - 1)
         END IF

      END FUNCTION tau_plane_product_shifted

   END SUBROUTINE collocate_pgf_product_gspace_buffered

END MODULE gapw_gspace_reconstruction
