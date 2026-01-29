!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!> \author JGH (01.2026)
! **************************************************************************************************
MODULE accint_weights_forces
   USE ao_util,                         ONLY: exp_radius_very_extended
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE grid_api,                        ONLY: integrate_pgf_product
   USE input_constants,                 ONLY: sic_none,&
                                              xc_none
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE memory_utilities,                ONLY: reallocate
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grids,                        ONLY: pw_grid_compare
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_multiply_with,&
                                              pw_scale,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_fxc,                          ONLY: qs_fxc_analytic
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                              transfer_pw2rs
   USE virial_types,                    ONLY: virial_type
   USE xc,                              ONLY: xc_exc_pw_create,&
                                              xc_vxc_pw_create
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'accint_weights_forces'

   PUBLIC :: accint_weight_force

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho ...
!> \param rho1 ...
!> \param order ...
!> \param xc_section ...
!> \param triplet ...
! **************************************************************************************************
   SUBROUTINE accint_weight_force(qs_env, rho, rho1, order, xc_section, triplet)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho1
      INTEGER, INTENT(IN)                                :: order
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN), OPTIONAL                      :: triplet

      CHARACTER(len=*), PARAMETER :: routineN = 'accint_weight_force'

      INTEGER                                            :: atom_a, handle, iatom, ikind, natom, &
                                                            natom_of_kind, nkind
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: lr_triplet, use_virial
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: calpha, cvalue
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aforce
      REAL(KIND=dp), DIMENSION(3, 3)                     :: avirial
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: e_rspace
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)

      IF (dft_control%qs_control%gapw_control%accurate_xcint) THEN

         CALL get_qs_env(qs_env=qs_env, force=force, virial=virial)
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

         CALL get_qs_env(qs_env, natom=natom, nkind=nkind)
         ALLOCATE (aforce(3, natom))
         ALLOCATE (calpha(nkind), cvalue(nkind))
         cvalue = 1.0_dp
         calpha(1:nkind) = dft_control%qs_control%gapw_control%aw(1:nkind)

         CALL get_qs_env(qs_env, ks_env=ks_env, pw_env=pw_env)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL auxbas_pw_pool%create_pw(e_rspace)

         lr_triplet = .FALSE.
         IF (PRESENT(triplet)) lr_triplet = triplet

         CALL xc_density(ks_env, rho, rho1, order, xc_section, lr_triplet, e_rspace)
         CALL pw_scale(e_rspace, e_rspace%pw_grid%dvol)

         CALL gauss_grid_force(e_rspace, qs_env, calpha, cvalue, aforce, avirial)

         CALL auxbas_pw_pool%give_back_pw(e_rspace)

         CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom_of_kind, atom_list=atom_list)
            DO iatom = 1, natom_of_kind
               atom_a = atom_list(iatom)
               force(ikind)%rho_elec(1:3, iatom) = &
                  force(ikind)%rho_elec(1:3, iatom) + aforce(1:3, atom_a)
            END DO
         END DO
         IF (use_virial) THEN
            virial%pv_exc = virial%pv_exc + avirial
            virial%pv_virial = virial%pv_virial + avirial
         END IF

         DEALLOCATE (aforce)

      END IF

      CALL timestop(handle)

   END SUBROUTINE accint_weight_force

! **************************************************************************************************
!> \brief computes the forces/virial due to atomic centered Gaussian functions
!> \param e_rspace Energy density
!> \param qs_env ...
!> \param calpha ...
!> \param cvalue ...
!> \param aforce ...
!> \param avirial ...
! **************************************************************************************************
   SUBROUTINE gauss_grid_force(e_rspace, qs_env, calpha, cvalue, aforce, avirial)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: e_rspace
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: calpha, cvalue
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: aforce
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: avirial

      CHARACTER(len=*), PARAMETER                        :: routineN = 'gauss_grid_force'

      INTEGER                                            :: atom_a, handle, iatom, ikind, j, &
                                                            natom_of_kind, npme
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: alpha, eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, ra
      REAL(KIND=dp), DIMENSION(3, 3)                     :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hab, pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(realspace_grid_type), POINTER                 :: rs_v

      CALL timeset(routineN, handle)

      ALLOCATE (cores(1))
      ALLOCATE (hab(1, 1))
      ALLOCATE (pab(1, 1))

      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_v)

      CALL transfer_pw2rs(rs_v, e_rspace)

      CALL get_qs_env(qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set)

      use_virial = .TRUE.
      avirial = 0.0_dp
      aforce = 0.0_dp

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      DO ikind = 1, SIZE(atomic_kind_set)

         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom_of_kind, atom_list=atom_list)

         alpha = calpha(ikind)
         pab(1, 1) = -cvalue(ikind)
         IF (alpha == 0.0_dp .OR. pab(1, 1) == 0.0_dp) CYCLE

         CALL reallocate(cores, 1, natom_of_kind)
         npme = 0
         cores = 0

         DO iatom = 1, natom_of_kind
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r, cell)
            IF (rs_v%desc%parallel .AND. .NOT. rs_v%desc%distributed) THEN
               ! replicated realspace grid, split the atoms up between procs
               IF (MODULO(iatom, rs_v%desc%group_size) == rs_v%desc%my_pos) THEN
                  npme = npme + 1
                  cores(npme) = iatom
               END IF
            ELSE
               npme = npme + 1
               cores(npme) = iatom
            END IF
         END DO

         DO j = 1, npme

            iatom = cores(j)
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r, cell)
            hab(1, 1) = 0.0_dp
            force_a(:) = 0.0_dp
            force_b(:) = 0.0_dp
            my_virial_a = 0.0_dp
            my_virial_b = 0.0_dp

            radius = exp_radius_very_extended(la_min=0, la_max=0, lb_min=0, lb_max=0, &
                                              ra=ra, rb=ra, rp=ra, &
                                              zetp=alpha, eps=eps_rho_rspace, &
                                              pab=pab, o1=0, o2=0, &
                                              prefactor=1.0_dp, cutoff=1.0_dp)

            CALL integrate_pgf_product(0, alpha, 0, &
                                       0, 0.0_dp, 0, ra, [0.0_dp, 0.0_dp, 0.0_dp], &
                                       rs_v, hab, pab=pab, o1=0, o2=0, &
                                       radius=radius, &
                                       calculate_forces=.TRUE., force_a=force_a, &
                                       force_b=force_b, use_virial=use_virial, my_virial_a=my_virial_a, &
                                       my_virial_b=my_virial_b, use_subpatch=.TRUE., subpatch_pattern=0)

            aforce(1:3, atom_a) = aforce(1:3, atom_a) + force_a(1:3)
            avirial = avirial + my_virial_a

         END DO

      END DO

      DEALLOCATE (hab, pab, cores)

      CALL timestop(handle)

   END SUBROUTINE gauss_grid_force

! **************************************************************************************************
!> \brief calculates the XC density:
!>        order=0: exc will contain the xc energy density E_xc(r)
!>        order=1: exc will contain V_xc(r) * rho1(r)
!>        order=2: exc will contain F_xc(r) * rho1(r) * rho1(r)
!> \param ks_env to get all the needed things
!> \param rho_struct density
!> \param rho1_struct response density
!> \param order requested derivative order
!> \param xc_section ...
!> \param triplet ...
!> \param exc ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE xc_density(ks_env, rho_struct, rho1_struct, order, xc_section, triplet, exc)

      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_struct, rho1_struct
      INTEGER, INTENT(IN)                                :: order
      TYPE(section_vals_type), POINTER                   :: xc_section
      LOGICAL, INTENT(IN)                                :: triplet
      TYPE(pw_r3d_rs_type)                               :: exc

      CHARACTER(len=*), PARAMETER                        :: routineN = 'xc_density'

      INTEGER                                            :: handle, ispin, myfun, nspins
      LOGICAL                                            :: uf_grid
      REAL(KIND=dp)                                      :: excint, factor
      REAL(KIND=dp), DIMENSION(3, 3)                     :: vdum
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_nlcc_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool, xc_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho1_r, rho_r, tau1_r, tau_r, vxc_rho, &
                                                            vxc_tau
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_nlcc, weights

      CALL timeset(routineN, handle)

      ! we always get true exc (not integration weighted)
      NULLIFY (weights)

      CALL get_ks_env(ks_env, &
                      dft_control=dft_control, &
                      pw_env=pw_env, &
                      cell=cell, &
                      rho_nlcc=rho_nlcc, &
                      rho_nlcc_g=rho_nlcc_g)

      CALL qs_rho_get(rho_struct, rho_r=rho_r, rho_g=rho_g, tau_r=tau_r)

      nspins = dft_control%nspins

      CALL section_vals_val_get(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", i_val=myfun)

      CALL pw_env_get(pw_env, xc_pw_pool=xc_pw_pool, auxbas_pw_pool=auxbas_pw_pool)
      uf_grid = .NOT. pw_grid_compare(auxbas_pw_pool%pw_grid, xc_pw_pool%pw_grid)
      IF (uf_grid) THEN
         CALL cp_warn(__LOCATION__, "Fine grid option not possible with energy density")
         CPABORT("Fine Grid in xc_density")
      END IF

      CALL pw_zero(exc)

      IF (myfun /= xc_none) THEN

         CPASSERT(ASSOCIATED(rho_struct))
         CPASSERT(dft_control%sic_method_id == sic_none)

         ! add the nlcc densities
         IF (ASSOCIATED(rho_nlcc) .AND. order <= 1) THEN
            factor = 1.0_dp
            DO ispin = 1, nspins
               CALL pw_axpy(rho_nlcc, rho_r(ispin), factor)
               CALL pw_axpy(rho_nlcc_g, rho_g(ispin), factor)
            END DO
         END IF

         NULLIFY (vxc_rho, vxc_tau)
         SELECT CASE (order)
         CASE (0)
            ! we could reduce to energy only here
            CALL xc_exc_pw_create(rho_r, rho_g, tau_r, xc_section, weights, xc_pw_pool, exc)
         CASE (1)
            CALL qs_rho_get(rho1_struct, rho_r=rho1_r, tau_r=tau1_r)
            CALL xc_vxc_pw_create(vxc_rho=vxc_rho, vxc_tau=vxc_tau, rho_r=rho_r, &
                                  rho_g=rho_g, tau=tau_r, exc=excint, &
                                  xc_section=xc_section, &
                                  weights=weights, pw_pool=xc_pw_pool, &
                                  compute_virial=.FALSE., &
                                  virial_xc=vdum)
         CASE (2)
            CALL qs_rho_get(rho1_struct, rho_r=rho1_r, tau_r=tau1_r)
            CALL qs_fxc_analytic(rho_struct, rho1_r, tau1_r, xc_section, weights, xc_pw_pool, &
                                 triplet, vxc_rho, vxc_tau)
         CASE DEFAULT
            CPABORT("Derivative order not available in xc_density")
         END SELECT

         ! remove the nlcc densities (keep stuff in original state)
         IF (ASSOCIATED(rho_nlcc) .AND. order <= 1) THEN
            factor = -1.0_dp
            DO ispin = 1, dft_control%nspins
               CALL pw_axpy(rho_nlcc, rho_r(ispin), factor)
               CALL pw_axpy(rho_nlcc_g, rho_g(ispin), factor)
            END DO
         END IF
         !
         SELECT CASE (order)
         CASE (0)
            !
         CASE (1, 2)
            CALL pw_zero(exc)
            IF (ASSOCIATED(vxc_rho)) THEN
               DO ispin = 1, nspins
                  CALL pw_multiply_with(vxc_rho(ispin), rho1_r(ispin))
                  CALL pw_axpy(vxc_rho(ispin), exc, 1.0_dp)
                  CALL vxc_rho(ispin)%release()
               END DO
               DEALLOCATE (vxc_rho)
            END IF
            IF (ASSOCIATED(vxc_tau)) THEN
               DO ispin = 1, nspins
                  CALL pw_multiply_with(vxc_tau(ispin), tau1_r(ispin))
                  CALL pw_axpy(vxc_tau(ispin), exc, 1.0_dp)
                  CALL vxc_tau(ispin)%release()
               END DO
               DEALLOCATE (vxc_tau)
            END IF
         CASE DEFAULT
            CPABORT("Derivative order not available in xc_density")
         END SELECT

         IF (order == 2) THEN
            CALL pw_scale(exc, 0.5_dp)
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE xc_density

END MODULE accint_weights_forces
