!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief qs_environement methods that use many other modules
!> \par History
!>      09.2002 created [fawzi]
!>      - local atom distribution (25.06.2003,MK)
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE qs_environment_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_dist2d_to_dist
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_distribution_type
   USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                              distribution_2d_type
   USE distribution_methods,            ONLY: distribute_molecules_2d
   USE ewald_environment_types,         ONLY: ewald_environment_type
   USE ewald_pw_methods,                ONLY: ewald_pw_grid_update
   USE ewald_pw_types,                  ONLY: ewald_pw_type
   USE input_constants,                 ONLY: do_ppl_grid
   USE kinds,                           ONLY: dp
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_types,                  ONLY: molecule_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_methods,                  ONLY: pw_env_create,&
                                              pw_env_rebuild
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_release,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_release,&
                                              pw_type
   USE qs_charges_types,                ONLY: qs_charges_create,&
                                              qs_charges_release,&
                                              qs_charges_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              has_nlcc,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_matrix_pools,                 ONLY: mpools_rebuild_fm_pools
   USE qs_outer_scf,                    ONLY: outer_loop_variables_count
   USE qs_rho0_types,                   ONLY: rho0_mpole_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_environment_methods'

   PUBLIC :: qs_env_rebuild_pw_env, &
             qs_env_setup
!***
CONTAINS

! **************************************************************************************************
!> \brief initializes various components of the qs_env, that need only
!>      atomic_kind_set, cell, dft_control, scf_control, c(i)%nmo,
!>      c(i)%nao, and particle_set to be initialized.
!>      The previous components of qs_env must be valid.
!>      Initializes pools, charges and pw_env.
!> \param qs_env the qs_env to set up
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE qs_env_setup(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_setup', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, n_ao_aux_fit, nhistory, &
                                                            nvariables
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: outer_scf_history
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_charges_type), POINTER                     :: qs_charges
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, atomic_kind_set, dft_control, scf_control, qs_charges, para_env, &
               distribution_2d, molecule_kind_set, molecule_set, particle_set, cell, &
               ks_env, blacs_env)

      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      dft_control=dft_control, &
                      molecule_kind_set=molecule_kind_set, &
                      molecule_set=molecule_set, &
                      particle_set=particle_set, &
                      scf_control=scf_control, &
                      para_env=para_env, &
                      blacs_env=blacs_env, &
                      cell=cell, &
                      ks_env=ks_env)

      CPASSERT(ASSOCIATED(qs_kind_set))
      CPASSERT(ASSOCIATED(atomic_kind_set))
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(scf_control))
      ! allocate qs_charges
      CALL qs_charges_create(qs_charges, nspins=dft_control%nspins)
      CALL set_qs_env(qs_env, qs_charges=qs_charges)
      CALL qs_charges_release(qs_charges)

      ! outer scf setup
      IF (scf_control%outer_scf%have_scf) THEN
         nvariables = outer_loop_variables_count(scf_control)
         nhistory = scf_control%outer_scf%extrapolation_order
         ALLOCATE (outer_scf_history(nvariables, nhistory))
         CALL set_qs_env(qs_env, outer_scf_history=outer_scf_history)
         CALL set_qs_env(qs_env, outer_scf_ihistory=0)
      ENDIF

      ! set up pw_env
      CALL qs_env_rebuild_pw_env(qs_env)

      ! rebuilds fm_pools

      ! XXXX should get rid of the mpools
      IF (ASSOCIATED(qs_env%mos)) THEN
         CALL mpools_rebuild_fm_pools(qs_env%mpools, mos=qs_env%mos, &
                                      blacs_env=blacs_env, para_env=para_env)
      ENDIF

      ! If we use auxiliary density matrix methods rebuild fm_pools
      IF (dft_control%do_admm) THEN
         CALL get_qs_kind_set(qs_kind_set, nsgf=n_ao_aux_fit, &
                              basis_type="AUX_FIT")
         CALL mpools_rebuild_fm_pools(qs_env%mpools_aux_fit, mos=qs_env%mos_aux_fit, &
                                      blacs_env=blacs_env, para_env=para_env)
      END IF

      ! create 2d distribution

      CALL distribute_molecules_2d(cell=cell, &
                                   atomic_kind_set=atomic_kind_set, &
                                   qs_kind_set=qs_kind_set, &
                                   particle_set=particle_set, &
                                   molecule_kind_set=molecule_kind_set, &
                                   molecule_set=molecule_set, &
                                   distribution_2d=distribution_2d, &
                                   blacs_env=blacs_env, &
                                   force_env_section=qs_env%input)

      ! and use it to create the dbcsr_dist, which should be the sole user of distribution_2d by now.
      ALLOCATE (dbcsr_dist)
      CALL cp_dbcsr_dist2d_to_dist(distribution_2d, dbcsr_dist)
      CALL set_ks_env(ks_env, dbcsr_dist=dbcsr_dist)

      ! also keep distribution_2d in qs_env
      CALL set_ks_env(ks_env, distribution_2d=distribution_2d)
      CALL distribution_2d_release(distribution_2d)

      CALL timestop(handle)

   END SUBROUTINE qs_env_setup

! **************************************************************************************************
!> \brief rebuilds the pw_env in the given qs_env, allocating it if necessary
!> \param qs_env the qs_env whose pw_env has to be rebuilt
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE qs_env_rebuild_pw_env(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_rebuild_pw_env', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: nlcc
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(pw_env_type), POINTER                         :: new_pw_env
      TYPE(pw_p_type), POINTER                           :: external_vxc, rho0_gs, rho0_rs, &
                                                            rho_core, rho_nlcc, rho_nlcc_g, vee, &
                                                            vppl
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type), POINTER                             :: v_hartree_rspace
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(rho0_mpole_type), POINTER                     :: rho0_mpole

      CALL timeset(routineN, handle)
      ! rebuild pw_env
      NULLIFY (dft_control, cell, ks_env, v_hartree_rspace, auxbas_pw_pool)
      NULLIFY (rho0_mpole, rho0_gs, rho0_rs)
      NULLIFY (ewald_env, ewald_pw, new_pw_env, external_vxc, rho_core, rho_nlcc, rho_nlcc_g, vee, vppl)

      CALL get_qs_env(qs_env, ks_env=ks_env, pw_env=new_pw_env)
      IF (.NOT. ASSOCIATED(new_pw_env)) THEN
         CALL pw_env_create(new_pw_env)
         CALL set_ks_env(ks_env, pw_env=new_pw_env)
         CALL pw_env_release(new_pw_env)
      ENDIF

      CALL get_qs_env(qs_env, pw_env=new_pw_env, dft_control=dft_control, &
                      cell=cell)

      IF (ANY(new_pw_env%cell_hmat /= cell%hmat)) THEN
         ! only rebuild if nessecary
         new_pw_env%cell_hmat = cell%hmat
         CALL pw_env_rebuild(new_pw_env, qs_env=qs_env)

         ! reallocate rho_core
         CALL get_qs_env(qs_env, pw_env=new_pw_env, rho_core=rho_core)
         CPASSERT(ASSOCIATED(new_pw_env))
         IF (dft_control%qs_control%gapw) THEN
            IF (ASSOCIATED(rho_core)) THEN
               CALL pw_release(rho_core%pw)
               DEALLOCATE (rho_core)
            END IF
            IF (dft_control%qs_control%gapw_control%nopaw_as_gpw) THEN
               ALLOCATE (rho_core)
               CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool)
               CALL pw_pool_create_pw(auxbas_pw_pool, rho_core%pw, &
                                      use_data=COMPLEXDATA1D)
               rho_core%pw%in_space = RECIPROCALSPACE
               CALL set_ks_env(ks_env, rho_core=rho_core)
            END IF
            CALL get_qs_env(qs_env=qs_env, rho0_mpole=rho0_mpole)
            IF (ASSOCIATED(rho0_mpole)) THEN
               rho0_rs => rho0_mpole%rho0_s_rs
               rho0_gs => rho0_mpole%rho0_s_gs
               IF (ASSOCIATED(rho0_rs)) THEN
                  CALL pw_release(rho0_rs%pw)
                  DEALLOCATE (rho0_rs)
               END IF
               ALLOCATE (rho0_rs)
               CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool)
               CALL pw_pool_create_pw(auxbas_pw_pool, rho0_rs%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               rho0_mpole%rho0_s_rs => rho0_rs

               IF (ASSOCIATED(rho0_gs)) THEN
                  CALL pw_release(rho0_gs%pw)
                  DEALLOCATE (rho0_gs)
               END IF
               ALLOCATE (rho0_gs)
               CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool)
               CALL pw_pool_create_pw(auxbas_pw_pool, rho0_gs%pw, &
                                      use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
               rho0_mpole%rho0_s_gs => rho0_gs
            END IF
         ELSE IF (dft_control%qs_control%semi_empirical) THEN
            IF (dft_control%qs_control%se_control%do_ewald .OR. &
                dft_control%qs_control%se_control%do_ewald_gks) THEN
               ! rebuild Ewald environment
               CALL get_qs_env(qs_env=qs_env, ewald_env=ewald_env, ewald_pw=ewald_pw)
               CALL ewald_pw_grid_update(ewald_pw, ewald_env, cell%hmat)
            END IF
         ELSE IF (dft_control%qs_control%dftb) THEN
            IF (dft_control%qs_control%dftb_control%do_ewald) THEN
               ! rebuild Ewald environment
               CALL get_qs_env(qs_env=qs_env, ewald_env=ewald_env, ewald_pw=ewald_pw)
               CALL ewald_pw_grid_update(ewald_pw, ewald_env, cell%hmat)
            END IF
         ELSE IF (dft_control%qs_control%scptb) THEN
            IF (dft_control%qs_control%scptb_control%do_ewald) THEN
               ! rebuild Ewald environment
               CALL get_qs_env(qs_env=qs_env, ewald_env=ewald_env, ewald_pw=ewald_pw)
               CALL ewald_pw_grid_update(ewald_pw, ewald_env, cell%hmat)
            END IF
         ELSE
            IF (ASSOCIATED(rho_core)) THEN
               CALL pw_release(rho_core%pw)
               DEALLOCATE (rho_core)
            ENDIF
            ALLOCATE (rho_core)
            CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool)
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_core%pw, &
                                   use_data=COMPLEXDATA1D)
            rho_core%pw%in_space = RECIPROCALSPACE
            CALL set_ks_env(ks_env, rho_core=rho_core)
         END IF

         ! reallocate vppl (realspace grid of local pseudopotential
         IF (dft_control%qs_control%do_ppl_method == do_ppl_grid) THEN
            NULLIFY (vppl)
            CALL get_qs_env(qs_env, pw_env=new_pw_env, vppl=vppl)
            IF (ASSOCIATED(vppl)) THEN
               CALL pw_release(vppl%pw)
               DEALLOCATE (vppl)
            ENDIF
            ALLOCATE (vppl)
            CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool)
            CALL pw_pool_create_pw(auxbas_pw_pool, vppl%pw, use_data=REALDATA3D)
            vppl%pw%in_space = REALSPACE
            CALL set_ks_env(ks_env, vppl=vppl)
         END IF

         ! reallocate rho_nlcc
         CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
         nlcc = has_nlcc(qs_kind_set)
         IF (nlcc) THEN
            ! the realspace version
            NULLIFY (rho_nlcc)
            CALL get_qs_env(qs_env, pw_env=new_pw_env, rho_nlcc=rho_nlcc)
            IF (ASSOCIATED(rho_nlcc)) THEN
               CALL pw_release(rho_nlcc%pw)
               DEALLOCATE (rho_nlcc)
            ENDIF
            ALLOCATE (rho_nlcc)
            CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool)
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_nlcc%pw, use_data=REALDATA3D)
            rho_nlcc%pw%in_space = REALSPACE
            CALL set_ks_env(ks_env, rho_nlcc=rho_nlcc)
            ! the g-space version
            NULLIFY (rho_nlcc_g)
            CALL get_qs_env(qs_env, pw_env=new_pw_env, rho_nlcc_g=rho_nlcc_g)
            IF (ASSOCIATED(rho_nlcc_g)) THEN
               CALL pw_release(rho_nlcc_g%pw)
               DEALLOCATE (rho_nlcc_g)
            ENDIF
            ALLOCATE (rho_nlcc_g)
            CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool)
            CALL pw_pool_create_pw(auxbas_pw_pool, rho_nlcc_g%pw, use_data=COMPLEXDATA1D)
            rho_nlcc_g%pw%in_space = RECIPROCALSPACE
            CALL set_ks_env(ks_env, rho_nlcc_g=rho_nlcc_g)
         END IF

         ! reallocate vee: external electrostatic potential
         IF (dft_control%apply_external_potential) THEN
            NULLIFY (vee)
            CALL get_qs_env(qs_env, pw_env=new_pw_env, vee=vee)
            IF (ASSOCIATED(vee)) THEN
               CALL pw_release(vee%pw)
               DEALLOCATE (vee)
            ENDIF
            ALLOCATE (vee)
            CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool)
            CALL pw_pool_create_pw(auxbas_pw_pool, vee%pw, use_data=REALDATA3D)
            vee%pw%in_space = REALSPACE
            CALL set_ks_env(ks_env, vee=vee)
            dft_control%eval_external_potential = .TRUE.
         END IF

         ! ZMP Reallocate external_vxc: external vxc potential
         IF (dft_control%apply_external_vxc) THEN
            NULLIFY (external_vxc)
            CALL get_qs_env(qs_env, pw_env=new_pw_env, external_vxc=external_vxc)
            IF (ASSOCIATED(external_vxc)) THEN
               CALL pw_release(external_vxc%pw)
               DEALLOCATE (external_vxc)
            ENDIF
            ALLOCATE (external_vxc)
            CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool)
            CALL pw_pool_create_pw(auxbas_pw_pool, external_vxc%pw, use_data=REALDATA3D)
            external_vxc%pw%in_space = REALSPACE
            CALL set_qs_env(qs_env, external_vxc=external_vxc)
            dft_control%read_external_vxc = .TRUE.
         END IF

         CALL get_ks_env(ks_env, v_hartree_rspace=v_hartree_rspace)
         IF (ASSOCIATED(v_hartree_rspace)) &
            CALL pw_release(v_hartree_rspace)
         CALL get_qs_env(qs_env, pw_env=new_pw_env)
         CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL set_ks_env(ks_env, v_hartree_rspace=v_hartree_rspace)
      ENDIF

      !update the time in the poisson environement, to update time dependant constraints
      new_pw_env%poisson_env%parameters%dbc_params%time = qs_env%sim_time

      CALL timestop(handle)

   END SUBROUTINE qs_env_rebuild_pw_env

END MODULE qs_environment_methods
