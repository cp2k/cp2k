!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Methods to performs a path integral run
!> \author fawzi
!> \par History
!>      02.2005 created [fawzi]
!>           11.2006 modified so it might actually work [hforbert]
!>           10.2015 added RPMD propagator
!>           10.2015 added exact harmonic integrator [Felix Uhl]
!> \note   quick & dirty rewrite of my python program
! **************************************************************************************************
MODULE pint_methods

   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE bibliography,                    ONLY: Brieuc2016,&
                                              Ceriotti2010,&
                                              Ceriotti2012,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE constraint,                      ONLY: rattle_control,&
                                              shake_control,&
                                              shake_update_targets
   USE constraint_util,                 ONLY: getold
   USE cp_external_control,             ONLY: external_control
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              cp_rm_iter_level
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k,&
                                              cp_unit_to_cp2k
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE f77_interface,                   ONLY: f_env_add_defaults,&
                                              f_env_rm_defaults,&
                                              f_env_type
   USE force_env_types,                 ONLY: force_env_get
   USE gle_system_dynamics,             ONLY: gle_cholesky_stab,&
                                              gle_matrix_exp,&
                                              restart_gle
   USE gle_system_types,                ONLY: gle_dealloc,&
                                              gle_init,&
                                              gle_thermo_create
   USE global_types,                    ONLY: global_environment_type
   USE helium_interactions,             ONLY: helium_intpot_scan
   USE helium_io,                       ONLY: helium_write_cubefile
   USE helium_methods,                  ONLY: helium_create,&
                                              helium_init,&
                                              helium_release
   USE helium_sampling,                 ONLY: helium_do_run,&
                                              helium_step
   USE helium_types,                    ONLY: helium_solvent_p_type
   USE input_constants,                 ONLY: integrate_exact,&
                                              integrate_numeric,&
                                              propagator_rpmd,&
                                              transformation_normal,&
                                              transformation_stage
   USE input_cp2k_restarts,             ONLY: write_restart
   USE input_section_types,             ONLY: &
        section_type, section_vals_get, section_vals_get_subs_vals, section_vals_release, &
        section_vals_retain, section_vals_type, section_vals_val_get, section_vals_val_set, &
        section_vals_val_unset
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathconstants,                   ONLY: twopi
   USE mathlib,                         ONLY: gcd
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_comm_self
   USE molecule_kind_list_types,        ONLY: molecule_kind_list_type
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_list_types,             ONLY: molecule_list_type
   USE molecule_types,                  ONLY: global_constraint_type,&
                                              molecule_type
   USE parallel_rng_types,              ONLY: GAUSSIAN,&
                                              rng_stream_type
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE pint_gle,                        ONLY: pint_calc_gle_energy,&
                                              pint_gle_init,&
                                              pint_gle_step
   USE pint_io,                         ONLY: pint_write_action,&
                                              pint_write_centroids,&
                                              pint_write_com,&
                                              pint_write_ener,&
                                              pint_write_line,&
                                              pint_write_rgyr,&
                                              pint_write_step_info,&
                                              pint_write_trajectory
   USE pint_normalmode,                 ONLY: normalmode_calc_uf_h,&
                                              normalmode_env_create,&
                                              normalmode_init_masses,&
                                              normalmode_release
   USE pint_piglet,                     ONLY: pint_calc_piglet_energy,&
                                              pint_piglet_create,&
                                              pint_piglet_init,&
                                              pint_piglet_release,&
                                              pint_piglet_step
   USE pint_pile,                       ONLY: pint_calc_pile_energy,&
                                              pint_pile_init,&
                                              pint_pile_release,&
                                              pint_pile_step
   USE pint_public,                     ONLY: pint_levy_walk
   USE pint_qtb,                        ONLY: pint_calc_qtb_energy,&
                                              pint_qtb_init,&
                                              pint_qtb_release,&
                                              pint_qtb_step
   USE pint_staging,                    ONLY: staging_calc_uf_h,&
                                              staging_env_create,&
                                              staging_init_masses,&
                                              staging_release
   USE pint_transformations,            ONLY: pint_f2uf,&
                                              pint_u2x,&
                                              pint_x2u
   USE pint_types,                      ONLY: &
        e_conserved_id, e_kin_thermo_id, e_kin_virial_id, e_potential_id, pint_env_type, &
        thermostat_gle, thermostat_none, thermostat_nose, thermostat_piglet, thermostat_pile, &
        thermostat_qtb
   USE replica_methods,                 ONLY: rep_env_calc_e_f,&
                                              rep_env_create
   USE replica_types,                   ONLY: rep_env_release,&
                                              replica_env_type
   USE simpar_types,                    ONLY: create_simpar_type,&
                                              release_simpar_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_methods'
   INTEGER, SAVE, PRIVATE :: last_pint_id = 0

   PUBLIC :: do_pint_run

CONTAINS

! ***************************************************************************
!> \brief  Create a path integral environment
!> \param pint_env ...
!> \param input ...
!> \param input_declaration ...
!> \param para_env ...
!> \par    History
!>           Fixed some bugs [hforbert]
!>           Added normal mode transformation [hforbert]
!>           10.2015 Added RPMD propagator and harmonic integrator [Felix Uhl]
!>           10.2018 Added centroid constraints [cschran+rperez]
!> \author fawzi
!> \note   Might return an unassociated pointer in parallel on the processors
!>         that are not needed.
! **************************************************************************************************
   SUBROUTINE pint_create(pint_env, input, input_declaration, para_env)

      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pint_create'

      CHARACTER(len=2*default_string_length)             :: msg
      CHARACTER(len=default_path_length)                 :: output_file_name, project_name
      INTEGER                                            :: handle, iat, ibead, icont, idim, idir, &
                                                            ierr, ig, itmp, nrep, prep, stat
      LOGICAL                                            :: explicit, ltmp
      REAL(kind=dp)                                      :: dt, mass, omega
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(f_env_type), POINTER                          :: f_env
      TYPE(global_constraint_type), POINTER              :: gci
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(replica_env_type), POINTER                    :: rep_env
      TYPE(section_vals_type), POINTER :: constraint_section, gle_section, nose_section, &
         piglet_section, pile_section, pint_section, qtb_section, transform_section

      CALL timeset(routineN, handle)

      NULLIFY (f_env, subsys, particles, nose_section, gle_section, gci)

      CPASSERT(.NOT. ASSOCIATED(pint_env))
      CPASSERT(ASSOCIATED(input))
      CPASSERT(input%ref_count > 0)
      NULLIFY (rep_env)
      pint_section => section_vals_get_subs_vals(input, "MOTION%PINT")
      CALL section_vals_val_get(pint_section, "p", i_val=nrep)
      CALL section_vals_val_get(pint_section, "proc_per_replica", &
                                i_val=prep)
      ! Maybe let the user have his/her way as long as prep is
      ! within the bounds of number of CPUs??
      IF ((prep < 1) .OR. (prep > para_env%num_pe) .OR. &
          (MOD(prep*nrep, para_env%num_pe) /= 0)) THEN
         prep = para_env%num_pe/gcd(para_env%num_pe, nrep)
         IF (para_env%ionode) THEN
            WRITE (UNIT=msg, FMT=*) "PINT WARNING: Adjusting number of processors per replica to ", prep
            CPWARN(msg)
         END IF
      END IF

      ! replica_env modifies the global input structure - which is wrong - one
      ! of the side effects is the inifite adding of the -r-N string to the
      ! project name and the output file name, which corrupts restart files.
      ! For now: save the project name and output file name and restore them
      ! after the rep_env_create has executed - the initialization of the
      ! replicas will run correctly anyways.
      ! TODO: modify rep_env so that it behaves better
      CALL section_vals_val_get(input, "GLOBAL%PROJECT_NAME", c_val=project_name)
      CALL section_vals_val_get(input, "GLOBAL%OUTPUT_FILE_NAME", c_val=output_file_name)
      CALL rep_env_create(rep_env, para_env=para_env, input=input, &
                          input_declaration=input_declaration, nrep=nrep, prep=prep, row_force=.TRUE.)
      CALL section_vals_val_set(input, "GLOBAL%PROJECT_NAME", c_val=TRIM(project_name))
      IF (LEN_TRIM(output_file_name) .GT. 0) THEN
         CALL section_vals_val_set(input, "GLOBAL%OUTPUT_FILE_NAME", c_val=TRIM(output_file_name))
      ELSE
         CALL section_vals_val_unset(input, "GLOBAL%OUTPUT_FILE_NAME")
      END IF
      IF (.NOT. ASSOCIATED(rep_env)) RETURN

      ALLOCATE (pint_env)

      NULLIFY (pint_env%logger)
      pint_env%logger => cp_get_default_logger()
      CALL cp_add_iter_level(pint_env%logger%iter_info, "PINT")

      last_pint_id = last_pint_id + 1
      pint_env%id_nr = last_pint_id
      pint_env%ref_count = 1
      NULLIFY (pint_env%replicas, pint_env%input, pint_env%staging_env, &
               pint_env%normalmode_env, pint_env%propagator)
      pint_env%p = nrep
      pint_env%replicas => rep_env
      pint_env%ndim = rep_env%ndim
      pint_env%input => input

      CALL section_vals_retain(pint_env%input)

      ! get first step, last step, number of steps, etc
      CALL section_vals_val_get(input, "MOTION%PINT%ITERATION", &
                                i_val=itmp)
      pint_env%first_step = itmp
      CALL section_vals_val_get(input, "MOTION%PINT%MAX_STEP", &
                                explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(input, "MOTION%PINT%MAX_STEP", &
                                   i_val=itmp)
         pint_env%last_step = itmp
         pint_env%num_steps = pint_env%last_step - pint_env%first_step
      ELSE
         CALL section_vals_val_get(input, "MOTION%PINT%NUM_STEPS", &
                                   i_val=itmp)
         pint_env%num_steps = itmp
         pint_env%last_step = pint_env%first_step + pint_env%num_steps
      END IF

      CALL section_vals_val_get(pint_section, "DT", &
                                r_val=pint_env%dt)
      pint_env%t = pint_env%first_step*pint_env%dt

      CALL section_vals_val_get(pint_section, "nrespa", i_val=pint_env%nrespa)
      CALL section_vals_val_get(pint_section, "Temp", r_val=pint_env%kT)
      CALL section_vals_val_get(pint_section, "T_TOL", &
                                r_val=pint_env%t_tol)

      CALL section_vals_val_get(pint_section, "HARM_INT", i_val=pint_env%harm_integrator)

      ALLOCATE (pint_env%propagator)
      CALL section_vals_val_get(pint_section, "propagator", &
                                i_val=pint_env%propagator%prop_kind)
      !Initialize simulation temperature depending on the propagator
      !As well as the scaling factor for the physical potential
      IF (pint_env%propagator%prop_kind == propagator_rpmd) THEN
         pint_env%propagator%temp_phys2sim = REAL(pint_env%p, dp)
         pint_env%propagator%physpotscale = 1.0_dp
      ELSE
         pint_env%propagator%temp_phys2sim = 1.0_dp
         pint_env%propagator%physpotscale = 1.0_dp/REAL(pint_env%p, dp)
      END IF
      pint_env%propagator%temp_sim2phys = 1.0_dp/pint_env%propagator%temp_phys2sim
      pint_env%kT = pint_env%kT*pint_env%propagator%temp_phys2sim

      CALL section_vals_val_get(pint_section, "transformation", &
                                i_val=pint_env%transform)

      NULLIFY (pint_env%tx, pint_env%tv, pint_env%tv_t, pint_env%tv_old, pint_env%tv_new, pint_env%tf)

      pint_env%nnos = 0
      pint_env%pimd_thermostat = thermostat_none
      nose_section => section_vals_get_subs_vals(input, "MOTION%PINT%NOSE")
      CALL section_vals_get(nose_section, explicit=explicit)
      IF (explicit) THEN
         IF (pint_env%propagator%prop_kind == propagator_rpmd) THEN
            CPABORT("RPMD Propagator with Nose-thermostat not implemented!")
         END IF
         CALL section_vals_val_get(nose_section, "nnos", i_val=pint_env%nnos)
         IF (pint_env%nnos > 0) THEN
            pint_env%pimd_thermostat = thermostat_nose
            ALLOCATE ( &
               pint_env%tx(pint_env%nnos, pint_env%p, pint_env%ndim), &
               pint_env%tv(pint_env%nnos, pint_env%p, pint_env%ndim), &
               pint_env%tv_t(pint_env%nnos, pint_env%p, pint_env%ndim), &
               pint_env%tv_old(pint_env%nnos, pint_env%p, pint_env%ndim), &
               pint_env%tv_new(pint_env%nnos, pint_env%p, pint_env%ndim), &
               pint_env%tf(pint_env%nnos, pint_env%p, pint_env%ndim))
            pint_env%tx = 0._dp
            pint_env%tv = 0._dp
            pint_env%tv_t = 0._dp
            pint_env%tv_old = 0._dp
            pint_env%tv_new = 0._dp
            pint_env%tf = 0._dp
         END IF
      END IF

      pint_env%beta = 1._dp/(pint_env%kT*pint_env%propagator%temp_sim2phys)
!TODO
! v_tol not in current input structure
! should also probably be part of nose_section
!       CALL section_vals_val_get(transform_section,"v_tol_nose",r_val=pint_env%v_tol)
!MK ... but we have to initialise v_tol
      pint_env%v_tol = 0.0_dp ! to be fixed

      pint_env%randomG = rng_stream_type( &
                         name="pint_randomG", &
                         distribution_type=GAUSSIAN, &
                         extended_precision=.TRUE.)

      ALLOCATE (pint_env%e_pot_bead(pint_env%p))
      pint_env%e_pot_bead = 0._dp
      pint_env%e_pot_h = 0._dp
      pint_env%e_kin_beads = 0._dp
      pint_env%e_pot_t = 0._dp
      pint_env%e_gle = 0._dp
      pint_env%e_pile = 0._dp
      pint_env%e_piglet = 0._dp
      pint_env%e_qtb = 0._dp
      pint_env%e_kin_t = 0._dp
      pint_env%energy(:) = 0.0_dp

!TODO: rearrange to use standard nose hoover chain functions/data types

      ALLOCATE ( &
         pint_env%x(pint_env%p, pint_env%ndim), &
         pint_env%v(pint_env%p, pint_env%ndim), &
         pint_env%f(pint_env%p, pint_env%ndim), &
         pint_env%external_f(pint_env%p, pint_env%ndim), &
         pint_env%ux(pint_env%p, pint_env%ndim), &
         pint_env%ux_t(pint_env%p, pint_env%ndim), &
         pint_env%uv(pint_env%p, pint_env%ndim), &
         pint_env%uv_t(pint_env%p, pint_env%ndim), &
         pint_env%uv_new(pint_env%p, pint_env%ndim), &
         pint_env%uf(pint_env%p, pint_env%ndim), &
         pint_env%uf_h(pint_env%p, pint_env%ndim), &
         pint_env%centroid(pint_env%ndim), &
         pint_env%rtmp_ndim(pint_env%ndim), &
         pint_env%rtmp_natom(pint_env%ndim/3), &
         STAT=stat)
      CPASSERT(stat == 0)
      pint_env%x = 0._dp
      pint_env%v = 0._dp
      pint_env%f = 0._dp
      pint_env%external_f = 0._dp
      pint_env%ux = 0._dp
      pint_env%ux_t = 0._dp
      pint_env%uv = 0._dp
      pint_env%uv_t = 0._dp
      pint_env%uv_new = 0._dp
      pint_env%uf = 0._dp
      pint_env%uf_h = 0._dp
      pint_env%centroid(:) = 0.0_dp
      pint_env%rtmp_ndim = 0._dp
      pint_env%rtmp_natom = 0._dp
      pint_env%time_per_step = 0.0_dp

      IF (pint_env%transform == transformation_stage) THEN
         transform_section => section_vals_get_subs_vals(input, &
                                                         "MOTION%PINT%STAGING")
         CALL staging_env_create(pint_env%staging_env, transform_section, &
                                 p=pint_env%p, kT=pint_env%kT)
      ELSE
         transform_section => section_vals_get_subs_vals(input, &
                                                         "MOTION%PINT%NORMALMODE")
         CALL normalmode_env_create(pint_env%normalmode_env, &
                                    transform_section, p=pint_env%p, kT=pint_env%kT, propagator=pint_env%propagator%prop_kind)
         IF (para_env%ionode) THEN
            IF (pint_env%harm_integrator == integrate_numeric) THEN
               IF (10.0_dp*pint_env%dt/REAL(pint_env%nrespa, dp) > &
                   twopi/(pint_env%p*SQRT(MAXVAL(pint_env%normalmode_env%lambda))* &
                          pint_env%normalmode_env%modefactor)) THEN
                  msg = "PINT WARNING| Number of RESPA steps to small "// &
                        "to integrate the harmonic springs."
                  CPWARN(msg)
               END IF
            END IF
         END IF
      END IF
      ALLOCATE (pint_env%mass(pint_env%ndim))
      CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id, &
                              f_env=f_env)
      CALL force_env_get(force_env=f_env%force_env, subsys=subsys)
      CALL cp_subsys_get(subsys, particles=particles, gci=gci)

!TODO length of pint_env%mass is redundant
      idim = 0
      DO iat = 1, pint_env%ndim/3
         CALL get_atomic_kind(particles%els(iat)%atomic_kind, mass=mass)
         DO idir = 1, 3
            idim = idim + 1
            pint_env%mass(idim) = mass
         END DO
      END DO
      CALL f_env_rm_defaults(f_env, ierr)
      CPASSERT(ierr == 0)

      ALLOCATE (pint_env%Q(pint_env%p), &
                pint_env%mass_beads(pint_env%p, pint_env%ndim), &
                pint_env%mass_fict(pint_env%p, pint_env%ndim))
      IF (pint_env%transform == transformation_stage) THEN
         CALL staging_init_masses(pint_env%staging_env, mass=pint_env%mass, &
                                  mass_beads=pint_env%mass_beads, mass_fict=pint_env%mass_fict, &
                                  Q=pint_env%Q)
      ELSE
         CALL normalmode_init_masses(pint_env%normalmode_env, &
                                     mass=pint_env%mass, mass_beads=pint_env%mass_beads, &
                                     mass_fict=pint_env%mass_fict, Q=pint_env%Q)
      END IF

      NULLIFY (pint_env%gle)
      gle_section => section_vals_get_subs_vals(input, "MOTION%PINT%GLE")
      CALL section_vals_get(gle_section, explicit=explicit)
      IF (explicit) THEN
         ALLOCATE (pint_env%gle)
         CALL gle_init(pint_env%gle, dt=pint_env%dt/pint_env%nrespa, temp=pint_env%kT, &
                       section=gle_section)
         IF (pint_env%pimd_thermostat == thermostat_none .AND. pint_env%gle%ndim .GT. 0) THEN
            pint_env%pimd_thermostat = thermostat_gle

            ! initialize a GLE with ALL degrees of freedom on node 0,
            ! as it seems to me that here everything but force eval is replicated
            pint_env%gle%loc_num_gle = pint_env%p*pint_env%ndim
            pint_env%gle%glob_num_gle = pint_env%gle%loc_num_gle
            ALLOCATE (pint_env%gle%map_info%index(pint_env%gle%loc_num_gle))
            CPASSERT(stat == 0)
            DO itmp = 1, pint_env%gle%loc_num_gle
               pint_env%gle%map_info%index(itmp) = itmp
            ENDDO
            CALL gle_thermo_create(pint_env%gle, pint_env%gle%loc_num_gle)

            ! here we should have read a_mat and c_mat;
            !we can therefore compute the matrices needed for the propagator
            ! deterministic part of the propagator
            CALL gle_matrix_exp((-pint_env%dt/pint_env%nrespa*0.5_dp)*pint_env%gle%a_mat, &
                                pint_env%gle%ndim, 15, 15, pint_env%gle%gle_t)
            ! stochastic part
            CALL gle_cholesky_stab(pint_env%gle%c_mat - MATMUL(pint_env%gle%gle_t, &
                                                               MATMUL(pint_env%gle%c_mat, TRANSPOSE(pint_env%gle%gle_t))), &
                                   pint_env%gle%gle_s, pint_env%gle%ndim)
            ! and initialize the additional momenta
            CALL pint_gle_init(pint_env)
         END IF
      END IF

      !Setup pile thermostat
      NULLIFY (pint_env%pile_therm)
      pile_section => section_vals_get_subs_vals(input, "MOTION%PINT%PILE")
      CALL section_vals_get(pile_section, explicit=explicit)
      IF (explicit) THEN
         CALL cite_reference(Ceriotti2010)
         CALL section_vals_val_get(pint_env%input, &
                                   "MOTION%PINT%INIT%THERMOSTAT_SEED", &
                                   i_val=pint_env%thermostat_rng_seed)
         IF (pint_env%pimd_thermostat == thermostat_none) THEN
            pint_env%pimd_thermostat = thermostat_pile
            CALL pint_pile_init(pile_therm=pint_env%pile_therm, &
                                pint_env=pint_env, &
                                normalmode_env=pint_env%normalmode_env, &
                                section=pile_section)
         ELSE
            CPABORT("PILE thermostat can't be used with another thermostat.")
         END IF
      END IF

      !Setup PIGLET thermostat
      NULLIFY (pint_env%piglet_therm)
      piglet_section => section_vals_get_subs_vals(input, "MOTION%PINT%PIGLET")
      CALL section_vals_get(piglet_section, explicit=explicit)
      IF (explicit) THEN
         CALL cite_reference(Ceriotti2012)
         CALL section_vals_val_get(pint_env%input, &
                                   "MOTION%PINT%INIT%THERMOSTAT_SEED", &
                                   i_val=pint_env%thermostat_rng_seed)
         IF (pint_env%pimd_thermostat == thermostat_none) THEN
            pint_env%pimd_thermostat = thermostat_piglet
            CALL pint_piglet_create(pint_env%piglet_therm, &
                                    pint_env, &
                                    piglet_section)
            CALL pint_piglet_init(pint_env%piglet_therm, &
                                  pint_env, &
                                  piglet_section, &
                                  dt=pint_env%dt, para_env=para_env)
         ELSE
            CPABORT("PILE thermostat can't be used with another thermostat.")
         END IF
      END IF

      !Setup qtb thermostat
      NULLIFY (pint_env%qtb_therm)
      qtb_section => section_vals_get_subs_vals(input, "MOTION%PINT%QTB")
      CALL section_vals_get(qtb_section, explicit=explicit)
      IF (explicit) THEN
         CALL cite_reference(Brieuc2016)
         CALL section_vals_val_get(pint_env%input, &
                                   "MOTION%PINT%INIT%THERMOSTAT_SEED", &
                                   i_val=pint_env%thermostat_rng_seed)
         IF (pint_env%pimd_thermostat == thermostat_none) THEN
            pint_env%pimd_thermostat = thermostat_qtb
            CALL pint_qtb_init(qtb_therm=pint_env%qtb_therm, &
                               pint_env=pint_env, &
                               normalmode_env=pint_env%normalmode_env, &
                               section=qtb_section)
         ELSE
            CPABORT("QTB thermostat can't be used with another thermostat.")
         END IF
      END IF

      !Initialize integrator scheme
      CALL section_vals_val_get(pint_section, "HARM_INT", i_val=pint_env%harm_integrator)
      IF (pint_env%harm_integrator == integrate_exact) THEN
         IF (pint_env%pimd_thermostat == thermostat_nose) THEN
            WRITE (UNIT=msg, FMT=*) "PINT WARNING| Nose Thermostat only available in "// &
               "the numeric harmonic integrator. Switching to numeric harmonic integrator."
            CPWARN(msg)
            pint_env%harm_integrator = integrate_numeric
         END IF
         IF (pint_env%pimd_thermostat == thermostat_gle) THEN
            WRITE (UNIT=msg, FMT=*) "PINT WARNING| GLE Thermostat only available in "// &
               "the numeric harmonic integrator. Switching to numeric harmonic integrator."
            CPWARN(msg)
            pint_env%harm_integrator = integrate_numeric
         END IF
      ELSE IF (pint_env%harm_integrator == integrate_numeric) THEN
         IF (pint_env%pimd_thermostat == thermostat_pile) THEN
            WRITE (UNIT=msg, FMT=*) "PINT WARNING| PILE Thermostat only available in "// &
               "the exact harmonic integrator. Switching to exact harmonic integrator."
            CPWARN(msg)
            pint_env%harm_integrator = integrate_exact
         END IF
         IF (pint_env%pimd_thermostat == thermostat_piglet) THEN
            WRITE (UNIT=msg, FMT=*) "PINT WARNING| PIGLET Thermostat only available in "// &
               "the exact harmonic integrator. Switching to exact harmonic integrator."
            CPWARN(msg)
            pint_env%harm_integrator = integrate_exact
         END IF
         IF (pint_env%pimd_thermostat == thermostat_qtb) THEN
            WRITE (UNIT=msg, FMT=*) "PINT WARNING| QTB Thermostat only available in "// &
               "the exact harmonic integrator. Switching to exact harmonic integrator."
            CPWARN(msg)
            pint_env%harm_integrator = integrate_exact
         END IF
      END IF

      IF (pint_env%harm_integrator == integrate_exact) THEN
         IF (pint_env%nrespa /= 1) THEN
            pint_env%nrespa = 1
            WRITE (UNIT=msg, FMT=*) "PINT WARNING| Adjusting NRESPA to 1 for exact harmonic integration."
            CPWARN(msg)
         END IF
         NULLIFY (pint_env%wsinex)
         ALLOCATE (pint_env%wsinex(pint_env%p))
         NULLIFY (pint_env%iwsinex)
         ALLOCATE (pint_env%iwsinex(pint_env%p))
         NULLIFY (pint_env%cosex)
         ALLOCATE (pint_env%cosex(pint_env%p))
         dt = pint_env%dt/REAL(pint_env%nrespa, KIND=dp)
         !Centroid mode shoud not be propagated
         pint_env%wsinex(1) = 0.0_dp
         pint_env%iwsinex(1) = dt
         pint_env%cosex(1) = 1.0_dp
         DO ibead = 2, pint_env%p
            omega = SQRT(pint_env%normalmode_env%lambda(ibead))
            pint_env%wsinex(ibead) = SIN(omega*dt)*omega
            pint_env%iwsinex(ibead) = SIN(omega*dt)/omega
            pint_env%cosex(ibead) = COS(omega*dt)
         END DO
      END IF

      CALL section_vals_val_get(pint_section, "FIX_CENTROID_POS", &
                                l_val=ltmp)
      IF (ltmp .AND. (pint_env%transform .EQ. transformation_normal)) THEN
         pint_env%first_propagated_mode = 2
      ELSE
         pint_env%first_propagated_mode = 1
      END IF

      ! Constraint information:
      NULLIFY (pint_env%simpar)
      constraint_section => section_vals_get_subs_vals(pint_env%input, &
                                                       "MOTION%CONSTRAINT")
      CALL section_vals_get(constraint_section, explicit=explicit)
      CALL create_simpar_type(pint_env%simpar)
      pint_env%simpar%constraint = explicit
      pint_env%kTcorr = 1.0_dp
      IF (explicit) THEN
         ! Staging not supported yet, since lowest mode is assumed
         ! to be related to centroid
         IF (pint_env%transform == transformation_stage) THEN
            CPABORT("Centroid constraints are not supported for staging transformation")
         END IF

         ! Check thermostats that are not supported:
         IF (pint_env%pimd_thermostat == thermostat_gle) THEN
            WRITE (UNIT=msg, FMT=*) "GLE Thermostat not supported for "// &
               "centroid constraints. Switch to NOSE for numeric integration."
            CPABORT(msg)
         END IF
         ! Warn for NOSE
         IF (pint_env%pimd_thermostat == thermostat_nose) THEN
            WRITE (UNIT=msg, FMT=*) "PINT WARNING| Nose Thermostat set to "// &
               "zero for constrained atoms. Careful interpretation of temperature."
            CPWARN(msg)
            WRITE (UNIT=msg, FMT=*) "PINT WARNING| Lagrange multipliers are "// &
               "are printed every RESPA step and need to be treated carefully."
            CPWARN(msg)
         END IF

         CALL section_vals_val_get(constraint_section, "SHAKE_TOLERANCE", &
                                   r_val=pint_env%simpar%shake_tol)
         pint_env%simpar%info_constraint = cp_print_key_unit_nr(pint_env%logger, &
                                                                constraint_section, &
                                                                "CONSTRAINT_INFO", &
                                                                extension=".shakeLog", &
                                                                log_filename=.FALSE.)
         pint_env%simpar%lagrange_multipliers = cp_print_key_unit_nr(pint_env%logger, &
                                                                     constraint_section, &
                                                                     "LAGRANGE_MULTIPLIERS", &
                                                                     extension=".LagrangeMultLog", &
                                                                     log_filename=.FALSE.)
         pint_env%simpar%dump_lm = &
            BTEST(cp_print_key_should_output(pint_env%logger%iter_info, &
                                             constraint_section, &
                                             "LAGRANGE_MULTIPLIERS"), cp_p_file)

         ! Determine constrained atoms:
         pint_env%n_atoms_constraints = 0
         DO ig = 1, gci%ncolv%ntot
            ! Double counts, if the same atom is involved in different collective variables
            pint_env%n_atoms_constraints = pint_env%n_atoms_constraints + SIZE(gci%colv_list(ig)%i_atoms)
         END DO

         ALLOCATE (pint_env%atoms_constraints(pint_env%n_atoms_constraints))
         icont = 0
         DO ig = 1, gci%ncolv%ntot
            DO iat = 1, SIZE(gci%colv_list(ig)%i_atoms)
               icont = icont + 1
               pint_env%atoms_constraints(icont) = gci%colv_list(ig)%i_atoms(iat)
            END DO
         END DO

         ! Set the correction to the temperature due to the frozen degrees of freedom in NOSE:
         CALL section_vals_val_get(pint_section, "kT_CORRECTION", &
                                   l_val=ltmp)
         IF (ltmp) THEN
            pint_env%kTcorr = 1.0_dp + REAL(3*pint_env%n_atoms_constraints, dp)/(REAL(pint_env%ndim, dp)*REAL(pint_env%p, dp))
         END IF
      END IF

      CALL timestop(handle)

      RETURN
   END SUBROUTINE pint_create

! ***************************************************************************
!> \brief Retain a path integral environment
!> \param pint_env the pint_env to retain
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pint_retain(pint_env)
      TYPE(pint_env_type), POINTER                       :: pint_env

      CPASSERT(ASSOCIATED(pint_env))
      CPASSERT(pint_env%ref_count > 0)
      pint_env%ref_count = pint_env%ref_count + 1
      RETURN
   END SUBROUTINE pint_retain

! ***************************************************************************
!> \brief Release a path integral environment
!> \param pint_env the pint_env to release
!> \par History
!>      Added normal mode transformation [hforbert]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pint_release(pint_env)
      TYPE(pint_env_type), POINTER                       :: pint_env

      IF (ASSOCIATED(pint_env)) THEN
         CPASSERT(pint_env%ref_count > 0)
         pint_env%ref_count = pint_env%ref_count - 1
         IF (pint_env%ref_count == 0) THEN
            CALL rep_env_release(pint_env%replicas)
            CALL section_vals_release(pint_env%input)
            IF (ASSOCIATED(pint_env%staging_env)) THEN
               CALL staging_release(pint_env%staging_env)
            END IF
            IF (ASSOCIATED(pint_env%normalmode_env)) THEN
               CALL normalmode_release(pint_env%normalmode_env)
            END IF

            DEALLOCATE (pint_env%mass)
            DEALLOCATE (pint_env%e_pot_bead)

            DEALLOCATE (pint_env%x)
            DEALLOCATE (pint_env%v)
            DEALLOCATE (pint_env%f)
            DEALLOCATE (pint_env%external_f)
            DEALLOCATE (pint_env%mass_beads)
            DEALLOCATE (pint_env%mass_fict)
            DEALLOCATE (pint_env%ux)
            DEALLOCATE (pint_env%ux_t)
            DEALLOCATE (pint_env%uv)
            DEALLOCATE (pint_env%uv_t)
            DEALLOCATE (pint_env%uv_new)
            DEALLOCATE (pint_env%uf)
            DEALLOCATE (pint_env%uf_h)
            DEALLOCATE (pint_env%centroid)
            DEALLOCATE (pint_env%rtmp_ndim)
            DEALLOCATE (pint_env%rtmp_natom)
            DEALLOCATE (pint_env%propagator)

            IF (pint_env%simpar%constraint) THEN
               DEALLOCATE (pint_env%atoms_constraints)
            END IF
            CALL release_simpar_type(pint_env%simpar)

            IF (pint_env%harm_integrator == integrate_exact) THEN
               DEALLOCATE (pint_env%wsinex)
               DEALLOCATE (pint_env%iwsinex)
               DEALLOCATE (pint_env%cosex)
            END IF

            SELECT CASE (pint_env%pimd_thermostat)
            CASE (thermostat_nose)
               DEALLOCATE (pint_env%tx)
               DEALLOCATE (pint_env%tv)
               DEALLOCATE (pint_env%tv_t)
               DEALLOCATE (pint_env%tv_old)
               DEALLOCATE (pint_env%tv_new)
               DEALLOCATE (pint_env%tf)
            CASE (thermostat_gle)
               CALL gle_dealloc(pint_env%gle)
            CASE (thermostat_pile)
               CALL pint_pile_release(pint_env%pile_therm)
            CASE (thermostat_piglet)
               CALL pint_piglet_release(pint_env%piglet_therm)
            CASE (thermostat_qtb)
               CALL pint_qtb_release(pint_env%qtb_therm)
            END SELECT

            DEALLOCATE (pint_env%Q)

            DEALLOCATE (pint_env)

         END IF
      END IF

      NULLIFY (pint_env)

      RETURN
   END SUBROUTINE pint_release

! ***************************************************************************
!> \brief Tests the path integral methods
!> \param para_env parallel environment
!> \param input the input to test
!> \param input_declaration ...
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE pint_test(para_env, input, input_declaration)
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(section_type), POINTER                        :: input_declaration

      INTEGER                                            :: i, ib, idim, unit_nr
      REAL(kind=dp)                                      :: c, e_h, err
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: x1
      TYPE(pint_env_type), POINTER                       :: pint_env

      CPASSERT(ASSOCIATED(para_env))
      CPASSERT(ASSOCIATED(input))
      CPASSERT(para_env%ref_count > 0)
      CPASSERT(input%ref_count > 0)
      NULLIFY (pint_env)
      unit_nr = cp_logger_get_default_io_unit()
      CALL pint_create(pint_env, input, input_declaration, para_env)
      IF (ASSOCIATED(pint_env)) THEN
         ALLOCATE (x1(pint_env%ndim, pint_env%p))
         x1(:, :) = pint_env%x
         CALL pint_x2u(pint_env)
         pint_env%x = 0._dp
         CALL pint_u2x(pint_env)
         err = 0._dp
         DO i = 1, pint_env%ndim
            err = MAX(err, ABS(x1(1, i) - pint_env%x(1, i)))
         END DO
         IF (unit_nr > 0) WRITE (unit_nr, *) "diff_r1="//cp_to_string(err)

         CALL pint_calc_uf_h(pint_env, e_h=e_h)
         c = -pint_env%staging_env%w_p**2
         pint_env%f = 0._dp
         DO idim = 1, pint_env%ndim
            DO ib = 1, pint_env%p
               pint_env%f(ib, idim) = pint_env%f(ib, idim) + &
                                      c*(2._dp*pint_env%x(ib, idim) &
                                         - pint_env%x(MODULO(ib - 2, pint_env%p) + 1, idim) &
                                         - pint_env%x(MODULO(ib, pint_env%p) + 1, idim))
            END DO
         END DO
         CALL pint_f2uf(pint_env)
         err = 0._dp
         DO idim = 1, pint_env%ndim
            DO ib = 1, pint_env%p
               err = MAX(err, ABS(pint_env%uf(ib, idim) - pint_env%uf_h(ib, idim)))
            END DO
         END DO
         IF (unit_nr > 0) WRITE (unit_nr, *) "diff_f_h="//cp_to_string(err)
      END IF
      RETURN
   END SUBROUTINE pint_test

! ***************************************************************************
!> \brief  Perform a path integral simulation
!> \param  para_env parallel environment
!> \param  input the input to test
!> \param input_declaration ...
!> \param globenv ...
!> \par    History
!>         2003-11 created [fawzi]
!>         2009-12-14 globenv parameter added to handle soft exit
!>           requests [lwalewski]
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE do_pint_run(para_env, input, input_declaration, globenv)
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(section_type), POINTER                        :: input_declaration
      TYPE(global_environment_type), POINTER             :: globenv

      CHARACTER(len=*), PARAMETER                        :: routineN = 'do_pint_run'
      INTEGER, PARAMETER                                 :: helium_only_mid = 1, &
                                                            int_pot_scan_mid = 4, &
                                                            solute_only_mid = 2, &
                                                            solute_with_helium_mid = 3

      CHARACTER(len=default_string_length)               :: stmp
      INTEGER                                            :: handle, mode
      LOGICAL                                            :: explicit, helium_only, int_pot_scan, &
                                                            solvent_present
      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(section_vals_type), POINTER                   :: helium_section

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(para_env))
      CPASSERT(ASSOCIATED(input))
      CPASSERT(para_env%ref_count > 0)
      CPASSERT(input%ref_count > 0)

      ! check if helium solvent is present
      NULLIFY (helium_section)
      helium_section => section_vals_get_subs_vals(input, &
                                                   "MOTION%PINT%HELIUM")
      CALL section_vals_get(helium_section, explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(helium_section, "_SECTION_PARAMETERS_", &
                                   l_val=solvent_present)
      ELSE
         solvent_present = .FALSE.
      END IF

      ! check if there is anything but helium
      IF (solvent_present) THEN
         CALL section_vals_val_get(helium_section, "HELIUM_ONLY", &
                                   l_val=helium_only)
      ELSE
         helium_only = .FALSE.
      END IF

      ! check wheather to perform solute-helium interaction pot scan
      IF (solvent_present) THEN
         CALL section_vals_val_get(helium_section, "INTERACTION_POT_SCAN", &
                                   l_val=int_pot_scan)
      ELSE
         int_pot_scan = .FALSE.
      END IF

      ! input consistency check
      IF (helium_only .AND. int_pot_scan) THEN
         stmp = "Options HELIUM_ONLY and INTERACTION_POT_SCAN are exclusive"
         CPABORT(stmp)
      END IF

      ! select mode of operation
      mode = 0
      IF (solvent_present) THEN
         IF (helium_only) THEN
            mode = helium_only_mid
         ELSE
            IF (int_pot_scan) THEN
               mode = int_pot_scan_mid
            ELSE
               mode = solute_with_helium_mid
            END IF
         END IF
      ELSE
         mode = solute_only_mid
      END IF

      NULLIFY (pint_env)

      ! perform the simulation according to the chosen mode
      SELECT CASE (mode)

      CASE (helium_only_mid)
         CALL helium_create(helium_env, input)
         CALL helium_init(helium_env, pint_env)
         CALL helium_do_run(helium_env, globenv)
         CALL helium_release(helium_env)

      CASE (solute_only_mid)
         CALL pint_create(pint_env, input, input_declaration, para_env)
         CALL pint_init(pint_env)
         CALL pint_do_run(pint_env, globenv)
         CALL pint_release(pint_env)

      CASE (int_pot_scan_mid)
         CALL pint_create(pint_env, input, input_declaration, para_env)
! TODO only initialization of positions is necessary, but rep_env_calc_e_f called
! from within pint_init_f does something to the replica environments which can not be
! avoided (has something to do with f_env_add_defaults) so leaving for now..
         CALL pint_init(pint_env)
         CALL helium_create(helium_env, input, solute=pint_env)
         CALL pint_run_scan(pint_env, helium_env)
         CALL helium_release(helium_env)
         CALL pint_release(pint_env)

      CASE (solute_with_helium_mid)
         CALL pint_create(pint_env, input, input_declaration, para_env)
         ! init pint without helium forces (they are not yet initialized)
         CALL pint_init(pint_env)
         ! init helium with solute's positions (they are already initialized)
         CALL helium_create(helium_env, input, solute=pint_env)
         CALL helium_init(helium_env, pint_env)
         ! reinit pint forces with helium forces (they are now initialized)
         CALL pint_init_f(pint_env, helium_env=helium_env)

         CALL pint_do_run(pint_env, globenv, helium_env=helium_env)
         CALL helium_release(helium_env)
         CALL pint_release(pint_env)

      CASE DEFAULT
         CPABORT("Unknown mode ("//TRIM(ADJUSTL(cp_to_string(mode)))//")")
      END SELECT

      CALL timestop(handle)

      RETURN
   END SUBROUTINE do_pint_run

! ***************************************************************************
!> \brief  Reads the restart, initializes the beads, etc.
!> \param pint_env ...
!> \par    History
!>           11.2003 created [fawzi]
!>           actually ASSIGN input pointer [hforbert]
!>           2010-12-16 turned into a wrapper routine [lwalewski]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pint_init(pint_env)

      TYPE(pint_env_type), POINTER                       :: pint_env

      CPASSERT(ASSOCIATED(pint_env))
      CPASSERT(pint_env%ref_count > 0)

      CALL pint_init_x(pint_env)
      CALL pint_init_v(pint_env)
      CALL pint_init_t(pint_env)
      CALL pint_init_f(pint_env)

      RETURN
   END SUBROUTINE pint_init

! ***************************************************************************
!> \brief  Assign initial postions to the beads.
!> \param pint_env ...
!> \date   2010-12-15
!> \author Lukasz Walewski
!> \note  Initialization is done in the following way:
!>           1. assign all beads with the same classical positions from
!>              FORCE_EVAL (hot start)
!>           2. spread the beads around classical positions as if they were
!>              free particles (if requested)
!>           3. replace positions generated in steps 1-2 with the explicit
!>              ones if they are explicitly given in the input structure
!>           4. apply Gaussian noise to the positions generated so far (if
!>              requested)
! **************************************************************************************************
   SUBROUTINE pint_init_x(pint_env)

      TYPE(pint_env_type), POINTER                       :: pint_env

      CHARACTER(len=default_string_length)               :: msg, tmp
      INTEGER                                            :: ia, ib, ic, idim, input_seed, n_rep_val
      LOGICAL                                            :: done_init, done_levy, done_rand, &
                                                            explicit, levycorr, ltmp
      REAL(kind=dp)                                      :: tcorr, var
      REAL(kind=dp), DIMENSION(3)                        :: x0
      REAL(kind=dp), DIMENSION(3, 2)                     :: seed
      REAL(kind=dp), DIMENSION(:), POINTER               :: bx, r_vals
      TYPE(rng_stream_type), ALLOCATABLE                 :: rng_gaussian
      TYPE(section_vals_type), POINTER                   :: input_section

      CPASSERT(ASSOCIATED(pint_env))
      CPASSERT(pint_env%ref_count > 0)

      DO idim = 1, pint_env%ndim
         DO ib = 1, pint_env%p
            pint_env%x(ib, idim) = pint_env%replicas%r(idim, ib)
         END DO
      END DO

      done_levy = .FALSE.
      CALL section_vals_val_get(pint_env%input, &
                                "MOTION%PINT%INIT%LEVY_POS_SAMPLE", &
                                l_val=ltmp)
      CALL section_vals_val_get(pint_env%input, &
                                "MOTION%PINT%INIT%LEVY_TEMP_FACTOR", &
                                r_val=tcorr)
      IF (ltmp) THEN

         NULLIFY (bx)
         ALLOCATE (bx(3*pint_env%p))
         CALL section_vals_val_get(pint_env%input, &
                                   "MOTION%PINT%INIT%LEVY_SEED", i_val=input_seed)
         seed(:, :) = REAL(input_seed, KIND=dp)
!      seed(:,:) = next_rng_seed()
         rng_gaussian = rng_stream_type( &
                        name="tmp_rng_gaussian", &
                        distribution_type=GAUSSIAN, &
                        extended_precision=.TRUE., &
                        seed=seed)

         CALL section_vals_val_get(pint_env%input, &
                                   "MOTION%PINT%INIT%LEVY_CORRELATED", &
                                   l_val=levycorr)

         IF (levycorr) THEN

            ! correlated Levy walk - the same path for all atoms
            x0 = (/0.0_dp, 0.0_dp, 0.0_dp/)
            CALL pint_levy_walk(x0, pint_env%p, 1.0_dp, bx, rng_gaussian)
            idim = 0
            DO ia = 1, pint_env%ndim/3
               var = SQRT(1.0_dp/(pint_env%kT*tcorr*pint_env%mass(3*ia)))
               DO ic = 1, 3
                  idim = idim + 1
                  DO ib = 1, pint_env%p
                     pint_env%x(ib, idim) = pint_env%x(ib, idim) + bx(3*(ib - 1) + ic)*var
                  END DO
               END DO
            END DO

         ELSE

            ! uncorrelated bead initialization - distinct Levy walk for each atom
            idim = 0
            DO ia = 1, pint_env%ndim/3
               x0(1) = pint_env%x(1, 3*(ia - 1) + 1)
               x0(2) = pint_env%x(1, 3*(ia - 1) + 2)
               x0(3) = pint_env%x(1, 3*(ia - 1) + 3)
               var = SQRT(1.0_dp/(pint_env%kT*tcorr*pint_env%mass(3*ia)))
               CALL pint_levy_walk(x0, pint_env%p, var, bx, rng_gaussian)
               DO ic = 1, 3
                  idim = idim + 1
                  DO ib = 1, pint_env%p
                     pint_env%x(ib, idim) = pint_env%x(ib, idim) + bx(3*(ib - 1) + ic)
                  END DO
               END DO
            END DO

         END IF

         DEALLOCATE (bx)
         done_levy = .TRUE.
      END IF

      done_init = .FALSE.
      NULLIFY (input_section)
      input_section => section_vals_get_subs_vals(pint_env%input, &
                                                  "MOTION%PINT%BEADS%COORD")
      CALL section_vals_get(input_section, explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(input_section, "_DEFAULT_KEYWORD_", &
                                   n_rep_val=n_rep_val)
         IF (n_rep_val > 0) THEN
            CPASSERT(n_rep_val == 1)
            CALL section_vals_val_get(input_section, "_DEFAULT_KEYWORD_", &
                                      r_vals=r_vals)
            IF (SIZE(r_vals) /= pint_env%p*pint_env%ndim) &
               CPABORT("Invalid size of MOTION%PINT%BEADS%COORD")
            ic = 0
            DO idim = 1, pint_env%ndim
               DO ib = 1, pint_env%p
                  ic = ic + 1
                  pint_env%x(ib, idim) = r_vals(ic)
               END DO
            END DO
            done_init = .TRUE.
         END IF
      END IF

      done_rand = .FALSE.
      CALL section_vals_val_get(pint_env%input, &
                                "MOTION%PINT%INIT%RANDOMIZE_POS", &
                                l_val=ltmp)
      IF (ltmp) THEN
         DO idim = 1, pint_env%ndim
            DO ib = 1, pint_env%p
               pint_env%x(ib, idim) = pint_env%x(ib, idim) + &
                                      pint_env%randomG%next(variance=pint_env%beta/ &
                                                            SQRT(12.0_dp*pint_env%mass(idim)))
            END DO
         END DO
         done_rand = .TRUE.
      END IF

      WRITE (tmp, '(A)') "Bead positions initialization:"
      IF (done_init) THEN
         WRITE (msg, '(A,A)') TRIM(tmp), " input structure"
      ELSE IF (done_levy) THEN
         WRITE (msg, '(A,A)') TRIM(tmp), " Levy random walk"
      ELSE
         WRITE (msg, '(A,A)') TRIM(tmp), " hot start"
      END IF
      CALL pint_write_line(msg)

      IF (done_levy) THEN
         WRITE (msg, '(A,F6.3)') "Levy walk at effective temperature: ", tcorr
      END IF

      IF (done_rand) THEN
         WRITE (msg, '(A)') "Added gaussian noise to the positions of the beads."
         CALL pint_write_line(msg)
      END IF

      RETURN
   END SUBROUTINE pint_init_x

! ***************************************************************************
!> \brief  Initialize velocities
!> \param  pint_env the pint env in which you should initialize the
!>         velocity
!> \par    History
!>         2010-12-16 gathered all velocity-init code here [lwalewski]
!>         2011-04-05 added centroid velocity initialization [lwalewski]
!>         2011-12-19 removed optional parameter kT, target temperature is
!>                    now determined from the input directly [lwalewski]
!> \author fawzi
!> \note   Initialization is done according to the following protocol:
!>         1. set all the velocities to FORCE_EVAL%SUBSYS%VELOCITY if present
!>         2. scale the velocities according to the actual temperature
!>            (has no effect if vels not present in 1.)
!>         3. draw vels for the remaining dof from MB distribution
!>            (all or non-centroid modes only depending on 1.)
!>         4. add random noise to the centroid vels if CENTROID_SPEED == T
!>         5. set the vels for all dof to 0.0 if VELOCITY_QUENCH == T
!>         6. set the vels according to the explicit values from the input
!>            if present
! **************************************************************************************************
   SUBROUTINE pint_init_v(pint_env)
      TYPE(pint_env_type), POINTER                       :: pint_env

      CHARACTER(len=default_string_length)               :: msg, stmp, stmp1, stmp2, unit_str
      INTEGER                                            :: first_mode, i, ia, ib, ic, idim, ierr, &
                                                            itmp, j, n_rep_val, nparticle, &
                                                            nparticle_kind
      LOGICAL                                            :: done_init, done_quench, done_scale, &
                                                            done_sped, explicit, ltmp, vels_present
      REAL(kind=dp)                                      :: actual_t, ek, factor, rtmp, target_t, &
                                                            unit_conv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: vel
      REAL(kind=dp), DIMENSION(:), POINTER               :: r_vals
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(distribution_1d_type), POINTER                :: local_molecules, local_particles
      TYPE(f_env_type), POINTER                          :: f_env
      TYPE(global_constraint_type), POINTER              :: gci
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: input_section

      CPASSERT(ASSOCIATED(pint_env))
      CPASSERT(pint_env%ref_count > 0)

      NULLIFY (logger)
      logger => cp_get_default_logger()

      ! Get centroid constraint info, if needed
      ! Create a force environment which will be identical to
      ! the bead that is being processed by the processor.
      IF (pint_env%simpar%constraint) THEN
         NULLIFY (subsys, cell)
         NULLIFY (atomic_kinds, local_particles, particles)
         NULLIFY (local_molecules, molecules, molecule_kinds, gci)
         NULLIFY (atomic_kind_set, molecule_kind_set, particle_set, molecule_set)

         CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id, f_env=f_env)
         CALL force_env_get(force_env=f_env%force_env, subsys=subsys)
         CALL f_env_rm_defaults(f_env, ierr)
         CPASSERT(ierr == 0)

         ! Get gci and more from subsys
         CALL cp_subsys_get(subsys=subsys, &
                            cell=cell, &
                            atomic_kinds=atomic_kinds, &
                            local_particles=local_particles, &
                            particles=particles, &
                            local_molecules=local_molecules, &
                            molecules=molecules, &
                            molecule_kinds=molecule_kinds, &
                            gci=gci)

         nparticle_kind = atomic_kinds%n_els
         atomic_kind_set => atomic_kinds%els
         molecule_kind_set => molecule_kinds%els
         nparticle = particles%n_els
         particle_set => particles%els
         molecule_set => molecules%els

         ! Allocate work storage
         ALLOCATE (vel(3, nparticle))
         vel(:, :) = 0.0_dp
         CALL getold(gci, local_molecules, molecule_set, &
                     molecule_kind_set, particle_set, cell)
      END IF

      ! read the velocities from the input file if they are given explicitly
      vels_present = .FALSE.
      NULLIFY (input_section)
      input_section => section_vals_get_subs_vals(pint_env%input, &
                                                  "FORCE_EVAL%SUBSYS%VELOCITY")
      CALL section_vals_get(input_section, explicit=explicit)
      IF (explicit) THEN

         CALL section_vals_val_get(input_section, "PINT_UNIT", &
                                   c_val=unit_str)
         unit_conv = cp_unit_to_cp2k(1.0_dp, TRIM(unit_str))

         ! assign all the beads with the same velocities from FORCE_EVAL%SUBSYS%VELOCITY
         NULLIFY (r_vals)
         CALL section_vals_val_get(input_section, "_DEFAULT_KEYWORD_", &
                                   n_rep_val=n_rep_val)
         stmp = ""
         WRITE (stmp, *) n_rep_val
         msg = "Invalid number of atoms in FORCE_EVAL%SUBSYS%VELOCITY ("// &
               TRIM(ADJUSTL(stmp))//")."
         IF (3*n_rep_val /= pint_env%ndim) &
            CPABORT(msg)
         DO ia = 1, pint_env%ndim/3
            CALL section_vals_val_get(input_section, "_DEFAULT_KEYWORD_", &
                                      i_rep_val=ia, r_vals=r_vals)
            itmp = SIZE(r_vals)
            stmp = ""
            WRITE (stmp, *) itmp
            msg = "Number of coordinates != 3 in FORCE_EVAL%SUBSYS%VELOCITY ("// &
                  TRIM(ADJUSTL(stmp))//")."
            IF (itmp /= 3) &
               CPABORT(msg)
            DO ib = 1, pint_env%p
               DO ic = 1, 3
                  idim = 3*(ia - 1) + ic
                  pint_env%v(ib, idim) = r_vals(ic)*unit_conv
               END DO
            END DO
         END DO

         vels_present = .TRUE.
      END IF

      ! set the actual temperature...
      IF (vels_present) THEN
         ! ...from the initial velocities
         ek = 0.0_dp
         DO ia = 1, pint_env%ndim/3
            rtmp = 0.0_dp
            DO ic = 1, 3
               idim = 3*(ia - 1) + ic
               rtmp = rtmp + pint_env%v(1, idim)*pint_env%v(1, idim)
            END DO
            ek = ek + 0.5_dp*pint_env%mass(idim)*rtmp
         END DO
         actual_t = 2.0_dp*ek/pint_env%ndim
      ELSE
         ! ...using the temperature value from the input
         actual_t = pint_env%kT
      END IF

      ! set the target temperature
      target_t = pint_env%kT
      CALL section_vals_val_get(pint_env%input, &
                                "MOTION%PINT%INIT%VELOCITY_SCALE", &
                                l_val=done_scale)
      IF (vels_present) THEN
         IF (done_scale) THEN
            ! rescale the velocities to match the target temperature
            rtmp = SQRT(target_t/actual_t)
            DO ia = 1, pint_env%ndim/3
               DO ib = 1, pint_env%p
                  DO ic = 1, 3
                     idim = 3*(ia - 1) + ic
                     pint_env%v(ib, idim) = rtmp*pint_env%v(ib, idim)
                  END DO
               END DO
            END DO
         ELSE
            target_t = actual_t
         END IF
      END IF

      ! draw velocities from the M-B distribution...
      IF (vels_present) THEN
         ! ...for non-centroid modes only
         CALL pint_x2u(pint_env, ux=pint_env%uv, x=pint_env%v)
         first_mode = 2
      ELSE
         ! ...for all the modes
         first_mode = 1
      END IF
      DO idim = 1, SIZE(pint_env%uv, 2)
         DO ib = first_mode, SIZE(pint_env%uv, 1)
            pint_env%uv(ib, idim) = &
               pint_env%randomG%next(variance=target_t/pint_env%mass_fict(ib, idim))
         END DO
      END DO

      ! add random component to the centroid velocity if requested
      done_sped = .FALSE.
      CALL section_vals_val_get(pint_env%input, &
                                "MOTION%PINT%INIT%CENTROID_SPEED", &
                                l_val=ltmp)
      IF (ltmp) THEN
         CALL pint_u2x(pint_env, ux=pint_env%uv, x=pint_env%v)
         DO idim = 1, pint_env%ndim
            rtmp = pint_env%randomG%next(variance=pint_env%mass(idim)*pint_env%kT) &
                   /pint_env%mass(idim)
            DO ib = 1, pint_env%p
               pint_env%v(ib, idim) = pint_env%v(ib, idim) + rtmp
            END DO
         END DO
         CALL pint_x2u(pint_env, ux=pint_env%uv, x=pint_env%v)
         done_sped = .TRUE.
      END IF

      ! quench (set to zero) velocities for all the modes if requested
      ! (disregard all the initialization done so far)
      done_quench = .FALSE.
      CALL section_vals_val_get(pint_env%input, &
                                "MOTION%PINT%INIT%VELOCITY_QUENCH", &
                                l_val=ltmp)
      IF (ltmp) THEN
         DO idim = 1, pint_env%ndim
            DO ib = 1, pint_env%p
               pint_env%v(ib, idim) = 0.0_dp
            END DO
         END DO
         CALL pint_x2u(pint_env, ux=pint_env%uv, x=pint_env%v)
         done_quench = .TRUE.
      END IF

      ! set the velocities to the values from the input if they are explicit
      ! (disregard all the initialization done so far)
      done_init = .FALSE.
      NULLIFY (input_section)
      input_section => section_vals_get_subs_vals(pint_env%input, &
                                                  "MOTION%PINT%BEADS%VELOCITY")
      CALL section_vals_get(input_section, explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(input_section, "_DEFAULT_KEYWORD_", &
                                   n_rep_val=n_rep_val)
         IF (n_rep_val > 0) THEN
            CPASSERT(n_rep_val == 1)
            CALL section_vals_val_get(input_section, "_DEFAULT_KEYWORD_", &
                                      r_vals=r_vals)
            IF (SIZE(r_vals) /= pint_env%p*pint_env%ndim) &
               CPABORT("Invalid size of MOTION%PINT%BEAD%VELOCITY")
            itmp = 0
            DO idim = 1, pint_env%ndim
               DO ib = 1, pint_env%p
                  itmp = itmp + 1
                  pint_env%v(ib, idim) = r_vals(itmp)
               END DO
            END DO
            CALL pint_x2u(pint_env, ux=pint_env%uv, x=pint_env%v)
            done_init = .TRUE.
         END IF
      END IF

      unit_conv = cp_unit_from_cp2k(1.0_dp, "K")
      WRITE (stmp1, '(F10.2)') target_t*pint_env%propagator%temp_sim2phys*unit_conv
      msg = "Bead velocities initialization:"
      IF (done_init) THEN
         msg = TRIM(msg)//" input structure"
      ELSE IF (done_quench) THEN
         msg = TRIM(msg)//" quenching (set to 0.0)"
      ELSE
         IF (vels_present) THEN
            msg = TRIM(ADJUSTL(msg))//" centroid +"
         END IF
         msg = TRIM(ADJUSTL(msg))//" Maxwell-Boltzmann at "//TRIM(ADJUSTL(stmp1))//" K."
      END IF
      CALL pint_write_line(msg)

      IF (done_init .AND. done_quench) THEN
         msg = "WARNING: exclusive options requested (velocity restart and quenching)"
         CPWARN(msg)
         msg = "WARNING: velocity restart took precedence"
         CPWARN(msg)
      END IF

      IF ((.NOT. done_init) .AND. (.NOT. done_quench)) THEN
         IF (vels_present .AND. done_scale) THEN
            WRITE (stmp1, '(F10.2)') actual_t*unit_conv
            WRITE (stmp2, '(F10.2)') target_t*unit_conv
            msg = "Scaled initial velocities from "//TRIM(ADJUSTL(stmp1))// &
                  " to "//TRIM(ADJUSTL(stmp2))//" K as requested."
            CPWARN(msg)
         END IF
         IF (done_sped) THEN
            msg = "Added random component to the initial centroid velocities."
            CPWARN(msg)
         END IF
      END IF

      ! Apply constraints to the initial velocities for centroid constraints
      IF (pint_env%simpar%constraint) THEN
         IF (pint_env%propagator%prop_kind == propagator_rpmd) THEN
            ! Multiply with 1/SQRT(n_beads) due to normal mode transformation in RPMD
            factor = SQRT(REAL(pint_env%p, dp))
         ELSE
            ! lowest NM is centroid
            factor = 1.0_dp
         END IF

         IF (pint_env%logger%para_env%ionode) THEN
            DO i = 1, nparticle
               DO j = 1, 3
                  vel(j, i) = pint_env%uv(1, j + (i - 1)*3)/factor
               END DO
            END DO

            ! Possibly update the target values
            CALL shake_update_targets(gci, local_molecules, molecule_set, &
                                      molecule_kind_set, pint_env%dt, &
                                      f_env%force_env%root_section)
            CALL rattle_control(gci, local_molecules, molecule_set, &
                                molecule_kind_set, particle_set, &
                                vel, pint_env%dt, pint_env%simpar%shake_tol, &
                                pint_env%simpar%info_constraint, &
                                pint_env%simpar%lagrange_multipliers, &
                                .FALSE., &
                                cell, mp_comm_self, &
                                local_particles)
         END IF

         ! Broadcast updated velocities to other nodes
         CALL mp_bcast(vel, &
                       pint_env%logger%para_env%source, &
                       pint_env%logger%para_env%group)

         DO i = 1, nparticle
            DO j = 1, 3
               pint_env%uv(1, j + (i - 1)*3) = vel(j, i)*factor
            END DO
         END DO

      END IF

      RETURN
   END SUBROUTINE pint_init_v

! ***************************************************************************
!> \brief  Assign initial postions and velocities to the thermostats.
!> \param pint_env ...
!> \param kT ...
!> \date   2010-12-15
!> \author Lukasz Walewski
!> \note   Extracted from pint_init
! **************************************************************************************************
   SUBROUTINE pint_init_t(pint_env, kT)

      TYPE(pint_env_type), POINTER                       :: pint_env
      REAL(kind=dp), INTENT(in), OPTIONAL                :: kT

      INTEGER                                            :: ib, idim, ii, inos, n_rep_val
      LOGICAL                                            :: explicit, gle_restart
      REAL(kind=dp)                                      :: mykt
      REAL(kind=dp), DIMENSION(:), POINTER               :: r_vals
      TYPE(section_vals_type), POINTER                   :: input_section

      CPASSERT(ASSOCIATED(pint_env))
      CPASSERT(pint_env%ref_count > 0)

      IF (pint_env%pimd_thermostat == thermostat_nose) THEN

         mykt = pint_env%kT
         IF (PRESENT(kT)) mykt = kT
         DO idim = 1, SIZE(pint_env%tv, 3)
            DO ib = 1, SIZE(pint_env%tv, 2)
               DO inos = 1, SIZE(pint_env%tv, 1)
                  pint_env%tv(inos, ib, idim) = &
                     pint_env%randomG%next(variance=mykt/pint_env%Q(ib))
               END DO
            END DO
         END DO

         NULLIFY (input_section)
         input_section => section_vals_get_subs_vals(pint_env%input, &
                                                     "MOTION%PINT%NOSE%COORD")
         CALL section_vals_get(input_section, explicit=explicit)
         IF (explicit) THEN
            CALL section_vals_val_get(input_section, "_DEFAULT_KEYWORD_", &
                                      n_rep_val=n_rep_val)
            IF (n_rep_val > 0) THEN
               CPASSERT(n_rep_val == 1)
               CALL section_vals_val_get(input_section, "_DEFAULT_KEYWORD_", &
                                         r_vals=r_vals)
               IF (SIZE(r_vals) /= pint_env%p*pint_env%ndim*pint_env%nnos) &
                  CPABORT("Invalid size of MOTION%PINT%NOSE%COORD")
               ii = 0
               DO idim = 1, pint_env%ndim
                  DO ib = 1, pint_env%p
                     DO inos = 1, pint_env%nnos
                        ii = ii + 1
                        pint_env%tx(inos, ib, idim) = r_vals(ii)
                     END DO
                  END DO
               END DO
            END IF
         END IF

         NULLIFY (input_section)
         input_section => section_vals_get_subs_vals(pint_env%input, &
                                                     "MOTION%PINT%NOSE%VELOCITY")
         CALL section_vals_get(input_section, explicit=explicit)
         IF (explicit) THEN
            CALL section_vals_val_get(input_section, "_DEFAULT_KEYWORD_", &
                                      n_rep_val=n_rep_val)
            IF (n_rep_val > 0) THEN
               CPASSERT(n_rep_val == 1)
               CALL section_vals_val_get(input_section, "_DEFAULT_KEYWORD_", &
                                         r_vals=r_vals)
               IF (SIZE(r_vals) /= pint_env%p*pint_env%ndim*pint_env%nnos) &
                  CPABORT("Invalid size of MOTION%PINT%NOSE%VELOCITY")
               ii = 0
               DO idim = 1, pint_env%ndim
                  DO ib = 1, pint_env%p
                     DO inos = 1, pint_env%nnos
                        ii = ii + 1
                        pint_env%tv(inos, ib, idim) = r_vals(ii)
                     END DO
                  END DO
               END DO
            END IF
         END IF

      ELSEIF (pint_env%pimd_thermostat == thermostat_gle) THEN
         NULLIFY (input_section)
         input_section => section_vals_get_subs_vals(pint_env%input, &
                                                     "MOTION%PINT%GLE")
         CALL section_vals_get(input_section, explicit=explicit)
         IF (explicit) THEN
            CALL restart_gle(pint_env%gle, input_section, save_mem=.FALSE., &
                             restart=gle_restart)
         END IF
      END IF

      RETURN
   END SUBROUTINE pint_init_t

! ***************************************************************************
!> \brief  Prepares the forces, etc. to perform an PIMD step
!> \param pint_env ...
!> \param helium_env ...
!> \par    History
!>           Added nh_energy calculation [hforbert]
!>           Bug fixes for no thermostats [hforbert]
!>           2016-07-14 Modified to work with independent helium_env [cschran]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE pint_init_f(pint_env, helium_env)
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: helium_env

      INTEGER                                            :: ib, idim, inos
      REAL(kind=dp)                                      :: e_h
      TYPE(cp_logger_type), POINTER                      :: logger

      CPASSERT(ASSOCIATED(pint_env))
      CPASSERT(pint_env%ref_count > 0)

      NULLIFY (logger)
      logger => cp_get_default_logger()

      ! initialize iteration info
      CALL cp_iterate(logger%iter_info, iter_nr=pint_env%first_step)
      CALL cp_iterate(pint_env%logger%iter_info, iter_nr=pint_env%first_step)

      CALL pint_x2u(pint_env)
      CALL pint_calc_uf_h(pint_env=pint_env, e_h=e_h)
      CALL pint_calc_f(pint_env)

      ! add helium forces to the solute's internal ones
      ! Assume that helium has been already initialized and helium_env(1)
      ! contains proper forces in force_avrg array at ionode
      IF (PRESENT(helium_env)) THEN
         IF (logger%para_env%ionode) THEN
            pint_env%f(:, :) = pint_env%f(:, :) + helium_env(1)%helium%force_avrg(:, :)
         END IF
         CALL mp_bcast(pint_env%f, &
                       logger%para_env%source, &
                       logger%para_env%group)
      END IF
      CALL pint_f2uf(pint_env)

      ! set the centroid forces to 0 if FIX_CENTROID_POS
      IF (pint_env%first_propagated_mode .EQ. 2) THEN
         pint_env%uf(1, :) = 0.0_dp
      END IF

      CALL pint_calc_e_kin_beads_u(pint_env)
      CALL pint_calc_e_vir(pint_env)
      DO idim = 1, SIZE(pint_env%uf_h, 2)
         DO ib = pint_env%first_propagated_mode, SIZE(pint_env%uf_h, 1)
            pint_env%uf(ib, idim) = REAL(pint_env%nrespa, dp)*pint_env%uf(ib, idim)
         END DO
      END DO

      IF (pint_env%nnos > 0) THEN
         DO idim = 1, SIZE(pint_env%uf_h, 2)
            DO ib = 1, SIZE(pint_env%uf_h, 1)
               pint_env%tf(1, ib, idim) = (pint_env%mass_fict(ib, idim)* &
                                           pint_env%uv(ib, idim)**2 - pint_env%kT)/pint_env%Q(ib)
            END DO
         END DO

         DO idim = 1, pint_env%ndim
            DO ib = 1, pint_env%p
               DO inos = 1, pint_env%nnos - 1
                  pint_env%tf(inos + 1, ib, idim) = pint_env%tv(inos, ib, idim)**2 - &
                                                    pint_env%kT/pint_env%Q(ib)
               END DO
               DO inos = 1, pint_env%nnos - 1
                  pint_env%tf(inos, ib, idim) = pint_env%tf(inos, ib, idim) &
                                                - pint_env%tv(inos, ib, idim)*pint_env%tv(inos + 1, ib, idim)
               END DO
            END DO
         END DO
         CALL pint_calc_nh_energy(pint_env)
      END IF

      RETURN
   END SUBROUTINE pint_init_f

! ***************************************************************************
!> \brief  Perform the PIMD simulation (main MD loop)
!> \param pint_env ...
!> \param globenv ...
!> \param helium_env ...
!> \par    History
!>         2003-11 created [fawzi]
!>         renamed from pint_run to pint_do_run because of conflicting name
!>           of pint_run in input_constants [hforbert]
!>         2009-12-14 globenv parameter added to handle soft exit
!>           requests [lwalewski]
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Fawzi Mohamed
!> \note   Everything should be read for an md step.
! **************************************************************************************************
   SUBROUTINE pint_do_run(pint_env, globenv, helium_env)
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(helium_solvent_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: helium_env

      INTEGER                                            :: k, step
      LOGICAL                                            :: should_stop
      REAL(kind=dp)                                      :: scal
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(f_env_type), POINTER                          :: f_env

      CPASSERT(ASSOCIATED(pint_env))

      ! initialize iteration info
      CALL cp_iterate(pint_env%logger%iter_info, iter_nr=pint_env%first_step)

      ! iterate replica pint counter by accessing the globally saved
      ! force environment error/logger variables and setting them
      ! explicitly to the pimd "PINT" step value
      CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id, &
                              f_env=f_env)
      NULLIFY (logger)
      logger => cp_get_default_logger()
      CALL cp_iterate(logger%iter_info, &
                      iter_nr=pint_env%first_step)
      CALL f_env_rm_defaults(f_env)

      pint_env%iter = pint_env%first_step

      IF (PRESENT(helium_env)) THEN
         IF (ASSOCIATED(helium_env)) THEN
            ! set the properties accumulated over the whole MC process to 0
            DO k = 1, SIZE(helium_env)
               helium_env(k)%helium%proarea%accu(:) = 0.0_dp
               helium_env(k)%helium%prarea2%accu(:) = 0.0_dp
               helium_env(k)%helium%wnmber2%accu(:) = 0.0_dp
               helium_env(k)%helium%mominer%accu(:) = 0.0_dp
               IF (helium_env(k)%helium%rho_present) THEN
                  helium_env(k)%helium%rho_accu(:, :, :, :) = 0.0_dp
               END IF
               IF (helium_env(k)%helium%rdf_present) THEN
                  helium_env(k)%helium%rdf_accu(:, :) = 0.0_dp
               END IF
            END DO
         END IF
      END IF

      ! write the properties at 0-th step
      CALL pint_calc_energy(pint_env)
      CALL pint_calc_total_action(pint_env)
      CALL pint_write_ener(pint_env)
      CALL pint_write_action(pint_env)
      CALL pint_write_centroids(pint_env)
      CALL pint_write_trajectory(pint_env)
      CALL pint_write_com(pint_env)
      CALL pint_write_rgyr(pint_env)

      ! main PIMD loop
      DO step = 1, pint_env%num_steps

         pint_env%iter = pint_env%iter + 1
         CALL cp_iterate(pint_env%logger%iter_info, &
                         last=(step == pint_env%num_steps), &
                         iter_nr=pint_env%iter)
         CALL cp_iterate(logger%iter_info, &
                         last=(step == pint_env%num_steps), &
                         iter_nr=pint_env%iter)
         pint_env%t = pint_env%t + pint_env%dt

         IF (pint_env%t_tol > 0.0_dp) THEN
            IF (ABS(2._dp*pint_env%e_kin_beads/(pint_env%p*pint_env%ndim) &
                    - pint_env%kT) > pint_env%t_tol) THEN
               scal = SQRT(pint_env%kT*(pint_env%p*pint_env%ndim)/(2.0_dp*pint_env%e_kin_beads))
               pint_env%uv = scal*pint_env%uv
               CALL pint_init_f(pint_env, helium_env=helium_env)
            END IF
         END IF
         CALL pint_step(pint_env, helium_env=helium_env)

         CALL pint_write_ener(pint_env)
         CALL pint_write_action(pint_env)
         CALL pint_write_centroids(pint_env)
         CALL pint_write_trajectory(pint_env)
         CALL pint_write_com(pint_env)
         CALL pint_write_rgyr(pint_env)

         CALL write_restart(root_section=pint_env%input, &
                            pint_env=pint_env, helium_env=helium_env)

         ! exit from the main loop if soft exit has been requested
         CALL external_control(should_stop, "PINT", globenv=globenv)
         IF (should_stop) EXIT

      END DO

      ! remove iteration level
      CALL cp_rm_iter_level(pint_env%logger%iter_info, "PINT")

      RETURN
   END SUBROUTINE pint_do_run

! ***************************************************************************
!> \brief  Performs a scan of the helium-solute interaction energy
!> \param pint_env ...
!> \param helium_env ...
!> \date   2013-11-26
!> \parm   History
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE pint_run_scan(pint_env, helium_env)
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env

      CHARACTER(len=default_string_length)               :: comment
      INTEGER                                            :: unit_nr
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: DATA
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (pint_env%logger, print_key)
      pint_env%logger => cp_get_default_logger()

      ! assume that ionode always has at least one helium_env
      IF (pint_env%logger%para_env%ionode) THEN
         print_key => section_vals_get_subs_vals(helium_env(1)%helium%input, &
                                                 "MOTION%PINT%HELIUM%PRINT%RHO")
      END IF

      ! perform the actual scan wrt the COM of the solute
      CALL helium_intpot_scan(pint_env, helium_env)

      ! output the interaction potential into a cubefile
      ! assume that ionode always has at least one helium_env
      IF (pint_env%logger%para_env%ionode) THEN

         unit_nr = cp_print_key_unit_nr( &
                   pint_env%logger, &
                   print_key, &
                   middle_name="helium-pot", &
                   extension=".cube", &
                   file_position="REWIND", &
                   do_backup=.FALSE.)

         comment = "Solute - helium interaction potential"
         NULLIFY (DATA)
         DATA => helium_env(1)%helium%rho_inst(1, :, :, :)
         CALL helium_write_cubefile( &
            unit_nr, &
            comment, &
            helium_env(1)%helium%center - 0.5_dp* &
            (helium_env(1)%helium%rho_maxr - helium_env(1)%helium%rho_delr), &
            helium_env(1)%helium%rho_delr, &
            helium_env(1)%helium%rho_nbin, &
            DATA)

         CALL m_flush(unit_nr)
         CALL cp_print_key_finished_output(unit_nr, pint_env%logger, print_key)

      END IF

      ! output solute positions
      CALL pint_write_centroids(pint_env)
      CALL pint_write_trajectory(pint_env)

      RETURN
   END SUBROUTINE pint_run_scan

! ***************************************************************************
!> \brief  Does an PINT step (and nrespa harmonic evaluations)
!> \param pint_env ...
!> \param helium_env ...
!> \par    History
!>           various bug fixes [hforbert]
!>           10.2015 Added RPMD propagator and harmonic integrator [Felix Uhl]
!>           04.2016 Changed to work with helium_env [cschran]
!>           10.2018 Added centroid constraints [cschran+rperez]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE pint_step(pint_env, helium_env)
      TYPE(pint_env_type), POINTER                       :: pint_env
      TYPE(helium_solvent_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: helium_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pint_step'

      INTEGER                                            :: handle, i, ia, ib, idim, ierr, inos, &
                                                            iresp, j, k, nparticle, nparticle_kind
      REAL(kind=dp)                                      :: dt_temp, dti, dti2, dti22, e_h, factor, &
                                                            rn, tdti, time_start, time_stop, tol
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: pos, vel
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: tmp
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(distribution_1d_type), POINTER                :: local_molecules, local_particles
      TYPE(f_env_type), POINTER                          :: f_env
      TYPE(global_constraint_type), POINTER              :: gci
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)
      time_start = m_walltime()

      CPASSERT(ASSOCIATED(pint_env))

      rn = REAL(pint_env%nrespa, dp)
      dti = pint_env%dt/rn
      dti2 = dti/2._dp
      tdti = 2.*dti
      dti22 = dti**2/2._dp

      ! Get centroid constraint info, if needed
      ! Create a force environment which will be identical to
      ! the bead that is being processed by the processor.
      IF (pint_env%simpar%constraint) THEN
         NULLIFY (subsys, cell)
         NULLIFY (atomic_kinds, local_particles, particles)
         NULLIFY (local_molecules, molecules, molecule_kinds, gci)
         NULLIFY (atomic_kind_set, molecule_kind_set, particle_set, molecule_set)

         CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id, f_env=f_env)
         CALL force_env_get(force_env=f_env%force_env, subsys=subsys)
         CALL f_env_rm_defaults(f_env, ierr)
         CPASSERT(ierr == 0)

         ! Get gci and more from subsys
         CALL cp_subsys_get(subsys=subsys, &
                            cell=cell, &
                            atomic_kinds=atomic_kinds, &
                            local_particles=local_particles, &
                            particles=particles, &
                            local_molecules=local_molecules, &
                            molecules=molecules, &
                            molecule_kinds=molecule_kinds, &
                            gci=gci)

         nparticle_kind = atomic_kinds%n_els
         atomic_kind_set => atomic_kinds%els
         molecule_kind_set => molecule_kinds%els
         nparticle = particles%n_els
         particle_set => particles%els
         molecule_set => molecules%els

         ! Allocate work storage
         ALLOCATE (pos(3, nparticle))
         ALLOCATE (vel(3, nparticle))
         pos(:, :) = 0.0_dp
         vel(:, :) = 0.0_dp

         IF (pint_env%propagator%prop_kind == propagator_rpmd) THEN
            ! Multiply with 1/SQRT(n_beads) due to normal mode transformation in RPMD
            factor = SQRT(REAL(pint_env%p, dp))
         ELSE
            factor = 1.0_dp
         END IF

         CALL getold(gci, local_molecules, molecule_set, &
                     molecule_kind_set, particle_set, cell)
      END IF

      SELECT CASE (pint_env%harm_integrator)
      CASE (integrate_numeric)

         DO iresp = 1, pint_env%nrespa

            ! integrate bead positions, first_propagated_mode = { 1, 2 }
            ! Nose needs an extra step
            IF (pint_env%pimd_thermostat == thermostat_nose) THEN

               !Set thermostat action of constrained DoF to zero:
               IF (pint_env%simpar%constraint) THEN
                  DO k = 1, pint_env%n_atoms_constraints
                     ia = pint_env%atoms_constraints(k)
                     DO j = 3*(ia - 1) + 1, 3*ia
                        pint_env%tv(:, 1, j) = 0.0_dp
                     END DO
                  END DO
               END IF

               DO i = pint_env%first_propagated_mode, pint_env%p
                  pint_env%ux(i, :) = pint_env%ux(i, :) - &
                                      dti22*pint_env%uv(i, :)*pint_env%tv(1, i, :)
               END DO
               pint_env%tx = pint_env%tx + dti*pint_env%tv + dti22*pint_env%tf

            END IF
            !Integrate position in harmonic springs (uf_h) and physical potential
            !(uf)
            DO i = pint_env%first_propagated_mode, pint_env%p
               pint_env%ux_t(i, :) = pint_env%ux(i, :) + &
                                     dti*pint_env%uv(i, :) + &
                                     dti22*(pint_env%uf_h(i, :) + &
                                            pint_env%uf(i, :))
            END DO

            ! apply thermostats to velocities
            SELECT CASE (pint_env%pimd_thermostat)
            CASE (thermostat_nose)
               pint_env%uv_t = pint_env%uv - dti2* &
                               pint_env%uv*pint_env%tv(1, :, :)
               tmp => pint_env%tv_t
               pint_env%tv_t => pint_env%tv
               pint_env%tv => tmp
               pint_env%tv = pint_env%tv_old + tdti*pint_env%tf
               pint_env%tv_old = pint_env%tv_t
               pint_env%tv_t = pint_env%tv_t + dti2*pint_env%tf
            CASE DEFAULT
               pint_env%uv_t = pint_env%uv
            END SELECT

            !Set thermostat action of constrained DoF to zero:
            IF (pint_env%simpar%constraint) THEN
               DO k = 1, pint_env%n_atoms_constraints
                  ia = pint_env%atoms_constraints(k)
                  DO j = 3*(ia - 1) + 1, 3*ia
                     pint_env%tv(:, 1, j) = 0.0_dp
                     pint_env%tv_t(:, 1, j) = 0.0_dp
                  END DO
               END DO
            END IF

            !Integrate harmonic velocities and physical velocities
            pint_env%uv_t = pint_env%uv_t + dti2*(pint_env%uf_h + pint_env%uf)

            ! physical forces are only applied in first respa step.
            pint_env%uf = 0.0_dp
            ! calc harmonic forces at new pos
            pint_env%ux = pint_env%ux_t

            ! Apply centroid constraints (SHAKE)
            IF (pint_env%simpar%constraint) THEN
               IF (pint_env%logger%para_env%ionode) THEN
                  DO i = 1, nparticle
                     DO j = 1, 3
                        pos(j, i) = pint_env%ux(1, j + (i - 1)*3)
                        vel(j, i) = pint_env%uv_t(1, j + (i - 1)*3)
                     END DO
                  END DO

                  ! Possibly update the target values
                  CALL shake_update_targets(gci, local_molecules, molecule_set, &
                                            molecule_kind_set, dti, &
                                            f_env%force_env%root_section)
                  CALL shake_control(gci, local_molecules, molecule_set, &
                                     molecule_kind_set, particle_set, &
                                     pos, vel, dti, pint_env%simpar%shake_tol, &
                                     pint_env%simpar%info_constraint, &
                                     pint_env%simpar%lagrange_multipliers, &
                                     pint_env%simpar%dump_lm, cell, &
                                     mp_comm_self, local_particles)
               END IF
               ! Positions and velocities of centroid were constrained by SHAKE
               CALL mp_bcast(pos, &
                             pint_env%logger%para_env%source, &
                             pint_env%logger%para_env%group)
               CALL mp_bcast(vel, &
                             pint_env%logger%para_env%source, &
                             pint_env%logger%para_env%group)
               ! Transform back to normal modes:
               DO i = 1, nparticle
                  DO j = 1, 3
                     pint_env%ux(1, j + (i - 1)*3) = pos(j, i)
                     pint_env%uv_t(1, j + (i - 1)*3) = vel(j, i)
                  END DO
               END DO

            END IF

            CALL pint_calc_uf_h(pint_env=pint_env, e_h=e_h)
            pint_env%uv_t = pint_env%uv_t + dti2*(pint_env%uf_h + pint_env%uf)

            ! For last respa step include integration of physical and helium
            ! forces
            IF (iresp == pint_env%nrespa) THEN
               CALL pint_u2x(pint_env)
               CALL pint_calc_f(pint_env)
               ! perform helium step and add helium forces
               IF (PRESENT(helium_env)) THEN
                  CALL helium_step(helium_env, pint_env)
                  !Update force of solute in pint_env
                  IF (pint_env%logger%para_env%ionode) THEN
                     pint_env%f(:, :) = pint_env%f(:, :) + helium_env(1)%helium%force_avrg(:, :)
                  END IF
                  CALL mp_bcast(pint_env%f, &
                                pint_env%logger%para_env%source, &
                                pint_env%logger%para_env%group)
               END IF

               CALL pint_f2uf(pint_env)
               ! set the centroid forces to 0 if FIX_CENTROID_POS
               IF (pint_env%first_propagated_mode .EQ. 2) THEN
                  pint_env%uf(1, :) = 0.0_dp
               END IF
               !Scale physical forces and integrate velocities with physical
               !forces
               pint_env%uf = pint_env%uf*rn
               pint_env%uv_t = pint_env%uv_t + dti2*pint_env%uf

            END IF

            ! Apply second half of thermostats
            SELECT CASE (pint_env%pimd_thermostat)
            CASE (thermostat_nose)
               DO i = 1, 6
                  tol = 0._dp
                  pint_env%uv_new = pint_env%uv_t/(1.+dti2*pint_env%tv(1, :, :))
                  DO idim = 1, pint_env%ndim
                     DO ib = 1, pint_env%p
                        pint_env%tf(1, ib, idim) = (pint_env%mass_fict(ib, idim)* &
                                                    pint_env%uv_new(ib, idim)**2 - pint_env%kT*pint_env%kTcorr)/ &
                                                   pint_env%Q(ib)
                     END DO
                  END DO

                  !Set thermostat action of constrained DoF to zero:
                  IF (pint_env%simpar%constraint) THEN
                     DO k = 1, pint_env%n_atoms_constraints
                        ia = pint_env%atoms_constraints(k)
                        DO j = 3*(ia - 1) + 1, 3*ia
                           pint_env%tf(:, 1, j) = 0.0_dp
                        END DO
                     END DO
                  END IF

                  DO idim = 1, pint_env%ndim
                     DO ib = 1, pint_env%p
                        DO inos = 1, pint_env%nnos - 1
                           pint_env%tv_new(inos, ib, idim) = &
                              (pint_env%tv_t(inos, ib, idim) + dti2*pint_env%tf(inos, ib, idim))/ &
                              (1._dp + dti2*pint_env%tv(inos + 1, ib, idim))
                           pint_env%tf(inos + 1, ib, idim) = &
                              (pint_env%tv_new(inos, ib, idim)**2 - &
                               pint_env%kT*pint_env%kTcorr/pint_env%Q(ib))
                           tol = MAX(tol, ABS(pint_env%tv(inos, ib, idim) &
                                              - pint_env%tv_new(inos, ib, idim)))
                        END DO
                        !Set thermostat action of constrained DoF to zero:
                        IF (pint_env%simpar%constraint) THEN
                           DO k = 1, pint_env%n_atoms_constraints
                              ia = pint_env%atoms_constraints(k)
                              DO j = 3*(ia - 1) + 1, 3*ia
                                 pint_env%tv_new(:, 1, j) = 0.0_dp
                                 pint_env%tf(:, 1, j) = 0.0_dp
                              END DO
                           END DO
                        END IF

                        pint_env%tv_new(pint_env%nnos, ib, idim) = &
                           pint_env%tv_t(pint_env%nnos, ib, idim) + &
                           dti2*pint_env%tf(pint_env%nnos, ib, idim)
                        tol = MAX(tol, ABS(pint_env%tv(pint_env%nnos, ib, idim) &
                                           - pint_env%tv_new(pint_env%nnos, ib, idim)))
                        tol = MAX(tol, ABS(pint_env%uv(ib, idim) &
                                           - pint_env%uv_new(ib, idim)))
                        !Set thermostat action of constrained DoF to zero:
                        IF (pint_env%simpar%constraint) THEN
                           DO k = 1, pint_env%n_atoms_constraints
                              ia = pint_env%atoms_constraints(k)
                              DO j = 3*(ia - 1) + 1, 3*ia
                                 pint_env%tv_new(:, 1, j) = 0.0_dp
                              END DO
                           END DO
                        END IF

                     END DO
                  END DO

                  pint_env%uv = pint_env%uv_new
                  pint_env%tv = pint_env%tv_new
                  IF (tol <= pint_env%v_tol) EXIT
               END DO

               ! Apply centroid constraints (RATTLE)
               IF (pint_env%simpar%constraint) THEN
                  IF (pint_env%logger%para_env%ionode) THEN
                     ! Reset particle r, due to force calc:
                     DO i = 1, nparticle
                        DO j = 1, 3
                           vel(j, i) = pint_env%uv(1, j + (i - 1)*3)
                           particle_set(i)%r(j) = pint_env%ux(1, j + (i - 1)*3)
                        END DO
                     END DO

                     ! Small time step for all small integrations steps
                     ! Big step for last RESPA
                     IF (iresp == pint_env%nrespa) THEN
                        dt_temp = dti
                     ELSE
                        dt_temp = dti*rn
                     END IF
                     CALL rattle_control(gci, local_molecules, molecule_set, &
                                         molecule_kind_set, particle_set, &
                                         vel, dt_temp, pint_env%simpar%shake_tol, &
                                         pint_env%simpar%info_constraint, &
                                         pint_env%simpar%lagrange_multipliers, &
                                         pint_env%simpar%dump_lm, cell, &
                                         mp_comm_self, local_particles)
                  END IF
                  ! Velocities of centroid were constrained by RATTLE
                  ! Broadcast updated velocities to other nodes
                  CALL mp_bcast(vel, &
                                pint_env%logger%para_env%source, &
                                pint_env%logger%para_env%group)

                  DO i = 1, nparticle
                     DO j = 1, 3
                        pint_env%uv(1, j + (i - 1)*3) = vel(j, i)
                     END DO
                  END DO
               END IF

               DO inos = 1, pint_env%nnos - 1
                  pint_env%tf(inos, :, :) = pint_env%tf(inos, :, :) &
                                            - pint_env%tv(inos, :, :)*pint_env%tv(inos + 1, :, :)
               END DO

            CASE (thermostat_gle)
               CALL pint_gle_step(pint_env)
               pint_env%uv = pint_env%uv_t
            CASE DEFAULT
               pint_env%uv = pint_env%uv_t
            END SELECT
         END DO

      CASE (integrate_exact)
         ! The Liouvillian splitting is as follows:
         ! 1. Thermostat
         ! 2. 0.5*physical integration
         ! 3. Exact harmonic integration ( + apply centroid constraints (SHAKE))
         ! 4. 0.5*physical integration
         ! 5. Thermostat
         ! 6. Apply centroid constraints (RATTLE)

         ! 1. Apply thermostats
         SELECT CASE (pint_env%pimd_thermostat)
         CASE (thermostat_pile)
            CALL pint_pile_step(vold=pint_env%uv, &
                                vnew=pint_env%uv_t, &
                                p=pint_env%p, &
                                ndim=pint_env%ndim, &
                                first_mode=pint_env%first_propagated_mode, &
                                masses=pint_env%mass_fict, &
                                pile_therm=pint_env%pile_therm)
         CASE (thermostat_piglet)
            CALL pint_piglet_step(vold=pint_env%uv, &
                                  vnew=pint_env%uv_t, &
                                  first_mode=pint_env%first_propagated_mode, &
                                  masses=pint_env%mass_fict, &
                                  piglet_therm=pint_env%piglet_therm)
         CASE (thermostat_qtb)
            CALL pint_qtb_step(vold=pint_env%uv, &
                               vnew=pint_env%uv_t, &
                               p=pint_env%p, &
                               ndim=pint_env%ndim, &
                               masses=pint_env%mass_fict, &
                               qtb_therm=pint_env%qtb_therm)
         CASE DEFAULT
            pint_env%uv_t = pint_env%uv
         END SELECT

         ! 2. 1/2*Physical integration
         pint_env%uv_t = pint_env%uv_t + dti2*pint_env%uf

         ! 3. Exact harmonic integration
         IF (pint_env%first_propagated_mode == 1) THEN
            ! The centroid is integrated via standard velocity-verlet
            ! Commented out code is only there to show similarities to
            ! Numeric integrator
            pint_env%ux_t(1, :) = pint_env%ux(1, :) + &
                                  dti*pint_env%uv_t(1, :) !+ &
            !                      dti22*pint_env%uf_h(1, :)
            !pint_env%uv_t(1, :) = pint_env%uv_t(1, :)+ &
            !                      dti2*pint_env%uf_h(1, :)

            ! Apply centroid constraints (SHAKE)
            IF (pint_env%simpar%constraint) THEN
               IF (pint_env%logger%para_env%ionode) THEN
                  ! Transform positions and velocities to Cartesian coordinates:
                  DO i = 1, nparticle
                     DO j = 1, 3
                        pos(j, i) = pint_env%ux_t(1, j + (i - 1)*3)/factor
                        vel(j, i) = pint_env%uv_t(1, j + (i - 1)*3)/factor
                     END DO
                  END DO

                  ! Possibly update the target values
                  CALL shake_update_targets(gci, local_molecules, molecule_set, &
                                            molecule_kind_set, dti, &
                                            f_env%force_env%root_section)
                  CALL shake_control(gci, local_molecules, molecule_set, &
                                     molecule_kind_set, particle_set, &
                                     pos, vel, dti, pint_env%simpar%shake_tol, &
                                     pint_env%simpar%info_constraint, &
                                     pint_env%simpar%lagrange_multipliers, &
                                     pint_env%simpar%dump_lm, cell, &
                                     mp_comm_self, local_particles)
               END IF
               ! Positions and velocities of centroid were constrained by SHAKE
               CALL mp_bcast(pos, &
                             pint_env%logger%para_env%source, &
                             pint_env%logger%para_env%group)
               CALL mp_bcast(vel, &
                             pint_env%logger%para_env%source, &
                             pint_env%logger%para_env%group)
               ! Transform back to normal modes:
               DO i = 1, nparticle
                  DO j = 1, 3
                     pint_env%ux_t(1, j + (i - 1)*3) = pos(j, i)*factor
                     pint_env%uv_t(1, j + (i - 1)*3) = vel(j, i)*factor
                  END DO
               END DO

            END IF

         ELSE
            ! set velocities to zero for fixed centroids
            pint_env%ux_t(1, :) = pint_env%ux(1, :)
            pint_env%uv_t(1, :) = 0.0_dp
         END IF

         DO i = 2, pint_env%p
            pint_env%ux_t(i, :) = pint_env%cosex(i)*pint_env%ux(i, :) &
                                  + pint_env%iwsinex(i)*pint_env%uv_t(i, :)
            pint_env%uv_t(i, :) = pint_env%cosex(i)*pint_env%uv_t(i, :) &
                                  - pint_env%wsinex(i)*pint_env%ux(i, :)
         END DO
         pint_env%ux = pint_env%ux_t

         ! 4. 1/2*Physical integration
         pint_env%uf = 0.0_dp
         CALL pint_u2x(pint_env)
         !TODO apply bead-wise constraints
         CALL pint_calc_f(pint_env)
         ! perform helium step and add helium forces
         IF (PRESENT(helium_env)) THEN
            CALL helium_step(helium_env, pint_env)
            !Update force of solute in pint_env
            IF (pint_env%logger%para_env%ionode) THEN
               pint_env%f(:, :) = pint_env%f(:, :) + helium_env(1)%helium%force_avrg(:, :)
            END IF
            CALL mp_bcast(pint_env%f, &
                          pint_env%logger%para_env%source, &
                          pint_env%logger%para_env%group)
         END IF
         CALL pint_f2uf(pint_env)
         ! set the centroid forces to 0 if FIX_CENTROID_POS
         IF (pint_env%first_propagated_mode .EQ. 2) THEN
            pint_env%uf(1, :) = 0.0_dp
         END IF
         pint_env%uv_t = pint_env%uv_t + dti2*pint_env%uf

         ! 5. Apply thermostats
         SELECT CASE (pint_env%pimd_thermostat)
         CASE (thermostat_pile)
            CALL pint_pile_step(vold=pint_env%uv_t, &
                                vnew=pint_env%uv, &
                                p=pint_env%p, &
                                ndim=pint_env%ndim, &
                                first_mode=pint_env%first_propagated_mode, &
                                masses=pint_env%mass_fict, &
                                pile_therm=pint_env%pile_therm)
         CASE (thermostat_piglet)
            CALL pint_piglet_step(vold=pint_env%uv_t, &
                                  vnew=pint_env%uv, &
                                  first_mode=pint_env%first_propagated_mode, &
                                  masses=pint_env%mass_fict, &
                                  piglet_therm=pint_env%piglet_therm)
         CASE (thermostat_qtb)
            CALL pint_qtb_step(vold=pint_env%uv_t, &
                               vnew=pint_env%uv, &
                               p=pint_env%p, &
                               ndim=pint_env%ndim, &
                               masses=pint_env%mass_fict, &
                               qtb_therm=pint_env%qtb_therm)
         CASE DEFAULT
            pint_env%uv = pint_env%uv_t
         END SELECT

         ! 6. Apply centroid constraints (RATTLE)
         IF (pint_env%simpar%constraint) THEN
            IF (pint_env%logger%para_env%ionode) THEN
               ! Transform positions and velocities to Cartesian coordinates:
               ! Reset particle r, due to force calc:
               DO i = 1, nparticle
                  DO j = 1, 3
                     vel(j, i) = pint_env%uv(1, j + (i - 1)*3)/factor
                     particle_set(i)%r(j) = pint_env%ux(1, j + (i - 1)*3)/factor
                  END DO
               END DO
               CALL rattle_control(gci, local_molecules, &
                                   molecule_set, molecule_kind_set, &
                                   particle_set, vel, dti, &
                                   pint_env%simpar%shake_tol, &
                                   pint_env%simpar%info_constraint, &
                                   pint_env%simpar%lagrange_multipliers, &
                                   pint_env%simpar%dump_lm, cell, &
                                   mp_comm_self, local_particles)
            END IF
            ! Velocities of centroid were constrained by RATTLE
            ! Broadcast updated velocities to other nodes
            CALL mp_bcast(vel, &
                          pint_env%logger%para_env%source, &
                          pint_env%logger%para_env%group)

            ! Transform back to normal modes:
            ! Multiply with SQRT(n_beads) due to normal mode transformation
            DO i = 1, nparticle
               DO j = 1, 3
                  pint_env%uv(1, j + (i - 1)*3) = vel(j, i)*factor
               END DO
            END DO

            !TODO apply bead-wise constraints
         END IF

      END SELECT

      IF (pint_env%simpar%constraint) DEALLOCATE (pos, vel)

      ! calculate the energy components
      CALL pint_calc_energy(pint_env)
      CALL pint_calc_total_action(pint_env)

      ! check that the number of PINT steps matches
      ! the number of force evaluations done so far
!TODO make this check valid if we start from ITERATION != 0
!     CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id,&
!          f_env=f_env,new_error=new_error)
!     NULLIFY(logger)
!     logger => cp_error_get_logger(new_error)
!     IF(logger%iter_info%iteration(2)/=pint_env%iter+1)&
!        CPABORT("md & force_eval lost sychro")
!     CALL f_env_rm_defaults(f_env,new_error,ierr)

      time_stop = m_walltime()
      pint_env%time_per_step = time_stop - time_start
      CALL pint_write_step_info(pint_env)
      CALL timestop(handle)

      RETURN
   END SUBROUTINE pint_step

! ***************************************************************************
!> \brief  Calculate the energy components (private wrapper function)
!> \param pint_env ...
!> \date   2011-01-07
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE pint_calc_energy(pint_env)

      TYPE(pint_env_type), POINTER                       :: pint_env

      REAL(KIND=dp)                                      :: e_h

      CALL pint_calc_e_kin_beads_u(pint_env)
      CALL pint_calc_e_vir(pint_env)

      CALL pint_calc_uf_h(pint_env, e_h=e_h)
      pint_env%e_pot_h = e_h

      SELECT CASE (pint_env%pimd_thermostat)
      CASE (thermostat_nose)
         CALL pint_calc_nh_energy(pint_env)
      CASE (thermostat_gle)
         CALL pint_calc_gle_energy(pint_env)
      CASE (thermostat_pile)
         CALL pint_calc_pile_energy(pint_env)
      CASE (thermostat_qtb)
         CALL pint_calc_qtb_energy(pint_env)
      CASE (thermostat_piglet)
         CALL pint_calc_piglet_energy(pint_env)
      END SELECT

      pint_env%energy(e_kin_thermo_id) = &
         (0.5_dp*REAL(pint_env%p, dp)*REAL(pint_env%ndim, dp)*pint_env%kT - &
          pint_env%e_pot_h)*pint_env%propagator%temp_sim2phys

      pint_env%energy(e_potential_id) = SUM(pint_env%e_pot_bead)

      pint_env%energy(e_conserved_id) = &
         pint_env%energy(e_potential_id)*pint_env%propagator%physpotscale + &
         pint_env%e_pot_h + &
         pint_env%e_kin_beads + &
         pint_env%e_pot_t + &
         pint_env%e_kin_t + &
         pint_env%e_gle + pint_env%e_pile + pint_env%e_piglet + pint_env%e_qtb

      pint_env%energy(e_potential_id) = &
         pint_env%energy(e_potential_id)/REAL(pint_env%p, dp)

      RETURN
   END SUBROUTINE pint_calc_energy

! ***************************************************************************
!> \brief  Calculate the harmonic force in the u basis
!> \param  pint_env the path integral environment in which the harmonic
!>         forces should be calculated
!> \param e_h ...
!> \par    History
!>           Added normal mode transformation [hforbert]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE pint_calc_uf_h(pint_env, e_h)
      TYPE(pint_env_type), POINTER                       :: pint_env
      REAL(KIND=dp), INTENT(OUT)                         :: e_h

      IF (pint_env%transform == transformation_stage) THEN
         CALL staging_calc_uf_h(pint_env%staging_env, &
                                pint_env%mass_beads, &
                                pint_env%ux, &
                                pint_env%uf_h, &
                                pint_env%e_pot_h)
      ELSE
         CALL normalmode_calc_uf_h(pint_env%normalmode_env, &
                                   pint_env%mass_beads, &
                                   pint_env%ux, &
                                   pint_env%uf_h, &
                                   pint_env%e_pot_h)
      END IF
      e_h = pint_env%e_pot_h
      pint_env%uf_h = pint_env%uf_h/pint_env%mass_fict
      RETURN
   END SUBROUTINE pint_calc_uf_h

! ***************************************************************************
!> \brief calculates the force (and energy) in each bead, returns the sum
!>      of the potential energy
!> \param pint_env path integral environment on which you want to calculate
!>        the forces
!> \param x positions at which you want to evaluate the forces
!> \param f the forces
!> \param e potential energy on each bead
!> \par    History
!>           2009-06-15 moved helium calls out from here [lwalewski]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE pint_calc_f(pint_env, x, f, e)
      TYPE(pint_env_type), POINTER                       :: pint_env
      REAL(kind=dp), DIMENSION(:, :), INTENT(in), &
         OPTIONAL, TARGET                                :: x
      REAL(kind=dp), DIMENSION(:, :), INTENT(out), &
         OPTIONAL, TARGET                                :: f
      REAL(kind=dp), DIMENSION(:), INTENT(out), &
         OPTIONAL, TARGET                                :: e

      INTEGER                                            :: ib, idim
      REAL(kind=dp), DIMENSION(:), POINTER               :: my_e
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: my_f, my_x

      CPASSERT(ASSOCIATED(pint_env))
      CPASSERT(pint_env%ref_count > 0)
      my_x => pint_env%x
      IF (PRESENT(x)) my_x => x
      my_f => pint_env%f
      IF (PRESENT(f)) my_f => f
      my_e => pint_env%e_pot_bead
      IF (PRESENT(e)) my_e => e
      DO idim = 1, pint_env%ndim
         DO ib = 1, pint_env%p
            pint_env%replicas%r(idim, ib) = my_x(ib, idim)
         END DO
      END DO
      CALL rep_env_calc_e_f(pint_env%replicas, calc_f=.TRUE.)
      DO idim = 1, pint_env%ndim
         DO ib = 1, pint_env%p
            !ljw: is that fine ? - idim <-> ib
            my_f(ib, idim) = pint_env%replicas%f(idim, ib)
         END DO
      END DO
      my_e = pint_env%replicas%f(SIZE(pint_env%replicas%f, 1), :)

      RETURN
   END SUBROUTINE pint_calc_f

! ***************************************************************************
!> \brief  Calculate the kinetic energy of the beads (in the u variables)
!> \param pint_env ...
!> \param uv ...
!> \param e_k ...
!> \par    History
!>         Bug fix to give my_uv a default location if not given in call [hforbert]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE pint_calc_e_kin_beads_u(pint_env, uv, e_k)
      TYPE(pint_env_type), POINTER                       :: pint_env
      REAL(kind=dp), DIMENSION(:, :), INTENT(in), &
         OPTIONAL, TARGET                                :: uv
      REAL(kind=dp), INTENT(out), OPTIONAL               :: e_k

      INTEGER                                            :: ib, idim
      REAL(kind=dp)                                      :: res
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: my_uv

      CPASSERT(ASSOCIATED(pint_env))
      CPASSERT(pint_env%ref_count > 0)
      res = -1.0_dp
      my_uv => pint_env%uv
      IF (PRESENT(uv)) my_uv => uv
      res = 0._dp
      DO idim = 1, pint_env%ndim
         DO ib = 1, pint_env%p
            res = res + pint_env%mass_fict(ib, idim)*my_uv(ib, idim)**2
         END DO
      END DO
      res = res*0.5
      IF (.NOT. PRESENT(uv)) pint_env%e_kin_beads = res
      IF (PRESENT(e_k)) e_k = res
      RETURN
   END SUBROUTINE pint_calc_e_kin_beads_u

! ***************************************************************************
!> \brief  Calculate the virial estimator of the real (quantum) kinetic energy
!> \param pint_env ...
!> \param e_vir ...
!> \author hforbert
!> \note   This subroutine modifies pint_env%energy(e_kin_virial_id) global
!>         variable [lwalewski]
! **************************************************************************************************
   SUBROUTINE pint_calc_e_vir(pint_env, e_vir)
      TYPE(pint_env_type), POINTER                       :: pint_env
      REAL(kind=dp), INTENT(out), OPTIONAL               :: e_vir

      INTEGER                                            :: ib, idim
      REAL(kind=dp)                                      :: res, xcentroid

      CPASSERT(ASSOCIATED(pint_env))
      CPASSERT(pint_env%ref_count > 0)
      res = -1.0_dp
      res = 0._dp
      DO idim = 1, pint_env%ndim
         ! calculate the centroid
         xcentroid = 0._dp
         DO ib = 1, pint_env%p
            xcentroid = xcentroid + pint_env%x(ib, idim)
         END DO
         xcentroid = xcentroid/REAL(pint_env%p, dp)
         DO ib = 1, pint_env%p
            res = res + (pint_env%x(ib, idim) - xcentroid)*pint_env%f(ib, idim)
         END DO
      END DO
      res = 0.5_dp*(REAL(pint_env%ndim, dp)* &
                    (pint_env%kT*pint_env%propagator%temp_sim2phys) - res/REAL(pint_env%p, dp))
      pint_env%energy(e_kin_virial_id) = res
      IF (PRESENT(e_vir)) e_vir = res
      RETURN
   END SUBROUTINE pint_calc_e_vir

! ***************************************************************************
!> \brief calculates the energy (potential and kinetic) of the Nose-Hoover
!>      chain thermostats
!> \param pint_env the path integral environment
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE pint_calc_nh_energy(pint_env)
      TYPE(pint_env_type), POINTER                       :: pint_env

      INTEGER                                            :: ib, idim, inos
      REAL(kind=dp)                                      :: ekin, epot

      CPASSERT(ASSOCIATED(pint_env))
      CPASSERT(pint_env%ref_count > 0)
      ekin = 0._dp
      DO idim = 1, pint_env%ndim
         DO ib = 1, pint_env%p
            DO inos = 1, pint_env%nnos
               ekin = ekin + pint_env%Q(ib)*pint_env%tv(inos, ib, idim)**2
            END DO
         END DO
      END DO
      pint_env%e_kin_t = 0.5_dp*ekin
      epot = 0._dp
      DO idim = 1, pint_env%ndim
         DO ib = 1, pint_env%p
            DO inos = 1, pint_env%nnos
               epot = epot + pint_env%tx(inos, ib, idim)
            END DO
         END DO
      END DO
      pint_env%e_pot_t = pint_env%kT*epot
      RETURN
   END SUBROUTINE pint_calc_nh_energy

! ***************************************************************************
!> \brief calculates the total link action of the PI system (excluding helium)
!> \param pint_env the path integral environment
!> \return ...
!> \author Felix Uhl
! **************************************************************************************************
   FUNCTION pint_calc_total_link_action(pint_env) RESULT(link_action)
      TYPE(pint_env_type), POINTER                       :: pint_env
      REAL(KIND=dp)                                      :: link_action

      INTEGER                                            :: iatom, ibead, idim, indx
      REAL(KIND=dp)                                      :: hb2m, tau, tmp_link_action
      REAL(KIND=dp), DIMENSION(3)                        :: r

      !tau = 1/(k_B T p)
      tau = pint_env%beta/REAL(pint_env%p, dp)

      CPASSERT(ASSOCIATED(pint_env))
      CPASSERT(pint_env%ref_count > 0)
      link_action = 0.0_dp
      DO iatom = 1, pint_env%ndim/3
         ! hbar / (2.0*m)
         hb2m = 1.0_dp/pint_env%mass((iatom - 1)*3 + 1)
         tmp_link_action = 0.0_dp
         DO ibead = 1, pint_env%p - 1
            DO idim = 1, 3
               indx = (iatom - 1)*3 + idim
               r(idim) = pint_env%x(ibead, indx) - pint_env%x(ibead + 1, indx)
            END DO
            tmp_link_action = tmp_link_action + (r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
         END DO
         DO idim = 1, 3
            indx = (iatom - 1)*3 + idim
            r(idim) = pint_env%x(pint_env%p, indx) - pint_env%x(1, indx)
         END DO
         tmp_link_action = tmp_link_action + (r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
         link_action = link_action + tmp_link_action/hb2m
      END DO

      link_action = link_action/(2.0_dp*tau)

      RETURN

   END FUNCTION pint_calc_total_link_action

! ***************************************************************************
!> \brief calculates the potential action of the PI system (excluding helium)
!> \param pint_env the path integral environment
!> \return ...
!> \author Felix Uhl
! **************************************************************************************************
   FUNCTION pint_calc_total_pot_action(pint_env) RESULT(pot_action)
      TYPE(pint_env_type), POINTER                       :: pint_env
      REAL(KIND=dp)                                      :: pot_action

      REAL(KIND=dp)                                      :: tau

      CPASSERT(ASSOCIATED(pint_env))
      CPASSERT(pint_env%ref_count > 0)

      tau = pint_env%beta/REAL(pint_env%p, dp)
      pot_action = tau*SUM(pint_env%e_pot_bead)

      RETURN

   END FUNCTION pint_calc_total_pot_action

! ***************************************************************************
!> \brief calculates the total action of the PI system (excluding helium)
!> \param pint_env the path integral environment
!> \author Felix Uhl
! **************************************************************************************************
   SUBROUTINE pint_calc_total_action(pint_env)
      TYPE(pint_env_type), POINTER                       :: pint_env

      CPASSERT(ASSOCIATED(pint_env))
      CPASSERT(pint_env%ref_count > 0)

      pint_env%pot_action = pint_calc_total_pot_action(pint_env)
      pint_env%link_action = pint_calc_total_link_action(pint_env)

      RETURN
   END SUBROUTINE pint_calc_total_action

END MODULE pint_methods
