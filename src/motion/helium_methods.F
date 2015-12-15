!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Methods dealing with helium_solvent_type
!> \author Lukasz Walewski
!> \date   2009-06-10
! *****************************************************************************
MODULE helium_methods

  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE bibliography,                    ONLY: Walewski2014,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                             cp_logger_type
  USE cp_output_handling,              ONLY: cp_printkey_is_on
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE f77_interface,                   ONLY: f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get
  USE helium_common,                   ONLY: helium_com,&
                                             helium_pbc
  USE helium_io,                       ONLY: helium_read_xyz,&
                                             helium_write_line,&
                                             helium_write_setup
  USE helium_sampling,                 ONLY: helium_sample
  USE helium_types,                    ONLY: &
       he_mass, helium_destroy_int_arr_ptr, helium_solvent_type, hid_carbon, &
       hid_chlorine, hid_hydrogen, hid_num, hid_oxygen, rho_atom_number, &
       rho_moment_of_inertia, rho_num, rho_projected_area, rho_winding_cycle, &
       rho_winding_number
  USE input_constants,                 ONLY: helium_cell_shape_cube,&
                                             helium_cell_shape_octahedron,&
                                             helium_solute_intpot_none
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp,&
                                             max_line_length
  USE mathconstants,                   ONLY: twopi
  USE message_passing,                 ONLY: mp_bcast
  USE parallel_rng_types,              ONLY: GAUSSIAN,&
                                             UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type,&
                                             set_rng_stream
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: a_mass,&
                                             angstrom,&
                                             boltzmann,&
                                             h_bar,&
                                             kelvin,&
                                             massunit
  USE pint_public,                     ONLY: pint_com_pos
  USE pint_types,                      ONLY: pint_env_type
  USE splines_methods,                 ONLY: init_spline,&
                                             init_splinexy
  USE splines_types,                   ONLY: spline_data_create,&
                                             spline_data_release,&
                                             spline_data_retain
#include "../base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_methods'
  INTEGER, SAVE, PRIVATE :: last_helium_id=0

  PUBLIC :: helium_create
  PUBLIC :: helium_init
  PUBLIC :: helium_release

  CONTAINS

! ***************************************************************************
!> \brief  Data-structure that holds all needed information about
!>         (superfluid) helium solvent
!> \param helium ...
!> \param input ...
!> \param solute ...
!> \author hforbert
! *****************************************************************************
  SUBROUTINE helium_create( helium, input, solute)
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(section_vals_type), POINTER         :: input
    TYPE(pint_env_type), OPTIONAL, POINTER   :: solute

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: msg_str, potential_file_name
    CHARACTER(len=default_string_length)     :: stmp
    INTEGER                                  :: handle, i, input_unit, isize, &
                                                itmp, j, nlines, ntab
    LOGICAL                                  :: expl_cell, expl_dens, &
                                                expl_nats, explicit, ltmp
    REAL(KIND=dp)                            :: cgeof, dx, mHe, rtmp, T, &
                                                tcheck, x1
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pot_transfer
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: helium_section

    CALL timeset(routineN,handle)

    CALL cite_reference(Walewski2014)
    NULLIFY(helium_section)
    helium_section => section_vals_get_subs_vals(input, &
                                           "MOTION%PINT%HELIUM")
    CALL section_vals_get(helium_section,explicit=explicit)
    ALLOCATE(helium)
    NULLIFY ( helium%input, &
              helium%ptable,  helium%permutation, &
              helium%iperm, &
              helium%itmp_atoms_1d, &
              helium%ltmp_atoms_1d, &
              helium%itmp_atoms_np_1d, &
              helium%pos, helium%work, &
              helium%force_avrg, &
              helium%force_inst, &
              helium%rtmp_3_np_1d, &
              helium%rtmp_p_ndim_1d, &
              helium%rtmp_p_ndim_np_1d, &
              helium%rtmp_3_atoms_beads_1d, &
              helium%rtmp_3_atoms_beads_np_1d, &
              helium%rtmp_p_ndim_2d, &
              helium%ltmp_3_atoms_beads_3d,&
              helium%tmatrix, helium%pmatrix,     &
              helium%nmatrix, helium%ipmatrix,    &
              helium%uij,     helium%eij, &
              helium%rdf_inst,&
              helium%plength_avrg, &
              helium%plength_inst, &
              helium%atom_plength, &
              helium%ename, &
              helium%eid&
             )

    helium%ref_count=1
    last_helium_id=last_helium_id+1
    helium%id_nr=last_helium_id
    helium%input => input
    helium%accepts = 0
    helium%relrot = 0

    NULLIFY(logger)
    logger => cp_get_default_logger()

    ! check if solute is present in our simulation
    helium%solute_present = .FALSE.
    helium%solute_atoms = 0
    helium%solute_beads = 0
    CALL section_vals_val_get( &
      helium_section, &
      "HELIUM_ONLY", &
      l_val=ltmp)
    IF (.NOT. ltmp) THEN
      IF (PRESENT(solute)) THEN
        IF (ASSOCIATED(solute)) THEN
          helium%solute_present = .TRUE.
          helium%solute_atoms = solute%ndim / 3
          helium%solute_beads = solute%p
        END IF
      END IF
    END IF

    ! get number of environments in the restart file (if present)
    CALL section_vals_val_get(helium_section,"NUM_ENV",&
         explicit=explicit)
    IF ( explicit ) THEN
      CALL section_vals_val_get(helium_section,"NUM_ENV",&
           i_val=itmp)
      helium%num_env_restart = itmp
    ELSE
      helium%num_env_restart = -1
    END IF

    ! set current number of environments
    helium%num_env = logger%para_env%num_pe
    CALL section_vals_val_set(helium%input,&
         "MOTION%PINT%HELIUM%NUM_ENV",&
         i_val=helium%num_env)

    ! if the restart file contains more environments than the runtime
    ! the additional environments (possibly with valuable permutation state)
    ! will be lost, we won't proceed in this case by default, the user might
    ! enforce this explicitly though
    IF ( helium%num_env_restart .GT. helium%num_env ) THEN
      stmp = ""
      WRITE(stmp,*) helium%num_env
      msg_str = "Number of He environments in the runtime (" // &
      TRIM(ADJUSTL(stmp)) // ") smaller than in the restart ("
      stmp = ""
      WRITE(stmp,*) helium%num_env_restart
      msg_str = TRIM(ADJUSTL(msg_str)) // TRIM(ADJUSTL(stmp)) // ")!"
      CALL helium_write_line(msg_str)
      CALL section_vals_val_get(helium_section,"DROP_UNUSED_ENVS",&
           l_val=ltmp)
      IF (ltmp) THEN
        msg_str = "DROP_UNUSED_ENVS set - proceeding anyways."
        CPWARN(msg_str)
      ELSE
        CPABORT(msg_str)
      END IF
    END IF

    CALL section_vals_val_get(helium_section,"NBEADS",&
           i_val=helium%beads)
    CALL section_vals_val_get(helium_section,"INOROT",&
           i_val=helium%iter_norot)
    CALL section_vals_val_get(helium_section,"IROT",&
           i_val=helium%iter_rot)

    ! get number of steps and current step number from PINT
    CALL section_vals_val_get(input,"MOTION%PINT%ITERATION",&
         i_val=itmp)
    helium%first_step = itmp
    CALL section_vals_val_get(input,"MOTION%PINT%MAX_STEP",&
         explicit=explicit)
    IF ( explicit ) THEN
      CALL section_vals_val_get(input,"MOTION%PINT%MAX_STEP",&
           i_val=itmp)
      helium%last_step = itmp
      helium%num_steps = helium%last_step - helium%first_step
    ELSE
      CALL section_vals_val_get(input,"MOTION%PINT%NUM_STEPS",&
           i_val=itmp)
      helium%num_steps = itmp
      helium%last_step = helium%first_step + helium%num_steps
    END IF

    ! boundary conditions
    CALL section_vals_val_get(helium_section,"PERIODIC",&
         l_val=helium%periodic)
    CALL section_vals_val_get(helium_section,"CELL_SHAPE",&
         i_val=helium%cell_shape)

    CALL section_vals_val_get(helium_section,"DROPLET_RADIUS",&
         r_val=helium%droplet_radius)


    ! Set density Rho, number of atoms N and volume V ( Rho = N / V ).
    ! Allow only 2 out of 3 values to be defined at the same time, calculate
    ! the third.
    ! Note, that DENSITY and NATOMS keywords have default values, while
    ! CELL_SIZE does not. Thus if CELL_SIZE is given explicitly then one and
    ! only one of the two remaining options must be give explicitly as well.
    ! If CELL_SIZE is not given explicitly then all four combinations of the
    ! two other options are valid.
    CALL section_vals_val_get(helium_section,"DENSITY",&
         explicit=expl_dens, r_val=helium%density)
    CALL section_vals_val_get(helium_section,"NATOMS",&
         explicit=expl_nats, i_val=helium%atoms)
    CALL section_vals_val_get(helium_section,"CELL_SIZE",&
         explicit=expl_cell)
    cgeof = 1.0_dp
    IF (helium%periodic) THEN
      IF ( helium%cell_shape .EQ. helium_cell_shape_octahedron ) cgeof = 2.0_dp
    END IF
    rtmp = ( cgeof * helium%atoms / helium%density )**(1.0_dp/3.0_dp)
    IF ( .NOT. expl_cell ) THEN
        helium%cell_size = rtmp
    ELSE
        CALL section_vals_val_get(helium_section,"CELL_SIZE",&
             r_val=helium%cell_size)
        ! only more work if not all three values are consistent:
        IF ( ABS(helium%cell_size-rtmp) .GT. 100.0_dp*EPSILON(0.0_dp)* &
                                       (ABS(helium%cell_size)+rtmp) ) THEN
           IF ( expl_dens .AND. expl_nats ) THEN
              msg_str = "DENSITY, NATOMS and CELL_SIZE options "//&
                        "contradict each other"
              CPWARN(msg_str)
           END IF
           !ok we have enough freedom to resolve the conflict:
           IF ( .NOT. expl_dens ) THEN
              helium%density = cgeof*helium%atoms / helium%cell_size**3.0_dp
              IF ( .NOT. expl_nats ) THEN
                 msg_str = "CELL_SIZE defined but neither "//&
                           "NATOMS nor DENSITY given, using default NATOMS."
                 CPWARN(msg_str)
              END IF
           ELSE ! ( expl_dens .AND. .NOT. expl_nats )
              ! calculate the nearest number of atoms for given conditions
              helium%atoms = NINT(helium%density * &
                             helium%cell_size**3.0_dp / cgeof)
              ! adjust cell size to maintain correct density
              ! (should be a small correction)
              rtmp = ( cgeof * helium%atoms / helium%density &
                     )**(1.0_dp/3.0_dp)
              IF ( ABS(helium%cell_size-rtmp) .GT. 100.0_dp*EPSILON(0.0_dp)&
                   * ( ABS(helium%cell_size)+rtmp ) ) THEN
                 msg_str = "Adjusting actual cell size "//&
                           "to maintain correct density."
                 CPWARN(msg_str)
                 helium%cell_size = rtmp
              END IF
           END IF
        END IF
    END IF
    helium%cell_size_inv = 1.0_dp / helium%cell_size
    ! From now on helium%density, helium%atoms and helium%cell_size are
    ! correctly defined.

    ! set the M matrix for winding number calculations
    SELECT CASE (helium%cell_shape)

      CASE (helium_cell_shape_octahedron)
        helium%cell_m(1,1) = helium%cell_size;
        helium%cell_m(2,1) = 0.0_dp;
        helium%cell_m(3,1) = 0.0_dp;
        helium%cell_m(1,2) = 0.0_dp;
        helium%cell_m(2,2) = helium%cell_size;
        helium%cell_m(3,2) = 0.0_dp;
        helium%cell_m(1,3) = helium%cell_size / 2.0_dp;
        helium%cell_m(2,3) = helium%cell_size / 2.0_dp;
        helium%cell_m(3,3) = helium%cell_size / 2.0_dp;

        helium%cell_m_inv(1,1) =  1.0_dp / helium%cell_size;
        helium%cell_m_inv(2,1) =  0.0_dp;
        helium%cell_m_inv(3,1) =  0.0_dp;
        helium%cell_m_inv(1,2) =  0.0_dp;
        helium%cell_m_inv(2,2) =  1.0_dp / helium%cell_size;
        helium%cell_m_inv(3,2) =  0.0_dp;
        helium%cell_m_inv(1,3) = -1.0_dp / helium%cell_size;
        helium%cell_m_inv(2,3) = -1.0_dp / helium%cell_size;
        helium%cell_m_inv(3,3) =  2.0_dp / helium%cell_size;

      CASE (helium_cell_shape_cube)

        helium%cell_m(1,1) = helium%cell_size;
        helium%cell_m(2,1) = 0.0_dp;
        helium%cell_m(3,1) = 0.0_dp;
        helium%cell_m(1,2) = 0.0_dp;
        helium%cell_m(2,2) = helium%cell_size;
        helium%cell_m(3,2) = 0.0_dp;
        helium%cell_m(1,3) = 0.0_dp;
        helium%cell_m(2,3) = 0.0_dp;
        helium%cell_m(3,3) = helium%cell_size;

        helium%cell_m_inv(1,1) =  1.0_dp / helium%cell_size;
        helium%cell_m_inv(2,1) =  0.0_dp;
        helium%cell_m_inv(3,1) =  0.0_dp;
        helium%cell_m_inv(1,2) =  0.0_dp;
        helium%cell_m_inv(2,2) =  1.0_dp / helium%cell_size;
        helium%cell_m_inv(3,2) =  0.0_dp;
        helium%cell_m_inv(1,3) =  0.0_dp;
        helium%cell_m_inv(2,3) =  0.0_dp;
        helium%cell_m_inv(3,3) =  1.0_dp / helium%cell_size;

      CASE DEFAULT
        helium%cell_m(:,:) = 0.0_dp
        helium%cell_m_inv(:,:) = 0.0_dp

    END SELECT

    ! check value of maxcycle
    CALL section_vals_val_get(helium_section,"MAX_PERM_CYCLE",&
           i_val=helium%maxcycle)
    i = helium%maxcycle
    i = helium%atoms - helium%maxcycle

    ! set m-distribution parameters
    CALL section_vals_val_get(helium_section,"M-SAMPLING%DISTRIBUTION-TYPE",&
         i_val=i)
    helium%m_dist_type = i
    CALL section_vals_val_get(helium_section,"M-SAMPLING%M-VALUE",&
         i_val=i)
    helium%m_value = i
    CALL section_vals_val_get(helium_section,"M-SAMPLING%M-RATIO",&
         r_val=rtmp)
    helium%m_ratio = rtmp

    CALL section_vals_val_get(helium_section,"BISECTION",&
           i_val=helium%bisection)
    ! precheck bisection value (not all invalids are filtered out here yet)
    i = helium%bisection
    i = helium%beads - helium%bisection
    !
    itmp = helium%bisection
    rtmp = 2.0_dp**(ANINT(LOG(REAL(itmp,dp))/LOG(2.0_dp)))
    tcheck=ABS(REAL(itmp,KIND=dp)-rtmp)
    IF (tcheck>100.0_dp*EPSILON(0.0_dp)) THEN
       msg_str = "BISECTION should be integer power of 2."
       CPABORT(msg_str)
    END IF
    helium%bisctlog2 = NINT(LOG(REAL(itmp,dp))/LOG(2.0_dp))

    CALL section_vals_val_get(helium_section,"SAMPLING_METHOD", &
    i_val=helium%sampling_method)

    ! hard coded He4 directly (mass so i get my original hb2m value)
    ! He4 mass defined as a constant in helium_types.F now [lwalewski]
    helium%hb2m = 1.0_dp/(he_mass*massunit)

    helium%pweight = 0.0_dp

    IF (logger%para_env%ionode) THEN
      CALL section_vals_val_get(helium_section,"POTENTIAL_FILE_NAME",&
          c_val=potential_file_name)
      CALL open_file(file_name=TRIM(potential_file_name), &
          file_action="READ", file_status="OLD",unit_number=input_unit)
      READ (input_unit,*) nlines, helium%pdx, helium%tau,&
                          x1,dx
      helium%tau = kelvin/helium%tau
      x1 = x1/angstrom
      dx = dx/angstrom
    END IF
    CALL mp_bcast(nlines,logger%para_env%source,&
                  logger%para_env%group)
    CALL mp_bcast(helium%pdx,logger%para_env%source,&
                  logger%para_env%group)
    CALL mp_bcast(helium%tau,logger%para_env%source,&
                  logger%para_env%group)
    CALL mp_bcast(x1,logger%para_env%source,&
                  logger%para_env%group)
    CALL mp_bcast(dx,logger%para_env%source,&
                  logger%para_env%group)

    ! boltzmann : Boltzmann constant [J/K]
    ! h_bar     : Planck constant [J*s]
    ! J = kg*m^2/s^2
    ! 4He mass in [kg]
    mHe = he_mass * a_mass
    ! physical temperature [K]
    T = kelvin / helium%tau / helium%beads
    ! prefactors for calculating superfluid fractions [Angstrom^-2]
    helium%wpref = (((1e-20/h_bar)*mHe)/h_bar)*boltzmann*T
    helium%apref = (((4e-20/h_bar)*mHe)/h_bar)*boltzmann*T

    isize = helium%pdx+1
    ALLOCATE(helium%uij(isize,isize))
    ALLOCATE(helium%eij(isize,isize))
    DO i = 1, isize
       DO j = 1, i
          CALL spline_data_create(helium%uij(i,j)%spline_data)
          CALL init_splinexy(helium%uij(i,j)%spline_data,nlines)
          helium%uij(i,j)%spline_data%x1 = x1
          CALL spline_data_create(helium%eij(i,j)%spline_data)
          CALL init_splinexy(helium%eij(i,j)%spline_data,nlines)
          helium%eij(i,j)%spline_data%x1 = x1
       END DO
    END DO
    DO i = 1, isize-1
       DO j = i+1, isize
          helium%uij(i,j) = helium%uij(j,i)
          CALL spline_data_retain(helium%uij(i,j)%spline_data)
          helium%eij(i,j) = helium%eij(j,i)
          CALL spline_data_retain(helium%eij(i,j)%spline_data)
       END DO
    END DO

    isize = (helium%pdx+1)*(helium%pdx+2)
    ALLOCATE(pot_transfer(nlines,isize))
    IF (logger%para_env%ionode) THEN
       DO i = 1, nlines
          READ (input_unit,*) pot_transfer(i,:)
       END DO
       CALL close_file(unit_number=input_unit)
    END IF
    CALL mp_bcast(pot_transfer,logger%para_env%source,&
                  logger%para_env%group)
    isize = helium%pdx+1
    ntab = 1
    DO i = 1, isize
       DO j = 1, i
          helium%uij(i,j)%spline_data%y(:)=pot_transfer(:,ntab)* &
                                                   angstrom**(2*i-2)
          CALL init_spline(helium%uij(i,j)%spline_data,dx=dx)
          ntab = ntab + 1
       END DO
    END DO
    DO i = 1, isize
       DO j = 1, i
          helium%eij(i,j)%spline_data%y(:)=pot_transfer(:,ntab)* &
                                                   angstrom**(2*i-2)/kelvin
          CALL init_spline(helium%eij(i,j)%spline_data,dx=dx)
          ntab = ntab + 1
       END DO
    END DO
    DEALLOCATE(pot_transfer)

    ! ALLOCATE helium-related arrays
    i = helium%atoms
    j = helium%beads
    ALLOCATE(helium%pos(3,i,j+1))
    helium%pos = 0.0_dp
    ALLOCATE(helium%work(3,i,j+1))
    ALLOCATE(helium%ptable(helium%maxcycle+1))
    ALLOCATE(helium%permutation(i))
    ALLOCATE(helium%iperm(i))
    ALLOCATE(helium%tmatrix(i,i))
    ALLOCATE(helium%nmatrix(i,2*i))
    ALLOCATE(helium%pmatrix(i,i))
    ALLOCATE(helium%ipmatrix(i,i))
    itmp = helium%bisctlog2 + 2
    ALLOCATE(helium%num_accepted(itmp,helium%maxcycle))
    ALLOCATE(helium%plength_avrg(helium%atoms))
    ALLOCATE(helium%plength_inst(helium%atoms))
    ALLOCATE(helium%atom_plength(helium%atoms))

    ! check whether rdfs should be calculated and printed
    helium%rdf_present = helium_property_active(helium,"RDF")
    IF (helium%rdf_present) THEN
      ! allocate & initialize rdf related data structures
      CALL helium_rdf_init(helium)
    END IF

    ! check whether densities should be calculated and printed
    helium%rho_present = helium_property_active(helium,"RHO")
    IF (helium%rho_present) THEN
      ! allocate & initialize density related data structures
      NULLIFY(helium%rho_property)
      CALL helium_rho_init(helium)
    END IF

    ! restore averages calculated in previous runs
    CALL helium_averages_restore(helium)

    ! RNG state create & init
    CALL helium_rng_init( helium )

    ! fill in the solute-related data structures
    helium%e_corr = 0.0_dp
    IF (helium%solute_present) THEN
      helium%bead_ratio = helium%beads / helium%solute_beads

      ! check if bead numbers are commensurate:
      i = helium%bead_ratio*helium%solute_beads - helium%beads
!TODO Adjust helium bead number if not comm. and if coords not given expl.

      ! check if tau, temperature and bead number are consistent:
      tcheck=ABS( (helium%tau*helium%beads-solute%beta) / solute%beta )
      IF (tcheck>1.0e-14_dp) THEN
         msg_str = "Tau, temperature and bead number are inconsistent."
         CPABORT(msg_str)
      END IF

      CALL helium_set_solute_indices(helium,solute)
      CALL helium_set_solute_cell(helium,solute)

      ! set the interaction potential type
      CALL section_vals_val_get(helium_section,"SOLUTE_INTERACTION",&
           i_val=helium%solute_interaction)
      IF (helium%solute_interaction .EQ. helium_solute_intpot_none) THEN
        WRITE(msg_str,'(A,I0,A)') &
         "Solute found but no helium-solute interaction selected "//&
         "(see SOLUTE_INTERACTION keyword)"
        CPABORT(msg_str)
      END IF

      ! ALLOCATE solute-related arrays
      ALLOCATE(helium%force_avrg(helium%solute_beads,&
        helium%solute_atoms*3))
      ALLOCATE(helium%force_inst(helium%solute_beads,&
        helium%solute_atoms*3))

      ALLOCATE(helium%rtmp_p_ndim_1d(solute%p*solute%ndim))
      ALLOCATE(helium%rtmp_p_ndim_np_1d(solute%p*solute%ndim*helium%num_env))
      ALLOCATE(helium%rtmp_p_ndim_2d(solute%p,solute%ndim))

    ELSE
      helium%bead_ratio = 0
      IF (helium%periodic) THEN
        ! this assumes a specific potential (and its ugly):
        x1 = angstrom*0.5_dp*helium%cell_size
        ! 10.8 is in Kelvin, x1 needs to be in Angstrom,
        ! since 2.9673 is in Angstrom
        helium%e_corr = (twopi* &
          helium%density/angstrom**3*10.8_dp*(544850.4_dp* &
          EXP(-13.353384_dp*x1/2.9673_dp)*(2.9673_dp/13.353384_dp)**3*&
          (2.0_dp+2.0_dp*13.353384_dp*x1/2.9673_dp+(13.353384_dp*&
          x1/2.9673_dp)**2)-(((0.1781_dp/7.0_dp*(2.9673_dp/x1)**2+&
          0.4253785_dp/5.0_dp)*(2.9673_dp/x1)**2+1.3732412_dp/3.0_dp)*&
          (2.9673_dp/x1)**3)*2.9673_dp**3))/kelvin
      END IF
    END IF

    ! ALLOCATE temporary arrays
    ALLOCATE(helium%itmp_atoms_1d(helium%atoms))
    ALLOCATE(helium%ltmp_atoms_1d(helium%atoms))
    ALLOCATE(helium%itmp_atoms_np_1d(helium%atoms*helium%num_env))
    ALLOCATE(helium%rtmp_3_np_1d(3*helium%num_env))
    ALLOCATE(helium%rtmp_3_atoms_beads_1d(3*helium%atoms*&
      helium%beads))
    ALLOCATE(helium%rtmp_3_atoms_beads_np_1d(3*helium%atoms*&
      helium%beads*helium%num_env))
    ALLOCATE(helium%ltmp_3_atoms_beads_3d(3,helium%atoms,&
      helium%beads))

    CALL helium_write_setup(helium)

    CALL timestop(handle)

    RETURN
  END SUBROUTINE helium_create

! ***************************************************************************
!> \brief  Releases helium_solvent_type
!> \param helium ...
!> \author hforbert
! *****************************************************************************
  SUBROUTINE helium_release(helium)
    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j

    IF (ASSOCIATED(helium)) THEN
      helium%ref_count=helium%ref_count-1
      IF (helium%ref_count<1) THEN

        ! DEALLOCATE temporary arrays
        DEALLOCATE ( &
          helium%ltmp_3_atoms_beads_3d, &
          helium%rtmp_3_atoms_beads_np_1d, &
          helium%rtmp_3_atoms_beads_1d, &
          helium%rtmp_3_np_1d, &
          helium%itmp_atoms_np_1d, &
          helium%ltmp_atoms_1d, &
          helium%itmp_atoms_1d)

        NULLIFY ( &
          helium%ltmp_3_atoms_beads_3d, &
          helium%rtmp_3_atoms_beads_np_1d, &
          helium%rtmp_3_atoms_beads_1d, &
          helium%rtmp_3_np_1d, &
          helium%itmp_atoms_np_1d, &
          helium%ltmp_atoms_1d, &
          helium%itmp_atoms_1d &
        )

        IF (helium%solute_present) THEN
          ! DEALLOCATE solute-related arrays
          DEALLOCATE ( &
            helium%rtmp_p_ndim_2d, &
            helium%rtmp_p_ndim_np_1d, &
            helium%rtmp_p_ndim_1d, &
            helium%force_inst, &
            helium%force_avrg)
          NULLIFY( &
            helium%rtmp_p_ndim_2d, &
            helium%rtmp_p_ndim_np_1d, &
            helium%rtmp_p_ndim_1d, &
            helium%force_inst, &
            helium%force_avrg )
        END IF

        IF ( helium%rho_present ) THEN
          DEALLOCATE( &
            helium%rho_rstr, &
            helium%rho_accu, &
            helium%rho_inst, &
            helium%rho_incr)
          NULLIFY(&
            helium%rho_rstr, &
            helium%rho_accu, &
            helium%rho_inst, &
            helium%rho_incr)
          ! DEALLOCATE everything
          DEALLOCATE(helium%rho_property(rho_atom_number)%filename_suffix)
          DEALLOCATE(helium%rho_property(rho_atom_number)%component_name)
          DEALLOCATE(helium%rho_property(rho_atom_number)%component_index)
          NULLIFY(helium%rho_property(rho_atom_number)%filename_suffix)
          NULLIFY(helium%rho_property(rho_atom_number)%component_name)
          NULLIFY(helium%rho_property(rho_atom_number)%component_index)
          DEALLOCATE(helium%rho_property(rho_winding_number)%filename_suffix)
          DEALLOCATE(helium%rho_property(rho_winding_number)%component_name)
          DEALLOCATE(helium%rho_property(rho_winding_number)%component_index)
          NULLIFY(helium%rho_property(rho_winding_number)%filename_suffix)
          NULLIFY(helium%rho_property(rho_winding_number)%component_name)
          NULLIFY(helium%rho_property(rho_winding_number)%component_index)
          DEALLOCATE(helium%rho_property(rho_winding_cycle)%filename_suffix)
          DEALLOCATE(helium%rho_property(rho_winding_cycle)%component_name)
          DEALLOCATE(helium%rho_property(rho_winding_cycle)%component_index)
          NULLIFY(helium%rho_property(rho_winding_cycle)%filename_suffix)
          NULLIFY(helium%rho_property(rho_winding_cycle)%component_name)
          NULLIFY(helium%rho_property(rho_winding_cycle)%component_index)
          DEALLOCATE(helium%rho_property(rho_projected_area)%filename_suffix)
          DEALLOCATE(helium%rho_property(rho_projected_area)%component_name)
          DEALLOCATE(helium%rho_property(rho_projected_area)%component_index)
          NULLIFY(helium%rho_property(rho_projected_area)%filename_suffix)
          NULLIFY(helium%rho_property(rho_projected_area)%component_name)
          NULLIFY(helium%rho_property(rho_projected_area)%component_index)
          DEALLOCATE(helium%rho_property(rho_moment_of_inertia)%filename_suffix)
          DEALLOCATE(helium%rho_property(rho_moment_of_inertia)%component_name)
          DEALLOCATE(helium%rho_property(rho_moment_of_inertia)%component_index)
          NULLIFY(helium%rho_property(rho_moment_of_inertia)%filename_suffix)
          NULLIFY(helium%rho_property(rho_moment_of_inertia)%component_name)
          NULLIFY(helium%rho_property(rho_moment_of_inertia)%component_index)
          DEALLOCATE(helium%rho_property)
          NULLIFY(helium%rho_property)
        END IF

        CALL helium_rdf_release (helium)

        ! DEALLOCATE helium-related arrays
        DEALLOCATE ( &
          helium%atom_plength, &
          helium%plength_inst, &
          helium%plength_avrg, &
          helium%num_accepted, &
          helium%ipmatrix, &
          helium%pmatrix, &
          helium%nmatrix, &
          helium%tmatrix, &
          helium%iperm, &
          helium%permutation, &
          helium%ptable, &
          helium%work, &
          helium%pos)
        NULLIFY( &
          helium%atom_plength, &
          helium%plength_inst, &
          helium%plength_avrg, &
          helium%num_accepted, &
          helium%ipmatrix, &
          helium%pmatrix, &
          helium%nmatrix, &
          helium%tmatrix, &
          helium%iperm, &
          helium%permutation, &
          helium%ptable, &
          helium%work, &
          helium%pos &
        )

          DO i = 1, SIZE ( helium%eij , 1 )
             DO j = 1, SIZE ( helium%eij , 1 )
                CALL spline_data_release(helium%eij(i,j)%spline_data)
                CALL spline_data_release(helium%uij(i,j)%spline_data)
                !TODO: shouldn't that be done in spline_data_release??
                NULLIFY(helium%eij(i,j)%spline_data, &
                        helium%uij(i,j)%spline_data)
             END DO
          END DO

          DEALLOCATE(helium%eij)
          NULLIFY(helium%eij)

          DEALLOCATE(helium%uij)
          NULLIFY(helium%uij)

          CALL delete_rng_stream(helium%rng_stream_uniform)
          CALL delete_rng_stream(helium%rng_stream_gaussian)

         ! deallocate solute-related arrays
          IF (helium%solute_present) THEN
            CALL helium_destroy_int_arr_ptr(helium%solute_i)
            DEALLOCATE(helium%solute_element)
            NULLIFY(helium%solute_element)
          END IF
        
          ! Deallocate everything from the helium_set_solute_indices
          IF (ASSOCIATED(helium%ename)) THEN
             DEALLOCATE(helium%ename)
             NULLIFY(helium%ename)
          END IF

          IF (ASSOCIATED(helium%eid)) THEN
             DEALLOCATE(helium%eid)
             NULLIFY(helium%eid)
          END IF

          DEALLOCATE( helium )

       END IF
    END IF
    RETURN
  END SUBROUTINE helium_release


! ***************************************************************************
!> \brief  Initialize helium data structures.
!> \param helium ...
!> \param pint_env ...
!> \par    History
!>         removed refereces to pint_env_type data structure [lwalewski]
!>         2009-11-10 init/restore coords, perm, RNG and forces [lwalewski]
!> \author hforbert
!> \note   Initializes helium coordinates either as random positions or from
!>         HELIUM%COORD section if it's present in the input file.
!>         Initializes helium permutation state as identity permutation or
!>         from HELIUM%PERM section if it's present in the input file.
! *****************************************************************************
  SUBROUTINE helium_init( helium, pint_env )

    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i
    LOGICAL                                  :: coords_presampled, explicit, &
                                                presample
    TYPE(section_vals_type), POINTER         :: helium_section, sec

    CALL timeset(routineN,handle)

    NULLIFY(helium_section)
    helium_section => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM")

    ! restore RNG state
    NULLIFY(sec)
    sec => section_vals_get_subs_vals(helium_section,"RNG_STATE")
    CALL section_vals_get(sec,explicit=explicit)
    IF ( explicit ) THEN
      CALL helium_rng_restore( helium )
    ELSE
      CALL helium_write_line("RNG state initialized as new.")
    END IF

    ! init/restore permutation state
    NULLIFY(sec)
    sec => section_vals_get_subs_vals(helium_section,"PERM")
    CALL section_vals_get(sec,explicit=explicit)
    IF ( explicit ) THEN
      CALL helium_perm_restore( helium )
    ELSE
      CALL helium_perm_init( helium )
      CALL helium_write_line("Permutation state initialized as identity.")
    END IF

    ! Specify if forces should be obtained as AVG or LAST
    CALL section_vals_val_get(helium_section, "GET_FORCES",&
                              i_val=helium%get_helium_forces)

    ! init/restore coordinates
    NULLIFY(sec)
    sec => section_vals_get_subs_vals(helium_section,"COORD")
    CALL section_vals_get(sec,explicit=explicit)
    IF ( explicit ) THEN
      CALL helium_coord_restore( helium )
    ELSE 
      CALL helium_coord_init (helium)
    END IF
    DO i=1,helium%atoms
      helium%pos(:,i,helium%beads+1)=helium%pos(:,helium%permutation(i),1)
    END DO
    helium%work=helium%pos
    
    ! init center of mass
    IF (helium%solute_present) THEN
      helium%center(:) = pint_com_pos(pint_env)
    ELSE
      IF (helium%periodic) THEN
        helium%center(:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
      ELSE
        helium%center(:) = helium_com(helium)
      END IF
    END IF

    ! Optional helium coordinate presampling:
    CALL section_vals_val_get(helium_section,"PRESAMPLE",&
         l_val=presample)
    coords_presampled = .FALSE.
    IF ( presample ) THEN
      helium%current_step = 0
      CALL helium_sample( helium, pint_env )
      IF (helium%solute_present) helium%force_avrg(:,:) = 0.0_dp
      helium%energy_avrg(:) = 0.0_dp
      helium%plength_avrg(:) = 0.0_dp
      helium%num_accepted(:,:) = 0.0_dp
      ! Reset properties accumulated over presample:
      helium%proarea%accu(:) = 0.0_dp
      helium%prarea2%accu(:) = 0.0_dp
      helium%wnmber2%accu(:) = 0.0_dp
      helium%mominer%accu(:) = 0.0_dp
      IF (helium%rho_present) THEN
        helium%rho_accu(:,:,:,:) = 0.0_dp
      END IF
      IF (helium%rdf_present) THEN
        helium%rdf_accu(:,:,:) = 0.0_dp
      END IF
      coords_presampled = .TRUE.
      CALL helium_write_line("Bead coordinates pre-sampled.")
    END IF

    IF ( helium%solute_present ) THEN
      ! restore helium forces
      NULLIFY(sec)
      sec => section_vals_get_subs_vals(helium_section,"FORCE")
      CALL section_vals_get(sec,explicit=explicit)
      IF ( explicit ) THEN
        IF ( .NOT. coords_presampled ) THEN
          CALL helium_force_restore( helium )
        END IF
      ELSE
        IF ( .NOT. coords_presampled ) THEN
          CALL helium_force_init (helium)
          CALL helium_write_line("Forces on the solute initialized as zero.")
        END IF
      END IF
    END IF

    CALL timestop(handle)

    RETURN
  END SUBROUTINE helium_init


! ***************************************************************************
! Data transfer functions.
!
! These functions manipulate and transfer data between the runtime
! environment and the input structure.
! ***************************************************************************

! ***************************************************************************
!> \brief  Initialize helium coordinates with random positions.
!> \param helium ...
!> \date   2009-11-09
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_coord_init(helium)
    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_coord_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ib, ic
    REAL(KIND=dp)                            :: r1, r2

    DO ia = 1, helium%atoms
      DO ic = 1, 3
        r1 = next_random_number(helium%rng_stream_uniform)
        r1 = r1 * helium%cell_size
        DO ib = 1, helium%beads
!TODO use thermal gaussian instead
          r2 = next_random_number(helium%rng_stream_uniform)
          helium%pos(ic,ia,ib) = r1+0.1_dp*r2
        END DO
      END DO
      DO ib = 1, helium%beads
        CALL helium_pbc( helium, helium%pos(:,ia,ib) )
      END DO
    END DO

    ! store positions at time slice nbeads+1 (rperez):
    DO ia = 1, helium%atoms
       helium%pos(:,ia,helium%beads+1) = helium%pos(:,ia,1)
    END DO
    ! initialize work array (rperez):
    helium%work=helium%pos

    RETURN
  END SUBROUTINE helium_coord_init

! ***************************************************************************
!> \brief  Restore coordinates from the input structure.
!> \param helium ...
!> \date   2009-11-09
!> \par    History
!>         2010-07-22 accomodate additional cpus in the runtime wrt the
!>                    restart [lwalewski]
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_coord_restore( helium )
    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_coord_restore', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: err_str, stmp
    INTEGER                                  :: actlen, msglen, &
                                                num_env_restart, offset
    LOGICAL, DIMENSION(:, :, :), POINTER     :: m
    REAL(KIND=dp), DIMENSION(:), POINTER     :: message
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: f
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_get_default_logger()

    ! assign the pointer to the memory location of the input structure, where
    ! the coordinates are stored
    NULLIFY(message)
    CALL section_vals_val_get( helium%input, &
         "MOTION%PINT%HELIUM%COORD%_DEFAULT_KEYWORD_", &
         r_vals=message)

    ! check that the number of values in the input match the current runtime
    actlen = SIZE(message)
    num_env_restart = actlen / helium%atoms / helium%beads / 3

    ! distribute coordinates over processors (no message passing)
    msglen = helium%atoms * helium%beads * 3
    offset = msglen * MOD( logger%para_env%mepos, num_env_restart )
    NULLIFY(m,f)
    ALLOCATE(m(3,helium%atoms,helium%beads))
    ALLOCATE(f(3,helium%atoms,helium%beads))
    m(:,:,:) = .TRUE.
    f(:,:,:) = 0.0_dp
    helium%pos(:,:,1:helium%beads) = UNPACK(message(offset+1:offset+msglen), MASK=m, FIELD=f )
    DEALLOCATE(f,m)

    IF ( num_env_restart .NE. logger%para_env%num_pe ) THEN
      err_str = "Reading bead coordinates from the input file."
      CALL helium_write_line(err_str)
      err_str = "Number of environments in the restart...: '"
      stmp = ""
      WRITE(stmp,*) num_env_restart
      err_str = TRIM(ADJUSTL(err_str)) // TRIM(ADJUSTL(stmp)) // "'."
      CALL helium_write_line(err_str)
      err_str = "Number of current run time environments.: '"
      stmp = ""
      WRITE(stmp,*) logger%para_env%num_pe
      err_str = TRIM(ADJUSTL(err_str)) // TRIM(ADJUSTL(stmp)) // "'."
      CALL helium_write_line(err_str)
      IF ( num_env_restart .LT. logger%para_env%num_pe ) THEN
        err_str = "Replicated bead coordinates from the restarted environments."
        CALL helium_write_line(err_str)
      END IF
      IF ( num_env_restart .GT. logger%para_env%num_pe ) THEN
        err_str = "Dropped bead coordinates from some restarted environments."
        CALL helium_write_line(err_str)
      END IF
      err_str = "Done."
      CALL helium_write_line(err_str)
    ELSE
      CALL helium_write_line("Bead coordinates read from the input file.")
    END IF

    NULLIFY(message)

    RETURN
  END SUBROUTINE helium_coord_restore

! ***************************************************************************
!> \brief  Initialize forces exerted on the solute
!> \param helium ...
!> \date   2009-11-10
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_force_init (helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_force_init', &
      routineP = moduleN//':'//routineN

    IF ( helium%solute_present ) THEN
      helium%force_avrg(:,:) = 0.0_dp
      helium%force_inst(:,:) = 0.0_dp
    END IF

    RETURN
  END SUBROUTINE helium_force_init


! ***************************************************************************
!> \brief  Restore forces from the input structure to the runtime environment.
!> \param helium ...
!> \date   2009-11-10
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_force_restore( helium )
    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_force_restore', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: err_str, stmp
    INTEGER                                  :: actlen, msglen
    LOGICAL, DIMENSION(:, :), POINTER        :: m
    REAL(KIND=dp), DIMENSION(:), POINTER     :: message
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: f

    ! assign the pointer to the memory location of the input structure, where
    ! the forces are stored
    NULLIFY(message)
    CALL section_vals_val_get( helium%input, &
         "MOTION%PINT%HELIUM%FORCE%_DEFAULT_KEYWORD_", &
         r_vals=message)

    ! check if the destination array has correct size
    msglen = helium%solute_atoms * helium%solute_beads * 3
    actlen = SIZE(helium%force_avrg)
    err_str = "Invalid size of helium%force_avrg array: actual '"
    stmp = ""
    WRITE(stmp,*) actlen
    err_str = TRIM(ADJUSTL(err_str)) // &
              TRIM(ADJUSTL(stmp)) // "' but expected '"
    stmp = ""
    WRITE(stmp,*) msglen
    IF (actlen/=msglen) THEN
       err_str = TRIM(ADJUSTL(err_str)) // &
                 TRIM(ADJUSTL(stmp)) // "'."
       CPABORT(err_str)
    END IF

    ! restore forces on all processors (no message passing)
    NULLIFY(m,f)
    ALLOCATE(m(helium%solute_beads,helium%solute_atoms*3))
    ALLOCATE(f(helium%solute_beads,helium%solute_atoms*3))
    m(:,:) = .TRUE.
    f(:,:) = 0.0_dp
    helium%force_avrg(:,:) = UNPACK(message(1:msglen), MASK=m, FIELD=f )
    helium%force_inst(:,:) = 0.0_dp
    DEALLOCATE(f,m)

    CALL helium_write_line("Forces on the solute read from the input file.")

    NULLIFY(message)

    RETURN
  END SUBROUTINE helium_force_restore

! ***************************************************************************
!> \brief  Initialize the permutation state.
!> \param helium ...
!> \date   2009-11-05
!> \author Lukasz Walewski
!> \note   Assign the identity permutation at each processor. Inverse
!>         permutation array gets assigned as well.
! *****************************************************************************
  SUBROUTINE helium_perm_init( helium )
    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_perm_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia

    DO ia = 1, helium%atoms
      helium%permutation(ia) = ia
      helium%iperm(ia) = ia
    END DO

    RETURN
  END SUBROUTINE helium_perm_init

! ***************************************************************************
!> \brief  Restore permutation state from the input structre.
!> \param helium ...
!> \date   2009-11-05
!> \par    History
!>         2010-07-22 accomodate additional cpus in the runtime wrt the
!>                    restart [lwalewski]
!> \author Lukasz Walewski
!> \note   Transfer permutation state from the input tree to the runtime
!>         data structures on each processor. Inverse permutation array is
!>         recalculated according to the restored permutation state.
! *****************************************************************************
  SUBROUTINE helium_perm_restore( helium )
    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_perm_restore', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: err_str, stmp
    INTEGER                                  :: actlen, ia, ic, msglen, &
                                                num_env_restart, offset
    INTEGER, DIMENSION(:), POINTER           :: message
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_get_default_logger()

    ! assign the pointer to the memory location of the input structure, where
    ! the permutation state is stored
    NULLIFY(message)
    CALL section_vals_val_get( helium%input, &
         "MOTION%PINT%HELIUM%PERM%_DEFAULT_KEYWORD_", &
         i_vals=message)

    ! check the number of environments presumably stored in the restart
    actlen = SIZE(message)
    num_env_restart = actlen / helium%atoms
!TODO maybe add some sanity checks here:
! is num_env_restart integer ?
! is num_env_restart == helium%num_env_restart ?

    ! distribute permutation state over processors (no message passing)
    msglen = helium%atoms
    offset = msglen * MOD( logger%para_env%mepos, num_env_restart )
    helium%permutation(:) = message(offset+1:offset+msglen)

    IF ( num_env_restart .NE. logger%para_env%num_pe ) THEN
      err_str = "Reading permutation state from the input file."
      CALL helium_write_line(err_str)
      err_str = "Number of environments in the restart...: '"
      stmp = ""
      WRITE(stmp,*) num_env_restart
      err_str = TRIM(ADJUSTL(err_str)) // TRIM(ADJUSTL(stmp)) // "'."
      CALL helium_write_line(err_str)
      err_str = "Number of current run time environments.: '"
      stmp = ""
      WRITE(stmp,*) logger%para_env%num_pe
      err_str = TRIM(ADJUSTL(err_str)) // TRIM(ADJUSTL(stmp)) // "'."
      CALL helium_write_line(err_str)
      IF ( num_env_restart .LT. logger%para_env%num_pe ) THEN
        err_str = "Replicated permutation state from the restarted environments."
        CALL helium_write_line(err_str)
      END IF
      IF ( num_env_restart .GT. logger%para_env%num_pe ) THEN
        err_str = "Dropped permutation state from some restarted environments."
        CALL helium_write_line(err_str)
      END IF
      err_str = "Done."
      CALL helium_write_line(err_str)
    ELSE
      CALL helium_write_line("Permutation state read from the input file.")
    END IF

    ! recalculate the inverse permutation array
    helium%iperm(:) = 0
    ic = 0
    DO ia = 1, msglen
      IF ((helium%permutation(ia)>0).AND.(helium%permutation(ia)<=msglen)) THEN
        helium%iperm(helium%permutation(ia)) = ia
        ic = ic + 1
      END IF
    END DO
    err_str = "Invalid HELIUM%PERM state: some numbers not within (1,"
    stmp = ""
    WRITE(stmp,*) msglen
    IF (ic/=msglen) THEN 
       err_str = TRIM(ADJUSTL(err_str)) // &
                 TRIM(ADJUSTL(stmp)) // ")."
       CPABORT(err_str)
    END IF

    NULLIFY(message)

    RETURN
  END SUBROUTINE helium_perm_restore


! ***************************************************************************
!> \brief  Restore averages from the input structure
!> \param helium ...
!> \date   2014-06-25
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_averages_restore(helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_averages_restore', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: msglen, num_env_restart, &
                                                offset
    LOGICAL                                  :: explicit
    REAL(KIND=dp), DIMENSION(:), POINTER     :: message
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_get_default_logger()

    ! restore projected area
    CALL section_vals_val_get( helium%input, &
      "MOTION%PINT%HELIUM%AVERAGES%PROJECTED_AREA", &
      explicit=explicit)
    IF (explicit) THEN
      NULLIFY(message)
      CALL section_vals_val_get( helium%input, &
        "MOTION%PINT%HELIUM%AVERAGES%PROJECTED_AREA", &
        r_vals=message)
      num_env_restart = SIZE(message) / 3 ! apparent number of environments
      msglen = 3
      offset = msglen * MOD( logger%para_env%mepos, num_env_restart )
      helium%proarea%rstr(:) = message(offset+1:offset+msglen)
    ELSE
      helium%proarea%rstr(:) = 0.0_dp
    END IF

    ! restore projected area squared
    CALL section_vals_val_get( helium%input, &
      "MOTION%PINT%HELIUM%AVERAGES%PROJECTED_AREA_2", &
      explicit=explicit)
    IF (explicit) THEN
      NULLIFY(message)
      CALL section_vals_val_get( helium%input, &
        "MOTION%PINT%HELIUM%AVERAGES%PROJECTED_AREA_2", &
        r_vals=message)
      num_env_restart = SIZE(message) / 3 ! apparent number of environments
      msglen = 3
      offset = msglen * MOD( logger%para_env%mepos, num_env_restart )
      helium%prarea2%rstr(:) = message(offset+1:offset+msglen)
    ELSE
      helium%prarea2%rstr(:) = 0.0_dp
    END IF

    ! restore winding number squared
    CALL section_vals_val_get( helium%input, &
      "MOTION%PINT%HELIUM%AVERAGES%WINDING_NUMBER_2", &
      explicit=explicit)
    IF (explicit) THEN
      NULLIFY(message)
      CALL section_vals_val_get( helium%input, &
        "MOTION%PINT%HELIUM%AVERAGES%WINDING_NUMBER_2", &
        r_vals=message)
      num_env_restart = SIZE(message) / 3 ! apparent number of environments
      msglen = 3
      offset = msglen * MOD( logger%para_env%mepos, num_env_restart )
      helium%wnmber2%rstr(:) = message(offset+1:offset+msglen)
    ELSE
      helium%wnmber2%rstr(:) = 0.0_dp
    END IF

    ! restore moment of inertia
    CALL section_vals_val_get( helium%input, &
      "MOTION%PINT%HELIUM%AVERAGES%MOMENT_OF_INERTIA", &
      explicit=explicit)
    IF (explicit) THEN
      NULLIFY(message)
      CALL section_vals_val_get( helium%input, &
        "MOTION%PINT%HELIUM%AVERAGES%MOMENT_OF_INERTIA", &
        r_vals=message)
      num_env_restart = SIZE(message) / 3 ! apparent number of environments
      msglen = 3
      offset = msglen * MOD( logger%para_env%mepos, num_env_restart )
      helium%mominer%rstr(:) = message(offset+1:offset+msglen)
    ELSE
      helium%mominer%rstr(:) = 0.0_dp
    END IF

    IF ( helium%rdf_present ) THEN
      CALL helium_rdf_restore( helium )
    END IF

    IF (helium%rho_present) THEN
      ! restore densities
      CALL helium_rho_restore( helium )
    END IF

    ! get the weighting factor
    CALL section_vals_val_get( &
      helium%input, &
      "MOTION%PINT%HELIUM%AVERAGES%IWEIGHT", &
      i_val=helium%averages_iweight)

    ! set the flag indicating whether the averages have been restarted
    CALL section_vals_val_get( &
      helium%input, &
      "EXT_RESTART%RESTART_HELIUM_AVERAGES", &
      l_val=helium%averages_restarted)

    RETURN
  END SUBROUTINE helium_averages_restore


! ***************************************************************************
!> \brief  Create RNG streams and initialize their state.
!> \param helium ...
!> \date   2009-11-04
!> \author Lukasz Walewski
!> \note   TODO: This function shouldn't create (allocate) objects! Only 
!>         initialization, i.e. setting the seed values etc, should be done 
!>         here, allocation should be moved to helium_create 
! *****************************************************************************
  SUBROUTINE helium_rng_init( helium )
    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_rng_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: rank
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(rng_stream_type), POINTER           :: next_rngs, prev_rngs

    NULLIFY(logger)
    logger => cp_get_default_logger()

    ! Create two RNG streams on each processor and initialize them so, that
    ! each processor gets unique RN sequences.
    NULLIFY(prev_rngs, next_rngs)
    NULLIFY(helium%rng_stream_uniform,helium%rng_stream_gaussian)

    ! Create two RNG strems at each processor by seeding one from the
    ! other. Then, on rank 0: save the pointers to both structures, while
    ! on all other ranks: delete the first structure and keep the second
    ! one to seed next RNG stream to be created.
    CALL create_rng_stream(prev_rngs,&
         name="helium_rns_uniform",&
         distribution_type=UNIFORM,&
         extended_precision=.TRUE.)
    IF (logger%para_env%mepos .EQ. 0) THEN
      helium%rng_stream_uniform => prev_rngs
    END IF
    CALL create_rng_stream(next_rngs,&
         name="helium_rns_gaussian",&
         last_rng_stream=prev_rngs,&
         distribution_type=GAUSSIAN,&
         extended_precision=.TRUE.)
    IF (logger%para_env%mepos .EQ. 0) THEN
      helium%rng_stream_gaussian => next_rngs
      NULLIFY(prev_rngs)
    ELSE
      CALL delete_rng_stream(prev_rngs)
      prev_rngs => next_rngs
    END IF
    NULLIFY(next_rngs)

    ! At all ranks higher than 0 keep creating new RNG streams one from
    ! the other (so they are all different) and on each rank cut this
    ! process at different stage.
    DO rank = 1, logger%para_env%mepos
      CALL create_rng_stream(next_rngs,&
           name="helium_rns_uniform",&
           last_rng_stream=prev_rngs,&
           distribution_type=UNIFORM,&
           extended_precision=.TRUE.)
      IF ( logger%para_env%mepos .EQ. rank ) THEN
        helium%rng_stream_uniform => next_rngs
      END IF
      CALL delete_rng_stream(prev_rngs)
      prev_rngs => next_rngs
      NULLIFY(next_rngs)
      CALL create_rng_stream(next_rngs,&
           name="helium_rns_gaussian",&
           last_rng_stream=prev_rngs,&
           distribution_type=GAUSSIAN,&
           extended_precision=.TRUE.)
      IF ( logger%para_env%mepos .EQ. rank ) THEN
        helium%rng_stream_gaussian => next_rngs
        NULLIFY(prev_rngs)
      ELSE
        CALL delete_rng_stream(prev_rngs)
        prev_rngs => next_rngs
      END IF
      NULLIFY(next_rngs)
    END DO

    RETURN
  END SUBROUTINE helium_rng_init

! ***************************************************************************
!> \brief  Restore RNG state from the input structure.
!> \param helium ...
!> \date   2009-11-04
!> \par    History
!>         2010-07-22 Create new rng streams if more cpus available in the
!>         runtime than in the restart [lwalewski]
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_rng_restore( helium )
    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_rng_restore', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: err_str, stmp
    INTEGER                                  :: actlen, msglen, &
                                                num_env_restart, offset, rank
    LOGICAL                                  :: lbf
    LOGICAL, DIMENSION(3, 2)                 :: m
    REAL(KIND=dp)                            :: bf, bu
    REAL(KIND=dp), DIMENSION(3, 2)           :: bg, cg, f, ig
    REAL(KIND=dp), DIMENSION(:), POINTER     :: message
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(rng_stream_type), POINTER           :: next_rngs, prev_rngs

    NULLIFY(logger)
    logger => cp_get_default_logger()

    ! assign the pointer to the memory location of the input structure
    ! where the RNG state is stored
    NULLIFY(message)
    CALL section_vals_val_get( helium%input, &
         "MOTION%PINT%HELIUM%RNG_STATE%_DEFAULT_KEYWORD_", &
         r_vals=message)

    ! check the number of environments presumably stored in the restart
    actlen = SIZE(message)
    num_env_restart = actlen / 40

    IF ( logger%para_env%mepos .LT. num_env_restart ) THEN

      ! unpack the buffer at each processor, set RNG state (no message passing)
      msglen = 40
      offset = msglen * logger%para_env%mepos
      m(:,:) = .TRUE.
      f(:,:) = 0.0_dp
      bg(:,:) = UNPACK(message(offset+1:offset+6), MASK=m, FIELD=f )
      cg(:,:) = UNPACK(message(offset+7:offset+12), MASK=m, FIELD=f )
      ig(:,:) = UNPACK(message(offset+13:offset+18), MASK=m, FIELD=f )
      bf = message(offset+19)
      bu = message(offset+20)
      IF ( bf .GT. 0) THEN
        lbf = .TRUE.
      ELSE
        lbf = .FALSE.
      END IF
      CALL set_rng_stream(helium%rng_stream_uniform,bg=bg,cg=cg,ig=ig,&
           buffer=bu,buffer_filled=lbf)
      bg(:,:) = UNPACK(message(offset+21:offset+26), MASK=m, FIELD=f )
      cg(:,:) = UNPACK(message(offset+27:offset+32), MASK=m, FIELD=f )
      ig(:,:) = UNPACK(message(offset+33:offset+38), MASK=m, FIELD=f )
      bf = message(offset+39)
      bu = message(offset+40)
      IF ( bf .GT. 0) THEN
        lbf = .TRUE.
      ELSE
        lbf = .FALSE.
      END IF
      CALL set_rng_stream(helium%rng_stream_gaussian,bg=bg,cg=cg,ig=ig,&
           buffer=bu,buffer_filled=lbf)

    ELSE
      ! On processors that did not receive rng state from the restart file
      ! delete rng streams created in helium_rng_init and create them
      ! anew, as they have been initialized with default initial state. Here
      ! the sequence of rng streams starts from the last stream from the
      ! restart file, each stream is initialized from the previously created
      ! one.

      CALL delete_rng_stream(helium%rng_stream_uniform)
      CALL delete_rng_stream(helium%rng_stream_gaussian)
      NULLIFY(prev_rngs, next_rngs)
      NULLIFY(helium%rng_stream_uniform,helium%rng_stream_gaussian)

      ! take the last uniform rng stream from the restart as a starting point
      msglen = 40
      offset = msglen * ( num_env_restart - 1 )
      m(:,:) = .TRUE.
      f(:,:) = 0.0_dp
      bg(:,:) = UNPACK(message(offset+1:offset+6), MASK=m, FIELD=f )
      cg(:,:) = UNPACK(message(offset+7:offset+12), MASK=m, FIELD=f )
      ig(:,:) = UNPACK(message(offset+13:offset+18), MASK=m, FIELD=f )
      bf = message(offset+19)
      bu = message(offset+20)
      IF ( bf .GT. 0) THEN
        lbf = .TRUE.
      ELSE
        lbf = .FALSE.
      END IF
      ! create new uniform rng from scratch
      CALL create_rng_stream(prev_rngs,&
           name="helium_rns_uniform",&
           distribution_type=UNIFORM,&
           extended_precision=.TRUE.)
      ! set it to the last uniform rng stream from the restart
      CALL set_rng_stream(prev_rngs,bg=bg,cg=cg,ig=ig,&
           buffer=bu,buffer_filled=lbf)
      ! use this on the first non-restarted environment as the new rng
      IF ( logger%para_env%mepos .EQ. num_env_restart ) THEN
        helium%rng_stream_uniform => prev_rngs
      END IF
      ! unpack the last gaussian rng stream from the restart
      bg(:,:) = UNPACK(message(offset+21:offset+26), MASK=m, FIELD=f )
      cg(:,:) = UNPACK(message(offset+27:offset+32), MASK=m, FIELD=f )
      ig(:,:) = UNPACK(message(offset+33:offset+38), MASK=m, FIELD=f )
      bf = message(offset+39)
      bu = message(offset+40)
      IF ( bf .GT. 0) THEN
        lbf = .TRUE.
      ELSE
        lbf = .FALSE.
      END IF
      ! create new gaussian rng stream from scratch
      CALL create_rng_stream(next_rngs,&
           name="helium_rns_gaussian",&
           last_rng_stream=prev_rngs,&
           distribution_type=GAUSSIAN,&
           extended_precision=.TRUE.)
      ! set it to the last gaussian rng stream from the restart
      CALL set_rng_stream(next_rngs,bg=bg,cg=cg,ig=ig,&
           buffer=bu,buffer_filled=lbf)
      ! use this on the first non-restarted environment as the new gaussian rng
      IF ( logger%para_env%mepos .EQ. num_env_restart ) THEN
        helium%rng_stream_gaussian => next_rngs
        NULLIFY(prev_rngs)
      ELSE
        CALL delete_rng_stream(prev_rngs)
        prev_rngs => next_rngs
      END IF
      NULLIFY(next_rngs)

      ! At all ranks higher than num_env_restart keep creating new RNG streams
      ! one from the other (so they are all different) and on each rank cut
      ! this process at different stage.
      DO rank = num_env_restart + 1, logger%para_env%mepos
        CALL create_rng_stream(next_rngs,&
             name="helium_rns_uniform",&
             last_rng_stream=prev_rngs,&
             distribution_type=UNIFORM,&
             extended_precision=.TRUE.)
        IF ( logger%para_env%mepos .EQ. rank ) THEN
          helium%rng_stream_uniform => next_rngs
        END IF
        CALL delete_rng_stream(prev_rngs)
        prev_rngs => next_rngs
        NULLIFY(next_rngs)
        CALL create_rng_stream(next_rngs,&
             name="helium_rns_gaussian",&
             last_rng_stream=prev_rngs,&
             distribution_type=GAUSSIAN,&
             extended_precision=.TRUE.)
        IF ( logger%para_env%mepos .EQ. rank ) THEN
          helium%rng_stream_gaussian => next_rngs
          NULLIFY(prev_rngs)
        ELSE
          CALL delete_rng_stream(prev_rngs)
          prev_rngs => next_rngs
        END IF
        NULLIFY(next_rngs)
      END DO

    END IF

    ! say what has been done
    IF ( num_env_restart .NE. logger%para_env%num_pe ) THEN
      err_str = "Reading RNG state from the input file."
      CALL helium_write_line(err_str)
      err_str = "Number of environments in the restart...: '"
      stmp = ""
      WRITE(stmp,*) num_env_restart
      err_str = TRIM(ADJUSTL(err_str)) // TRIM(ADJUSTL(stmp)) // "'."
      CALL helium_write_line(err_str)
      err_str = "Number of current run time environments.: '"
      stmp = ""
      WRITE(stmp,*) logger%para_env%num_pe
      err_str = TRIM(ADJUSTL(err_str)) // TRIM(ADJUSTL(stmp)) // "'."
      CALL helium_write_line(err_str)
      IF ( num_env_restart .LT. logger%para_env%num_pe ) THEN
        err_str = "Created some new RNGs from the restarted environments."
        CALL helium_write_line(err_str)
      END IF
      IF ( num_env_restart .GT. logger%para_env%num_pe ) THEN
        err_str = "Dropped RNG state from some restarted environments."
        CALL helium_write_line(err_str)
      END IF
      err_str = "Done."
      CALL helium_write_line(err_str)
    ELSE
      CALL helium_write_line("RNG state read from the input file.")
    END IF

    NULLIFY(message)

    RETURN
  END SUBROUTINE helium_rng_restore


! ***************************************************************************
!> \brief  Create the RDF related data structures and initialize
!> \param helium ...
!> \date   2014-02-25
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_rdf_init(helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_rdf_init', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: file_name
    CHARACTER(len=default_string_length)     :: err_str, stmp
    INTEGER                                  :: ifirst, ii, ij, ik
    LOGICAL                                  :: explicit
    REAL(KIND=dp), DIMENSION(:), POINTER     :: coords
    TYPE(cp_logger_type), POINTER            :: logger

    IF (helium%solute_present) THEN
      ! get number of centers from solute
      helium%rdf_num_ctr = helium%solute_atoms
    ELSE
      ! use center of the unit cell or COM of the He droplet
      helium%rdf_num_ctr = 1
    END IF
    ifirst = helium%rdf_num_ctr

    ! read auxiliary centers from file
    NULLIFY(logger)
    logger => cp_get_default_logger()
    CALL section_vals_val_get( &
      helium%input, &
      "MOTION%PINT%HELIUM%RDF%CENTERS_FILE_NAME",&
      c_val=file_name, &
      explicit=explicit)
    NULLIFY(coords)
    IF (explicit) THEN
      CALL section_vals_val_get( &
        helium%input, &
        "MOTION%PINT%HELIUM%RDF%CENTERS_FILE_NAME", &
        c_val=file_name)
      CALL helium_read_xyz(coords,file_name,logger%para_env)
    END IF

    ! increment number of centers if necessary
    IF (ASSOCIATED(coords)) THEN
      helium%rdf_num_ctr = helium%rdf_num_ctr + SIZE(coords)/3
    END IF

    ! set the flag for RDF and either proceed or return
    IF (helium%rdf_num_ctr>0) THEN
      helium%rdf_present = .TRUE.
    ELSE
      helium%rdf_present = .FALSE.
      RETURN
    END IF

    ! allocate & store auxiliary centers (align to the end of array)
    ALLOCATE(helium%rdf_centers(3*helium%rdf_num_ctr))
    helium%rdf_centers(:) = 0.0_dp
    IF (ASSOCIATED(coords)) THEN
      ii = 3*ifirst+1
      ij = 3*helium%rdf_num_ctr
      helium%rdf_centers(ii:ij) = coords(:)
      DEALLOCATE(coords)
      NULLIFY(coords)
    END IF

    helium%rdf_num = 4 ! see also helium%rho_num_def

    ! set the maximum RDF range
    CALL section_vals_val_get( &
      helium%input, &
      "MOTION%PINT%HELIUM%RDF%MAXR", &
      explicit=explicit)
    IF (explicit) THEN
      ! use the value explicitly set in the input
      CALL section_vals_val_get( &
        helium%input, &
        "MOTION%PINT%HELIUM%RDF%MAXR",&
        r_val=helium%rdf_maxr)
    ELSE
      ! use the default value
      CALL section_vals_val_get( &
        helium%input, &
        "MOTION%PINT%HELIUM%DROPLET_RADIUS", &
        explicit=explicit)
      IF (explicit) THEN
        ! use the droplet radius
        helium%rdf_maxr = helium%droplet_radius
      ELSE
        ! use cell_size and cell_shape
        ! (they are set regardless of us being periodic or not)
        SELECT CASE (helium%cell_shape)
          CASE (helium_cell_shape_cube)
            helium%rdf_maxr = helium%cell_size / 2.0_dp
          CASE (helium_cell_shape_octahedron)
            helium%rdf_maxr = helium%cell_size * SQRT(3.0_dp) / 4.0_dp
          CASE DEFAULT
            helium%rdf_maxr = 0.0_dp
            WRITE(stmp,*) helium%cell_shape
            err_str = "cell shape unknown (" // TRIM(ADJUSTL(stmp)) // ")"
            CPABORT(err_str)
        END SELECT
      END IF
    END IF

    ! get number of bins and set bin spacing
    CALL section_vals_val_get( &
      helium%input, &
      "MOTION%PINT%HELIUM%RDF%NBIN", &
      i_val=helium%rdf_nbin)
    helium%rdf_delr = helium%rdf_maxr / REAL(helium%rdf_nbin,dp)

    ! get the weighting factor
    CALL section_vals_val_get( &
      helium%input, &
      "MOTION%PINT%HELIUM%AVERAGES%IWEIGHT", &
      i_val=helium%rdf_iweight)

    ! allocate and initialize memory for RDF storage
    ii = helium%rdf_num
    ij = helium%rdf_nbin
    ik = helium%rdf_num_ctr
    ALLOCATE(helium%rdf_inst(ii,ij,ik))
    ALLOCATE(helium%rdf_accu(ii,ij,ik))
    ALLOCATE(helium%rdf_rstr(ii,ij,ik))
    helium%rdf_inst(:,:,:) = 0.0_dp
    helium%rdf_accu(:,:,:) = 0.0_dp
    helium%rdf_rstr(:,:,:) = 0.0_dp

    RETURN
  END SUBROUTINE helium_rdf_init


! ***************************************************************************
!> \brief  Restore the RDFs from the input structure
!> \param helium ...
!> \date   2011-06-22
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_rdf_restore (helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_rdf_restore', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: stmp1, stmp2
    CHARACTER(len=max_line_length)           :: err_str
    INTEGER                                  :: ii, ij, ik, itmp, msglen
    LOGICAL                                  :: explicit, ltmp
    LOGICAL, DIMENSION(:, :, :), POINTER     :: m
    REAL(KIND=dp), DIMENSION(:), POINTER     :: message
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: f

    CALL section_vals_val_get( helium%input, &
      "MOTION%PINT%HELIUM%AVERAGES%RDF", &
      explicit=explicit)
    IF (explicit) THEN
      NULLIFY(message)
      CALL section_vals_val_get( helium%input, &
        "MOTION%PINT%HELIUM%AVERAGES%RDF", &
        r_vals=message)
      msglen = SIZE(message)
      itmp = SIZE(helium%rdf_rstr)
      ltmp = ( msglen .EQ. itmp )
      IF ( .NOT. ltmp ) THEN
        stmp1 = ""
        WRITE(stmp1,*) msglen
        stmp2 = ""
        WRITE(stmp2,*) itmp
        err_str = "Size of the RDF array in the input (" // &
                  TRIM(ADJUSTL(stmp1)) // &
                  ") .NE. that in the runtime (" // &
                  TRIM(ADJUSTL(stmp2)) // ")."
        CPABORT(err_str)
      END IF
    ELSE
      RETURN
    END IF

    ii = helium%rdf_num
    ij = helium%rdf_nbin
    ik = helium%rdf_num_ctr
    NULLIFY(m,f)
    ALLOCATE(m(ii,ij,ik))
    ALLOCATE(f(ii,ij,ik))
    m(:,:,:) = .TRUE.
    f(:,:,:) = 0.0_dp

    helium%rdf_rstr(:,:,:) = UNPACK(message(1:msglen), MASK=m, FIELD=f)
    CALL section_vals_val_get( helium%input, &
        "MOTION%PINT%HELIUM%AVERAGES%IWEIGHT", &
        i_val=helium%rdf_iweight)

    DEALLOCATE(f,m)
    NULLIFY(message)

    RETURN
  END SUBROUTINE helium_rdf_restore


! ***************************************************************************
!> \brief  Release/deallocate RDF related data structures
!> \param helium ...
!> \date   2014-02-25
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_rdf_release (helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_rdf_release', &
      routineP = moduleN//':'//routineN

    IF (helium%rdf_present) THEN

      DEALLOCATE( &
        helium%rdf_centers, &
        helium%rdf_rstr, &
        helium%rdf_accu, &
        helium%rdf_inst)

      NULLIFY( &
        helium%rdf_centers, &
        helium%rdf_rstr, &
        helium%rdf_accu, &
        helium%rdf_inst )

    END IF

    RETURN
  END SUBROUTINE helium_rdf_release


! ***************************************************************************
!> \brief  Check whether property <prop> is activated in the input structure
!> \param helium ...
!> \param prop ...
!> \retval is_active ...
!> \date   2014-06-26
!> \author Lukasz Walewski
!> \note   The property is controlled by two items in the input structure,
!>         the printkey and the control section. Two settings result in
!>         the property being considered active:
!>         1. printkey is on at the given print level
!>         2. control section is explicit and on
!>         If the property is considered active it should be calculated
!>         and printed through out the run.
! *****************************************************************************
  FUNCTION helium_property_active(helium,prop) RESULT(is_active)

    TYPE(helium_solvent_type), POINTER       :: helium
    CHARACTER(len=*)                         :: prop
    LOGICAL                                  :: is_active

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_property_active', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: input_path
    INTEGER                                  :: print_level
    LOGICAL                                  :: explicit, is_on
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key, section

    NULLIFY(logger)
    logger => cp_get_default_logger()

    ! if the printkey is on at this runlevel consider prop to be active
    NULLIFY(print_key)
    input_path = "MOTION%PINT%HELIUM%PRINT%" // TRIM(ADJUSTL(prop))
    print_key => section_vals_get_subs_vals( &
      helium%input, &
      input_path)
    is_on = cp_printkey_is_on( &
      iteration_info=logger%iter_info, &
      print_key=print_key)
    IF (is_on) THEN
      is_active = .TRUE.
      RETURN
    END IF

    ! if the control section is explicit and on consider prop to be active
    ! and activate the printkey
    is_active = .FALSE.
    NULLIFY(section)
    input_path = "MOTION%PINT%HELIUM%" // TRIM(ADJUSTL(prop))
    section => section_vals_get_subs_vals( &
      helium%input, &
      input_path)
    CALL section_vals_get(section,explicit=explicit)
    IF (explicit) THEN
      ! control section explicitly present, continue checking
      CALL section_vals_val_get( &
        section, &
        "_SECTION_PARAMETERS_", &
        l_val=is_on)
      IF (is_on) THEN
        ! control section is explicit and on, activate the property
        is_active = .TRUE.
        ! activate the corresponding print_level as well
        print_level = logger%iter_info%print_level
        CALL section_vals_val_set( &
          print_key, &
          "_SECTION_PARAMETERS_", &
          i_val=print_level)
      END IF
    END IF

    RETURN
  END FUNCTION helium_property_active


! ***************************************************************************
!> \brief  Create the density related data structures and initialize
!> \param helium ...
!> \date   2014-02-25
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_rho_property_init(helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_rho_property_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nc

    ALLOCATE(helium%rho_property(rho_num))

    helium%rho_property(rho_atom_number)%name = 'Atom number density'
    nc = 1
    helium%rho_property(rho_atom_number)%num_components = nc
    ALLOCATE(helium%rho_property(rho_atom_number)%filename_suffix(nc))
    ALLOCATE(helium%rho_property(rho_atom_number)%component_name(nc))
    ALLOCATE(helium%rho_property(rho_atom_number)%component_index(nc))
    helium%rho_property(rho_atom_number)%filename_suffix(1) = 'an'
    helium%rho_property(rho_atom_number)%component_name(1) = ''
    helium%rho_property(rho_atom_number)%component_index(:) = 0

    helium%rho_property(rho_projected_area)%name = 'Projected area squared density, A*A(r)'
    nc = 3
    helium%rho_property(rho_projected_area)%num_components = nc
    ALLOCATE(helium%rho_property(rho_projected_area)%filename_suffix(nc))
    ALLOCATE(helium%rho_property(rho_projected_area)%component_name(nc))
    ALLOCATE(helium%rho_property(rho_projected_area)%component_index(nc))
    helium%rho_property(rho_projected_area)%filename_suffix(1) = 'pa_x'
    helium%rho_property(rho_projected_area)%filename_suffix(2) = 'pa_y'
    helium%rho_property(rho_projected_area)%filename_suffix(3) = 'pa_z'
    helium%rho_property(rho_projected_area)%component_name(1) = 'component x'
    helium%rho_property(rho_projected_area)%component_name(2) = 'component y'
    helium%rho_property(rho_projected_area)%component_name(3) = 'component z'
    helium%rho_property(rho_projected_area)%component_index(:) = 0

    helium%rho_property(rho_winding_number)%name = 'Winding number squared density, W*W(r)'
    nc = 3
    helium%rho_property(rho_winding_number)%num_components = nc
    ALLOCATE(helium%rho_property(rho_winding_number)%filename_suffix(nc))
    ALLOCATE(helium%rho_property(rho_winding_number)%component_name(nc))
    ALLOCATE(helium%rho_property(rho_winding_number)%component_index(nc))
    helium%rho_property(rho_winding_number)%filename_suffix(1) = 'wn_x'
    helium%rho_property(rho_winding_number)%filename_suffix(2) = 'wn_y'
    helium%rho_property(rho_winding_number)%filename_suffix(3) = 'wn_z'
    helium%rho_property(rho_winding_number)%component_name(1) = 'component x'
    helium%rho_property(rho_winding_number)%component_name(2) = 'component y'
    helium%rho_property(rho_winding_number)%component_name(3) = 'component z'
    helium%rho_property(rho_winding_number)%component_index(:) = 0

    helium%rho_property(rho_winding_cycle)%name = 'Winding number squared density, W^2(r)'
    nc = 3
    helium%rho_property(rho_winding_cycle)%num_components = nc
    ALLOCATE(helium%rho_property(rho_winding_cycle)%filename_suffix(nc))
    ALLOCATE(helium%rho_property(rho_winding_cycle)%component_name(nc))
    ALLOCATE(helium%rho_property(rho_winding_cycle)%component_index(nc))
    helium%rho_property(rho_winding_cycle)%filename_suffix(1) = 'wc_x'
    helium%rho_property(rho_winding_cycle)%filename_suffix(2) = 'wc_y'
    helium%rho_property(rho_winding_cycle)%filename_suffix(3) = 'wc_z'
    helium%rho_property(rho_winding_cycle)%component_name(1) = 'component x'
    helium%rho_property(rho_winding_cycle)%component_name(2) = 'component y'
    helium%rho_property(rho_winding_cycle)%component_name(3) = 'component z'
    helium%rho_property(rho_winding_cycle)%component_index(:) = 0

    helium%rho_property(rho_moment_of_inertia)%name = 'Moment of inertia'
    nc = 3
    helium%rho_property(rho_moment_of_inertia)%num_components = nc
    ALLOCATE(helium%rho_property(rho_moment_of_inertia)%filename_suffix(nc))
    ALLOCATE(helium%rho_property(rho_moment_of_inertia)%component_name(nc))
    ALLOCATE(helium%rho_property(rho_moment_of_inertia)%component_index(nc))
    helium%rho_property(rho_moment_of_inertia)%filename_suffix(1) = 'mi_x'
    helium%rho_property(rho_moment_of_inertia)%filename_suffix(2) = 'mi_y'
    helium%rho_property(rho_moment_of_inertia)%filename_suffix(3) = 'mi_z'
    helium%rho_property(rho_moment_of_inertia)%component_name(1) = 'component x'
    helium%rho_property(rho_moment_of_inertia)%component_name(2) = 'component y'
    helium%rho_property(rho_moment_of_inertia)%component_name(3) = 'component z'
    helium%rho_property(rho_moment_of_inertia)%component_index(:) = 0

    helium%rho_property(:)%is_calculated = .FALSE.

    RETURN
  END SUBROUTINE helium_rho_property_init


! ***************************************************************************
!> \brief  Create the density related data structures and initialize
!> \param helium ...
!> \date   2014-02-25
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_rho_init(helium)

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_rho_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii, itmp, jtmp
    LOGICAL                                  :: explicit, ltmp

    CALL helium_rho_property_init(helium)

    helium%rho_num_act = 0

    ! check for atom number density
    CALL section_vals_val_get(&
      helium%input,&
      "MOTION%PINT%HELIUM%RHO%ATOM_NUMBER",&
      l_val=ltmp)
    IF (ltmp) THEN
      helium%rho_property(rho_atom_number)%is_calculated = .TRUE.
      helium%rho_num_act = helium%rho_num_act + 1
      helium%rho_property(rho_atom_number)%component_index(1) = helium%rho_num_act
    END IF

    ! check for projected area density
    CALL section_vals_val_get(&
      helium%input,&
      "MOTION%PINT%HELIUM%RHO%PROJECTED_AREA_2",&
      l_val=ltmp)
    IF (ltmp) THEN
      helium%rho_property(rho_projected_area)%is_calculated = .TRUE.
      DO ii = 1, helium%rho_property(rho_projected_area)%num_components
        helium%rho_num_act = helium%rho_num_act + 1
        helium%rho_property(rho_projected_area)%component_index(ii) = helium%rho_num_act
      END DO
    END IF

    ! check for winding number density
    CALL section_vals_val_get(&
      helium%input,&
      "MOTION%PINT%HELIUM%RHO%WINDING_NUMBER_2",&
      l_val=ltmp)
    IF (ltmp) THEN
      helium%rho_property(rho_winding_number)%is_calculated = .TRUE.
      DO ii = 1, helium%rho_property(rho_winding_number)%num_components
        helium%rho_num_act = helium%rho_num_act + 1
        helium%rho_property(rho_winding_number)%component_index(ii) = helium%rho_num_act
      END DO
    END IF

    ! check for winding cycle density
    CALL section_vals_val_get(&
      helium%input,&
      "MOTION%PINT%HELIUM%RHO%WINDING_CYCLE_2",&
      l_val=ltmp)
    IF (ltmp) THEN
      helium%rho_property(rho_winding_cycle)%is_calculated = .TRUE.
      DO ii = 1, helium%rho_property(rho_winding_cycle)%num_components
        helium%rho_num_act = helium%rho_num_act + 1
        helium%rho_property(rho_winding_cycle)%component_index(ii) = helium%rho_num_act
      END DO
    END IF

    ! check for moment of inertia density
    CALL section_vals_val_get(&
      helium%input,&
      "MOTION%PINT%HELIUM%RHO%MOMENT_OF_INERTIA",&
      l_val=ltmp)
    IF (ltmp) THEN
      helium%rho_property(rho_moment_of_inertia)%is_calculated = .TRUE.
      DO ii = 1, helium%rho_property(rho_moment_of_inertia)%num_components
        helium%rho_num_act = helium%rho_num_act + 1
        helium%rho_property(rho_moment_of_inertia)%component_index(ii) = helium%rho_num_act
      END DO
    END IF

    ! set the cube dimensions, etc (common to all estimators)
    helium%rho_maxr = helium%cell_size
    CALL section_vals_val_get( &
      helium%input, &
      "MOTION%PINT%HELIUM%RHO%NBIN", &
      i_val=helium%rho_nbin)
    helium%rho_delr = helium%rho_maxr / REAL(helium%rho_nbin,dp)

    ! check for optional estimators based on winding paths
    helium%rho_num_min_len_wdg = 0
    CALL section_vals_val_get( &
      helium%input, &
      "MOTION%PINT%HELIUM%RHO%MIN_CYCLE_LENGTHS_WDG", &
      explicit=explicit)
    IF (explicit) THEN
      NULLIFY(helium%rho_min_len_wdg_vals)
      CALL section_vals_val_get( &
        helium%input, &
        "MOTION%PINT%HELIUM%RHO%MIN_CYCLE_LENGTHS_WDG", &
        i_vals=helium%rho_min_len_wdg_vals)
         itmp = SIZE(helium%rho_min_len_wdg_vals)
      IF (itmp .GT. 0) THEN
        helium%rho_num_min_len_wdg = itmp
        helium%rho_num_act = helium%rho_num_act + itmp
      END IF
    END IF

    ! check for optional estimators based on non-winding paths
    helium%rho_num_min_len_non = 0
    CALL section_vals_val_get( &
      helium%input, &
      "MOTION%PINT%HELIUM%RHO%MIN_CYCLE_LENGTHS_NON", &
      explicit=explicit)
    IF (explicit) THEN
      NULLIFY(helium%rho_min_len_non_vals)
      CALL section_vals_val_get( &
        helium%input, &
        "MOTION%PINT%HELIUM%RHO%MIN_CYCLE_LENGTHS_NON", &
        i_vals=helium%rho_min_len_non_vals)
      itmp = SIZE(helium%rho_min_len_non_vals)
      IF (itmp .GT. 0) THEN
        helium%rho_num_min_len_non = itmp
        helium%rho_num_act = helium%rho_num_act + itmp
      END IF
    END IF

    ! check for optional estimators based on all paths
    helium%rho_num_min_len_all = 0
    CALL section_vals_val_get( &
      helium%input, &
      "MOTION%PINT%HELIUM%RHO%MIN_CYCLE_LENGTHS_ALL",&
      explicit=explicit)
    IF (explicit) THEN
      NULLIFY(helium%rho_min_len_all_vals)
      CALL section_vals_val_get( &
        helium%input, &
        "MOTION%PINT%HELIUM%RHO%MIN_CYCLE_LENGTHS_ALL", &
        i_vals=helium%rho_min_len_all_vals)
      itmp = SIZE(helium%rho_min_len_all_vals)
      IF (itmp .GT. 0) THEN
        helium%rho_num_min_len_all = itmp
        helium%rho_num_act = helium%rho_num_act + itmp
      END IF
    END IF

    ! get the weighting factor
    CALL section_vals_val_get( &
      helium%input, &
      "MOTION%PINT%HELIUM%AVERAGES%IWEIGHT", &
      i_val=helium%rho_iweight)

    ! allocate and initialize memory for density storage
    itmp = helium%rho_nbin
    jtmp = helium%rho_num_act
    ALLOCATE(helium%rho_inst(jtmp,itmp,itmp,itmp))
    ALLOCATE(helium%rho_accu(jtmp,itmp,itmp,itmp))
    ALLOCATE(helium%rho_rstr(jtmp,itmp,itmp,itmp))
    ALLOCATE(helium%rho_incr(jtmp,helium%atoms,helium%beads))

    helium%rho_inst(:,:,:,:) = 0.0_dp
    helium%rho_accu(:,:,:,:) = 0.0_dp
    helium%rho_rstr(:,:,:,:) = 0.0_dp

    RETURN
  END SUBROUTINE helium_rho_init


! ***************************************************************************
!> \brief  Restore the densities from the input structure.
!> \param helium ...
!> \date   2011-06-22
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE helium_rho_restore( helium )

    TYPE(helium_solvent_type), POINTER       :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_rho_restore', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: stmp1, stmp2
    CHARACTER(len=max_line_length)           :: err_str
    INTEGER                                  :: itmp, msglen
    LOGICAL                                  :: explicit, ltmp
    LOGICAL, DIMENSION(:, :, :, :), POINTER  :: m
    REAL(KIND=dp), DIMENSION(:), POINTER     :: message
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: f

    CALL section_vals_val_get( helium%input, &
      "MOTION%PINT%HELIUM%AVERAGES%RHO", &
      explicit=explicit)
    IF (explicit) THEN
      NULLIFY(message)
      CALL section_vals_val_get( helium%input, &
        "MOTION%PINT%HELIUM%AVERAGES%RHO", &
        r_vals=message)
      msglen = SIZE(message)
      itmp = SIZE(helium%rho_rstr)
      ltmp = ( msglen .EQ. itmp )
      IF ( .NOT. ltmp ) THEN
        stmp1 = ""
        WRITE(stmp1,*) msglen
        stmp2 = ""
        WRITE(stmp2,*) itmp
        err_str = "Size of the S density array in the input (" // &
                  TRIM(ADJUSTL(stmp1)) // &
                  ") .NE. that in the runtime (" // &
                  TRIM(ADJUSTL(stmp2)) // ")."
        CPABORT(err_str)
      END IF
    ELSE
      RETURN
    END IF

    itmp = helium%rho_nbin
    NULLIFY(m,f)
    ALLOCATE(m(helium%rho_num_act,itmp,itmp,itmp))
    ALLOCATE(f(helium%rho_num_act,itmp,itmp,itmp))
    m(:,:,:,:) = .TRUE.
    f(:,:,:,:) = 0.0_dp

    helium%rho_rstr(:,:,:,:) = UNPACK(message(1:msglen), MASK=m, FIELD=f)

    DEALLOCATE(f,m)
    NULLIFY(message)

    RETURN
  END SUBROUTINE helium_rho_restore


! ***************************************************************************
!> \brief Count atoms of different types and store their global indices.
!> \param helium ...
!> \param pint_env ...
!> \author Lukasz Walewski
!> \note  Arrays ALLOCATEd here are (should be) DEALLOCATEd in
!>        helium_release.
! *****************************************************************************
  SUBROUTINE helium_set_solute_indices(helium, pint_env)
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_set_solute_indices', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: cur_ename
    CHARACTER(len=default_path_length)       :: msg_str
    INTEGER                                  :: cur_eid, i, j, mnum, natoms
    INTEGER, DIMENSION(hid_num)              :: el_counts
    LOGICAL                                  :: found
    REAL(KIND=dp)                            :: mass
    TYPE(cp_subsys_type), POINTER            :: my_subsys
    TYPE(f_env_type), POINTER                :: my_f_env
    TYPE(particle_list_type), POINTER        :: my_particles

    ! set up my_particles structure
    NULLIFY(my_f_env, my_subsys, my_particles)
    CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id, &
      f_env=my_f_env)
    CALL force_env_get(force_env=my_f_env%force_env, subsys=my_subsys)
    CALL cp_subsys_get(my_subsys, particles=my_particles)
    CALL f_env_rm_defaults(my_f_env)

    natoms  = helium%solute_atoms
    NULLIFY(helium%solute_element)
    ALLOCATE(helium%solute_element(natoms))

    NULLIFY(helium%ename, helium%eid)
    ALLOCATE(helium%ename(hid_num))
    ALLOCATE(helium%eid(hid_num))

    NULLIFY(helium%solute_i)
    ALLOCATE(helium%solute_i(hid_num))
    DO i = 1, hid_num
!TODO important to nullify here, otherwise helium_destroy_int_arr_ptr
!     in helium_release will fail
!     safer solution: change int_arr_ptr%iap from POINTER to ALLOCATABLE
      NULLIFY(helium%solute_i(i)%iap)
    END DO

    ! find out how many different atomic types are there
    helium%enum = 0
    el_counts(:) = 0
    DO i=1, natoms
      CALL get_atomic_kind( my_particles%els(i)%atomic_kind, mass=mass)
      mnum = NINT(cp_unit_from_cp2k(mass, "amu"))
      SELECT CASE (mnum)
        CASE (35)
          cur_ename = "CL"
          cur_eid = hid_chlorine
        CASE (16)
          cur_ename = "O "
          cur_eid = hid_oxygen
        CASE (1)
          cur_ename = "H "
          cur_eid = hid_hydrogen
        CASE (12)
          cur_ename = "C "
          cur_eid = hid_carbon
        CASE DEFAULT
          cur_ename = "XX"
          cur_eid = -1
          WRITE(msg_str,'(A,I0,A)') &
            "Atomic mass number '", &
            mnum, &
            "' not supported by the HELIUM code."
          CPABORT(msg_str)
      END SELECT
      helium%solute_element(i) = cur_ename
      el_counts(cur_eid) = el_counts(cur_eid) + 1
      ! check if this element symbol is already present in element table
      found = .FALSE.
      DO j=1, helium%enum
        IF ( helium%ename(j) == cur_ename ) THEN
          found = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. found) THEN
        ! increase current number of different elements
        helium%enum = helium%enum + 1
        helium%eid(helium%enum) = cur_eid
        helium%ename(helium%enum) = cur_ename
      END IF
    END DO

    DO i = 1, helium%enum
      ALLOCATE(helium%solute_i(helium%eid(i))%iap(el_counts(helium%eid(i))))
    END DO

    ! store atomic indices of different atomic kinds
    el_counts(:) = 0
    DO i=1, natoms
      SELECT CASE (helium%solute_element(i))
      CASE ("CL")
        el_counts(hid_chlorine) = el_counts(hid_chlorine) + 1
        helium%solute_i(hid_chlorine)%iap(el_counts(hid_chlorine)) = i
      CASE ("O ")
        el_counts(hid_oxygen) = el_counts(hid_oxygen) + 1
        helium%solute_i(hid_oxygen)%iap(el_counts(hid_oxygen)) = i
      CASE ("H ")
        el_counts(hid_hydrogen) = el_counts(hid_hydrogen) + 1
        helium%solute_i(hid_hydrogen)%iap(el_counts(hid_hydrogen)) = i
      CASE ("C ")
        el_counts(hid_carbon) = el_counts(hid_carbon) + 1
        helium%solute_i(hid_carbon)%iap(el_counts(hid_carbon)) = i
      CASE DEFAULT
          WRITE(msg_str,'(A)') &
            "Should never get here, check what happened!"
          CPABORT(msg_str)
      END SELECT
    END DO

    RETURN
  END SUBROUTINE helium_set_solute_indices


! ***************************************************************************
!> \brief Sets helium%solute_cell based on the solute's force_env.
!> \param helium ...
!> \param pint_env ...
!> \author Lukasz Walewski
!> \note  The simulation cell for the solvated molecule is taken from force_env
!>        which should assure that we get proper cell dimensions regardless of
!>        the method used for the solute (QS, FIST). Helium solvent needs the
!>        solute's cell dimensions to calculate the solute-solvent distances
!>        correctly.
!> \note  At the moment only orthorhombic cells are supported.
! *****************************************************************************
  SUBROUTINE helium_set_solute_cell(helium, pint_env)
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_set_solute_cell', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: my_orthorhombic
    TYPE(cell_type), POINTER                 :: my_cell
    TYPE(f_env_type), POINTER                :: my_f_env

    ! get the cell structure from pint_env
    NULLIFY(my_f_env, my_cell)
    CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id, &
      f_env=my_f_env)
    CALL force_env_get(force_env=my_f_env%force_env, cell=my_cell)
    CALL f_env_rm_defaults(my_f_env)

    CALL get_cell(my_cell, orthorhombic=my_orthorhombic)
    IF (.NOT. my_orthorhombic) THEN
      CPABORT("Helium solvent not implemented for non-orthorhombic cells.")
    ELSE
      helium%solute_cell => my_cell
    END IF

    RETURN
  END SUBROUTINE helium_set_solute_cell


END MODULE helium_methods
