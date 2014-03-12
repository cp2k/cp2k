!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Methods dealing with helium_solvent_type
!> \author Lukasz Walewski
!> \date   2009-06-10
! *****************************************************************************
MODULE helium_methods

  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE f77_interface,                   ONLY: f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get
  USE helium_common,                   ONLY: helium_eval_expansion,&
                                             helium_pbc
  USE helium_io,                       ONLY: helium_write_line,&
                                             helium_write_setup
  USE helium_sampling,                 ONLY: helium_sample
  USE helium_types,                    ONLY: he_mass,&
                                             helium_solvent_type
  USE input_constants,                 ONLY: helium_cell_shape_cube,&
                                             helium_cell_shape_none,&
                                             helium_cell_shape_octahedron,&
                                             silent_print_level
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
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
  USE physcon,                         ONLY: angstrom,&
                                             kelvin,&
                                             massunit
  USE pint_public,                     ONLY: pint_com_pos
  USE pint_types,                      ONLY: pint_env_type
  USE splines_methods,                 ONLY: init_spline,&
                                             init_splinexy
  USE splines_types,                   ONLY: spline_data_create,&
                                             spline_data_release,&
                                             spline_data_retain
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_methods'
  INTEGER, SAVE, PRIVATE :: last_helium_id=0

  PUBLIC :: helium_create
  PUBLIC :: helium_init
  PUBLIC :: helium_release
  PUBLIC :: helium_total_action

  CONTAINS

  ! ***************************************************************************
  !> \brief  Data-structure that holds all needed information about
  !>         (superfluid) helium solvent
  !> \author hforbert
  ! ***************************************************************************
  SUBROUTINE helium_create( helium, input, solute, error )
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(section_vals_type), POINTER         :: input
    TYPE(pint_env_type), OPTIONAL, POINTER   :: solute
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: potential_file_name
    CHARACTER(len=default_string_length)     :: msg_str, stmp
    INTEGER                                  :: handle, i, input_unit, isize, &
                                                itmp, j, jtmp, nlines, ntab, &
                                                stat
    LOGICAL                                  :: cell_shape_supported, &
                                                expl_cell, expl_dens, &
                                                expl_nats, explicit, failure, &
                                                ltmp
    REAL(KIND=dp)                            :: cgeof, dx, rtmp, tcheck, x1
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pot_transfer
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: helium_section

    CALL timeset(routineN,handle)

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(helium),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)
    CPPrecondition(input%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(helium_section)
       helium_section => section_vals_get_subs_vals(input, &
                                              "MOTION%PINT%HELIUM",error=error)
       CALL section_vals_get(helium_section,explicit=explicit,error=error)
       CPPostcondition(explicit,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       ALLOCATE(helium, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
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
                 helium%rdf_avrg,&
                 helium%rdf_inst,&
                 helium%plength_avrg, &
                 helium%plength_inst, &
                 helium%atom_plength &
                )

       helium%ref_count=1
       last_helium_id=last_helium_id+1
       helium%id_nr=last_helium_id
       helium%input => input
       helium%accepts = 0
       helium%relrot = 0

       NULLIFY(logger)
       logger => cp_error_get_logger(error)

       ! get number of environments in the restart file (if present)
       CALL section_vals_val_get(helium_section,"NUM_ENV",&
            explicit=explicit, error=error)
       IF ( explicit ) THEN
         CALL section_vals_val_get(helium_section,"NUM_ENV",&
              i_val=itmp, error=error)
         CPPostcondition(itmp>=0,cp_failure_level,routineP,error,failure)
         helium%num_env_restart = itmp
         ELSE
         helium%num_env_restart = -1
       END IF

       ! set current number of environments
       helium%num_env = logger%para_env%num_pe
       CALL section_vals_val_set(helium%input,&
            "MOTION%PINT%HELIUM%NUM_ENV",&
            i_val=helium%num_env, error=error)

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
         CALL helium_write_line(msg_str,error)
         CALL section_vals_val_get(helium_section,"DROP_UNUSED_ENVS",&
              l_val=ltmp, error=error)
         IF (ltmp) THEN
           msg_str = "DROP_UNUSED_ENVS set - proceeding anyways."
           CALL helium_write_line(msg_str,error)
         ELSE
           CALL cp_assert( .FALSE., cp_failure_level, &
                cp_assertion_failed, routineP, msg_str )
         END IF
       END IF

       CALL section_vals_val_get(helium_section,"NBEADS",&
              i_val=helium%beads, error=error)
       CALL section_vals_val_get(helium_section,"INOROT",&
              i_val=helium%iter_norot, error=error)
       CALL section_vals_val_get(helium_section,"IROT",&
              i_val=helium%iter_rot, error=error)

       ! get number of steps and current step number from PINT
       CALL section_vals_val_get(input,"MOTION%PINT%ITERATION",&
            i_val=itmp, error=error)
       helium%first_step = itmp
       CALL section_vals_val_get(input,"MOTION%PINT%MAX_STEP",&
            explicit=explicit, error=error)
       IF ( explicit ) THEN
         CALL section_vals_val_get(input,"MOTION%PINT%MAX_STEP",&
              i_val=itmp, error=error)
         helium%last_step = itmp
         helium%num_steps = helium%last_step - helium%first_step
       ELSE
         CALL section_vals_val_get(input,"MOTION%PINT%NUM_STEPS",&
              i_val=itmp, error=error)
         helium%num_steps = itmp
         helium%last_step = helium%first_step + helium%num_steps
       END IF

       ! If we should apply periodicity, check wheather we support the cell
       ! shape or not, and refuse to proceed if not. Otherwise just go on.
       CALL section_vals_val_get(helium_section,"PERIODIC",&
            l_val=helium%periodic, error=error)
       IF ( helium%periodic ) THEN
         CALL section_vals_val_get(helium_section,"CELL_SHAPE",&
              i_val=helium%cell_shape, error=error)
         cell_shape_supported = .FALSE.
         IF ( helium%cell_shape .EQ. helium_cell_shape_cube ) THEN
           cell_shape_supported = .TRUE.
         END IF
         IF ( helium%cell_shape .EQ. helium_cell_shape_octahedron ) THEN
           cell_shape_supported = .TRUE.
         END IF
         IF ( .NOT. cell_shape_supported ) THEN
           stmp = ""
           WRITE(stmp,*) helium%cell_shape
           msg_str = "PBC cell shape " // &
                     TRIM(ADJUSTL(stmp)) // &
                     " not supported for HELIUM"
           CALL cp_assert( .FALSE., cp_failure_level, &
             cp_assertion_failed, routineP, msg_str )
         END IF
       ELSE
         helium%cell_shape = helium_cell_shape_none
       END IF
       ! Further PBC code depends on the correct value of helium%cell_shape.

       ! Set density Rho, number of atoms N and volume V ( Rho = N / V ).
       ! Allow only 2 out of 3 values to be defined at the same time, calculate
       ! the third.
       ! Note, that DENSITY and NATOMS keywords have default values, while
       ! CELL_SIZE does not. Thus if CELL_SIZE is given explicitly then one and
       ! only one of the two remaining options must be give explicitly as well.
       ! If CELL_SIZE is not given explicitly then all four combinations of the
       ! two other options are valid.
       CALL section_vals_val_get(helium_section,"DENSITY",&
            explicit=expl_dens, r_val=helium%density, error=error)
       CALL section_vals_val_get(helium_section,"NATOMS",&
            explicit=expl_nats, i_val=helium%atoms, error=error)
       CALL section_vals_val_get(helium_section,"CELL_SIZE",&
            explicit=expl_cell, error=error)
       cgeof = 1.0_dp
       IF ( helium%cell_shape .EQ. helium_cell_shape_octahedron ) cgeof = 2.0_dp
       rtmp = ( cgeof * helium%atoms / helium%density )**(1.0_dp/3.0_dp)
       IF ( .NOT. expl_cell ) THEN
           helium%cell_size = rtmp
       ELSE
           CALL section_vals_val_get(helium_section,"CELL_SIZE",&
                r_val=helium%cell_size, error=error)
           ! only more work if not all three values are consistent:
           IF ( ABS(helium%cell_size-rtmp) .GT. 100.0_dp*EPSILON(0.0_dp)* &
                                          (ABS(helium%cell_size)+rtmp) ) THEN
              IF ( expl_dens .AND. expl_nats ) THEN
                 msg_str = "DENSITY, NATOMS and CELL_SIZE options "//&
                           "contradict each other"
                 CALL cp_assert( .FALSE., cp_failure_level, &
                      cp_assertion_failed, routineP, msg_str )
              END IF
              !ok we have enough freedom to resolve the conflict:
              IF ( .NOT. expl_dens ) THEN
                 helium%density = cgeof*helium%atoms / helium%cell_size**3.0_dp
                 IF ( .NOT. expl_nats ) THEN
                    msg_str = "Warning: CELL_SIZE defined but neither "//&
                              "NATOMS nor DENSITY given, using default NATOMS."
                    CALL helium_write_line(msg_str, error)
                 END IF
              ELSE ! ( expl_dens .AND. .NOT. expl_nats )
                 ! calculate the nearest number of atoms for given conditions
                 helium%atoms = ANINT(helium%density * &
                                helium%cell_size**3.0_dp / cgeof)
                 ! adjust cell size to maintain correct density
                 ! (should be a small correction)
                 rtmp = ( cgeof * helium%atoms / helium%density &
                        )**(1.0_dp/3.0_dp)
                 IF ( ABS(helium%cell_size-rtmp) .GT. 100.0_dp*EPSILON(0.0_dp)&
                      * ( ABS(helium%cell_size)+rtmp ) ) THEN
                    msg_str = "Warning: Adjusting actual cell size "//&
                              "to maintain correct density."
                    CALL helium_write_line(msg_str, error)
                    helium%cell_size = rtmp
                 END IF
              END IF
           END IF
       END IF
       helium%cell_size_inv = 1.0_dp / helium%cell_size
       ! From now on helium%density, helium%atoms and helium%cell_size are
       ! correctly defined.

       ! set the M matrix for winding number calculations
       IF ( helium%cell_shape .EQ. helium_cell_shape_octahedron ) THEN
         helium%cell_m_inv(1,1) =  1.0_dp / helium%cell_size;
         helium%cell_m_inv(2,1) =  0.0_dp;
         helium%cell_m_inv(3,1) =  0.0_dp;
         helium%cell_m_inv(1,2) =  0.0_dp;
         helium%cell_m_inv(2,2) =  1.0_dp / helium%cell_size;
         helium%cell_m_inv(3,2) =  0.0_dp;
         helium%cell_m_inv(1,3) = -1.0_dp / helium%cell_size;
         helium%cell_m_inv(2,3) = -1.0_dp / helium%cell_size;
         helium%cell_m_inv(3,3) =  2.0_dp / helium%cell_size;
       ELSE
         helium%cell_m_inv(1,1) =  1.0_dp / helium%cell_size;
         helium%cell_m_inv(2,1) =  0.0_dp;
         helium%cell_m_inv(3,1) =  0.0_dp;
         helium%cell_m_inv(1,2) =  0.0_dp;
         helium%cell_m_inv(2,2) =  1.0_dp / helium%cell_size;
         helium%cell_m_inv(3,2) =  0.0_dp;
         helium%cell_m_inv(1,3) =  0.0_dp;
         helium%cell_m_inv(2,3) =  0.0_dp;
         helium%cell_m_inv(3,3) =  1.0_dp / helium%cell_size;
       END IF

       ! check value of maxcycle
       CALL section_vals_val_get(helium_section,"MAX_PERM_CYCLE",&
              i_val=helium%maxcycle, error=error)
       i = helium%maxcycle
       CPPostcondition(i>=0,cp_failure_level,routineP,error,failure)
       i = helium%atoms - helium%maxcycle
       CPPostcondition(i>=0,cp_failure_level,routineP,error,failure)

       ! set m-distribution parameters
       CALL section_vals_val_get(helium_section,"M-SAMPLING%M-VALUE",&
            i_val=i, error=error)
       CPPostcondition(i>=1,cp_failure_level,routineP,error,failure)
       CPPostcondition(i<=helium%maxcycle,cp_failure_level,routineP,error,failure)
       helium%m_value = i
       CALL section_vals_val_get(helium_section,"M-SAMPLING%M-RATIO",&
            r_val=rtmp, error=error)
       CPPostcondition(rtmp>0.0_dp,cp_failure_level,routineP,error,failure)
       CPPostcondition(rtmp<=1.0_dp,cp_failure_level,routineP,error,failure)
       helium%m_ratio = rtmp

       CALL section_vals_val_get(helium_section,"BISECTION",&
              i_val=helium%bisection, error=error)
       ! precheck bisection value (not all invalids are filtered out here yet)
       i = helium%bisection
       CPPostcondition(i>1,cp_failure_level,routineP,error,failure)
       i = helium%beads - helium%bisection
       CPPostcondition(i>0,cp_failure_level,routineP,error,failure)
       !
       itmp = helium%bisection
       rtmp = 2**(ANINT(LOG(REAL(itmp))/LOG(2.0_dp)))
       tcheck=ABS(REAL(itmp,KIND=dp)-rtmp)
       msg_str = "BISECTION should be integer power of 2."
       CALL cp_assert(tcheck<100.0_dp*EPSILON(0.0_dp),cp_failure_level,&
            cp_assertion_failed, routineP, msg_str)
       helium%bisctlog2 = ANINT(LOG(REAL(itmp))/LOG(2.0_dp))

       ! hard coded He4 directly (mass so i get my original hb2m value)
       ! He4 mass defined as a constant in helium_types.F now [lwalewski]
       helium%hb2m = 1.0_dp/(he_mass*massunit)
       helium%pweight = 0.0_dp

       ! get the RDF parameters
       CALL section_vals_val_get(helium_section,"RDF%MAXR",&
            explicit=explicit, error=error)
       IF (explicit) THEN
         CALL section_vals_val_get(helium_section,"RDF%MAXR",&
              r_val=helium%rdf_maxr, error=error)
       ELSE
         helium%rdf_maxr = helium%cell_size
       END IF
       CALL section_vals_val_get(helium_section,"RDF%NBIN",&
            i_val=helium%rdf_nbin, error=error)
       helium%rdf_delr = helium%rdf_maxr / REAL(helium%rdf_nbin,KIND=dp)

       IF (logger%para_env%ionode) THEN
         CALL section_vals_val_get(helium_section,"POTENTIAL_FILE_NAME",&
             c_val=potential_file_name, error=error)
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

       isize = helium%pdx+1
       ALLOCATE(helium%uij(isize,isize),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%eij(isize,isize),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, isize
          DO j = 1, i
             CALL spline_data_create(helium%uij(i,j)%spline_data, &
                                          error=error)
             CALL init_splinexy(helium%uij(i,j)%spline_data,nlines)
             helium%uij(i,j)%spline_data%x1 = x1
             CALL spline_data_create(helium%eij(i,j)%spline_data, &
                                          error=error)
             CALL init_splinexy(helium%eij(i,j)%spline_data,nlines)
             helium%eij(i,j)%spline_data%x1 = x1
          END DO
       END DO
       DO i = 1, isize-1
          DO j = i+1, isize
             helium%uij(i,j) = helium%uij(j,i)
             CALL spline_data_retain(helium%uij(i,j)%spline_data, &
                                          error=error)
             helium%eij(i,j) = helium%eij(j,i)
             CALL spline_data_retain(helium%eij(i,j)%spline_data, &
                                          error=error)
          END DO
       END DO

       isize = (helium%pdx+1)*(helium%pdx+2)
       ALLOCATE(pot_transfer(nlines,isize), STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
       DEALLOCATE(pot_transfer, STAT=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

       ! ALLOCATE helium-related arrays
       i = helium%atoms
       j = helium%beads
       ALLOCATE(helium%pos(3,i,j),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%work(3,i,j),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%ptable(helium%maxcycle+1),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%permutation(i),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%iperm(i),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%tmatrix(i,i),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%nmatrix(i,2*i),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%pmatrix(i,i),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%ipmatrix(i,i),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       itmp = helium%bisctlog2 + 2
       ALLOCATE(helium%num_accepted(itmp,helium%maxcycle),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%rdf_avrg(helium%rdf_nbin),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%rdf_inst(helium%rdf_nbin),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%plength_avrg(helium%atoms),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%plength_inst(helium%atoms),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%atom_plength(helium%atoms),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ! density related data
       helium%rho_num = 1
       CALL section_vals_val_get(helium_section, &
            "RHO%_SECTION_PARAMETERS_", &
            l_val=helium%rho_present, error=error)
       IF ( helium%rho_present ) THEN
         helium%rho_maxr = helium%cell_size
         CALL section_vals_val_get(helium%input, &
              "MOTION%PINT%HELIUM%RHO%IWEIGHT",&
              i_val=helium%rho_iweight, error=error)
         CALL section_vals_val_get(helium_section,"RHO%NBIN",&
              i_val=helium%rho_nbin, error=error)
         helium%rho_delr = helium%rho_maxr / REAL(helium%rho_nbin,KIND=dp)
         helium%rho_minb = INT(-helium%cell_size / 2.0_dp / helium%rho_delr )
         itmp = helium%rho_nbin
         jtmp = helium%rho_num
         ALLOCATE(helium%rho_avrg(jtmp,itmp,itmp,itmp),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(helium%rho_inst(jtmp,itmp,itmp,itmp),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         helium%rho_avrg(:,:,:,:) = 0.0_dp
         helium%rho_inst(:,:,:,:) = 0.0_dp

         ! make sure that the density gets actually printed,
         ! as we will make so much effort to calculated it
         CALL section_vals_val_set(helium%input, &
              "MOTION%PINT%HELIUM%PRINT%RHO%_SECTION_PARAMETERS_",&
              i_val=silent_print_level-1, error=error)
       END IF

       CALL section_vals_val_get(helium%input, &
            "EXT_RESTART%RESTART_HELIUM_DENSITIES", &
            l_val=helium%rho_restart,error=error)
       IF ( helium%rho_restart ) THEN
         IF ( helium%rho_present) THEN
           CALL helium_densities_restore( helium, error )
         ELSE
           msg_str = "Inconsistent input: densities not calculated but restarted!"
           CALL cp_assert(.FALSE.,cp_failure_level,&
                cp_assertion_failed, routineP, msg_str)
         END IF
       END IF

       ! RNG state create & init
       CALL helium_rng_state_init( helium, error )

       ! check if solute is present in our simulation
       helium%solute_present = .FALSE.
       IF (PRESENT(solute)) THEN
         IF (ASSOCIATED(solute)) THEN
           helium%solute_present = .TRUE.
         END IF
       END IF

       ! fill in the solute-related data structures
       helium%e_corr = 0.0_dp
       IF (helium%solute_present) THEN
         helium%solute_atoms = solute%ndim / 3
         helium%solute_beads = solute%p
         helium%bead_ratio = helium%beads / helium%solute_beads

         ! check if bead numbers are commensurate:
         i = helium%bead_ratio*helium%solute_beads - helium%beads
!TODO Adjust helium bead number if not comm. and if coords not given expl.
         CPPostcondition(i==0,cp_failure_level,routineP,error,failure)

         ! check if tau, temperature and bead number are consistent:
         tcheck=ABS( (helium%tau*helium%beads-solute%beta) / solute%beta )
         msg_str = "Tau, temperature and bead number are inconsistent."
         CALL cp_assert(tcheck<1.0e-14_dp,cp_failure_level,&
           cp_assertion_failed, routineP, msg_str)

         CALL helium_set_solute_indices(helium,solute)
         CALL helium_set_solute_cell(helium,solute)
       ELSE
         helium%solute_atoms = 0
         helium%solute_beads = 0
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

       ! ALLOCATE solute-related arrays
       ALLOCATE(helium%force_avrg(helium%solute_beads,&
         helium%solute_atoms*3),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%force_inst(helium%solute_beads,&
         helium%solute_atoms*3),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ! ALLOCATE temporary arrays
       ALLOCATE(helium%itmp_atoms_1d(helium%atoms),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%ltmp_atoms_1d(helium%atoms),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%itmp_atoms_np_1d(helium%atoms*helium%num_env),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%rtmp_3_np_1d(3*helium%num_env),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%rtmp_p_ndim_1d(helium%solute_beads*&
         helium%solute_atoms*3),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%rtmp_p_ndim_np_1d(helium%solute_beads*&
         helium%solute_atoms*3*helium%num_env),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%rtmp_3_atoms_beads_1d(3*helium%atoms*&
         helium%beads),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%rtmp_3_atoms_beads_np_1d(3*helium%atoms*&
         helium%beads*helium%num_env),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%rtmp_p_ndim_2d(helium%solute_beads,&
         helium%solute_atoms*3),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(helium%ltmp_3_atoms_beads_3d(3,helium%atoms,&
         helium%beads),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       CALL helium_write_setup(helium,error)

    END IF

    CALL timestop(handle)

    RETURN
  END SUBROUTINE helium_create

  ! ***************************************************************************
  !> \brief  Releases helium_solvent_type
  !> \author hforbert
  ! ***************************************************************************
  SUBROUTINE helium_release(helium,error)
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(helium)) THEN
      CPPrecondition(helium%ref_count>0,cp_failure_level,routineP,error,failure)
      helium%ref_count=helium%ref_count-1
      IF (helium%ref_count<1) THEN

        ! DEALLOCATE temporary arrays
        DEALLOCATE ( &
          helium%ltmp_3_atoms_beads_3d, &
          helium%rtmp_p_ndim_2d, &
          helium%rtmp_3_atoms_beads_np_1d, &
          helium%rtmp_3_atoms_beads_1d, &
          helium%rtmp_p_ndim_np_1d, &
          helium%rtmp_p_ndim_1d, &
          helium%rtmp_3_np_1d, &
          helium%itmp_atoms_np_1d, &
          helium%ltmp_atoms_1d, &
          helium%itmp_atoms_1d, &
          STAT=stat )
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        NULLIFY ( &
          helium%ltmp_3_atoms_beads_3d, &
          helium%rtmp_p_ndim_2d, &
          helium%rtmp_3_atoms_beads_np_1d, &
          helium%rtmp_3_atoms_beads_1d, &
          helium%rtmp_p_ndim_np_1d, &
          helium%rtmp_p_ndim_1d, &
          helium%rtmp_3_np_1d, &
          helium%itmp_atoms_np_1d, &
          helium%ltmp_atoms_1d, &
          helium%itmp_atoms_1d &
        )

        ! DEALLOCATE solute-related arrays
        DEALLOCATE ( &
          helium%force_inst, &
          helium%force_avrg, &
          STAT=stat )
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        NULLIFY( &
          helium%force_inst, &
          helium%force_avrg &
        )

        IF ( helium%rho_present ) THEN
          IF ( helium%rho_restart ) THEN
            DEALLOCATE( helium%rho_rstr, STAT=stat )
            CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
            NULLIFY( helium%rho_rstr )
          END IF
          DEALLOCATE( &
            helium%rho_inst, &
            helium%rho_avrg, &
            STAT=stat &
          )
          NULLIFY(&
            helium%rho_inst, &
            helium%rho_avrg &
          )
        END IF

        ! DEALLOCATE helium-related arrays
        DEALLOCATE ( &
          helium%atom_plength, &
          helium%plength_inst, &
          helium%plength_avrg, &
          helium%rdf_inst, &
          helium%rdf_avrg, &
          helium%num_accepted, &
          helium%ipmatrix, &
          helium%pmatrix, &
          helium%nmatrix, &
          helium%tmatrix, &
          helium%iperm, &
          helium%permutation, &
          helium%ptable, &
          helium%work, &
          helium%pos, &
          STAT=stat )
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        NULLIFY( &
          helium%atom_plength, &
          helium%plength_inst, &
          helium%plength_avrg, &
          helium%rdf_inst, &
          helium%rdf_avrg, &
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
                CALL spline_data_release(helium%eij(i,j)%spline_data,&
                                                                 error=error)
                CALL spline_data_release(helium%uij(i,j)%spline_data,&
                                                                 error=error)
                !TODO: shouldn't that be done in spline_data_release??
                NULLIFY(helium%eij(i,j)%spline_data, &
                        helium%uij(i,j)%spline_data)
             END DO
          END DO

          DEALLOCATE( helium%eij, STAT=stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          NULLIFY(helium%eij)

          DEALLOCATE( helium%uij, STAT=stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          NULLIFY(helium%uij)

          CALL delete_rng_stream(helium%rng_stream_uniform,error=error)
          CALL delete_rng_stream(helium%rng_stream_gaussian,error=error)

         ! deallocate solute-related arrays
          IF (helium%solute_present) THEN
            DEALLOCATE(helium%solute_element, &
                       helium%solute_number, &
                       helium%solute_index, &
                       STAT=stat)
            CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
            NULLIFY(helium%solute_element, &
                    helium%solute_number, &
                    helium%solute_index)
          END IF

          DEALLOCATE( helium, STAT=stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

       END IF
    END IF
    RETURN
  END SUBROUTINE helium_release

  ! ***************************************************************************
  !> \brief  Retains helium_solvent_type
  !> \author hforbert
  ! ***************************************************************************
  SUBROUTINE helium_retain(helium,error)
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(helium%ref_count>0,cp_failure_level,routineP,error,failure)
       helium%ref_count=helium%ref_count+1
    END IF
    RETURN
  END SUBROUTINE helium_retain

  ! ***************************************************************************
  !> \brief  Initialize helium data structures.
  !> \author hforbert
  !> \descr  Initializes helium coordinates either as random positions or from
  !>         HELIUM%COORD section if it's present in the input file.
  !>         Initializes helium permutation state as identity permutation or
  !>         from HELIUM%PERM section if it's present in the input file.
  !> \par    History
  !>         removed refereces to pint_env_type data structure [lwalewski]
  !>         2009-11-10 init/restore coords, perm, RNG and forces [lwalewski]
  ! ***************************************************************************
  SUBROUTINE helium_init( helium, pint_env, error )

    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: coords_presampled, explicit, &
                                                explicit_presample, failure, &
                                                presample
    TYPE(section_vals_type), POINTER         :: helium_section, sec

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(helium%input),cp_failure_level,routineP,error,failure)

    NULLIFY(helium_section)
    helium_section => section_vals_get_subs_vals(helium%input, &
      "MOTION%PINT%HELIUM", error=error)

    ! set the origin for rdf and density calculations
    IF (helium%solute_present) THEN
      helium%origin = pint_com_pos(pint_env,error)
    ELSE
      helium%origin = (/0.0_dp, 0.0_dp, 0.0_dp/)
    END IF

    ! restore RNG state
    NULLIFY(sec)
    sec => section_vals_get_subs_vals(helium_section,"RNG_STATE",error=error)
    CALL section_vals_get(sec,explicit=explicit,error=error)
    IF ( explicit ) THEN
      CALL helium_rng_state_restore( helium, error )
    ELSE
      CALL helium_write_line("RNG state initialized as new.",error)
    END IF

    ! init/restore permutation state
    NULLIFY(sec)
    sec => section_vals_get_subs_vals(helium_section,"PERM",error=error)
    CALL section_vals_get(sec,explicit=explicit,error=error)
    IF ( explicit ) THEN
      CALL helium_perm_restore( helium, error )
    ELSE
      CALL helium_perm_init( helium, error )
      CALL helium_write_line("Permutation state initialized as identity.",error)
    END IF

    ! init/restore coordinates
    NULLIFY(sec)
    sec => section_vals_get_subs_vals(helium_section,"COORD",error=error)
    CALL section_vals_get(sec,explicit=explicit,error=error)

    CALL section_vals_val_get(helium_section,"PRESAMPLE",&
         l_val=presample, error=error)
    CALL section_vals_val_get(helium_section,"PRESAMPLE",&
         explicit=explicit_presample, error=error)

    coords_presampled = .FALSE.
    IF ( explicit ) THEN
      CALL helium_coord_restore( helium, error )
      IF ( presample ) THEN
        ! Default value of PINT%HELIUM%PRESAMPLE is .FALSE. so if the actual
        ! value is .TRUE. it must have been supplied explicitly.
        ! This might be used to force presampling even though the cooridnates
        ! are given in the restart file.
        CALL helium_sample( helium, pint_env, error )
        helium%force_avrg(:,:) = 0.0_dp
        helium%energy_avrg(:) = 0.0_dp
        helium%wnumber_avrg(:) = 0.0_dp
        helium%sdensity_avrg = 0.0_dp
        helium%plength_avrg(:) = 0.0_dp
        helium%rdf_avrg(:) = 0.0_dp
        helium%num_accepted(:,:) = 0.0_dp
        coords_presampled = .TRUE.
        CALL helium_write_line("Bead coordinates pre-sampled.",error)
      END IF
    ELSE
      CALL helium_coord_init( helium, pint_env, error )
      CALL helium_write_line("Bead coordinates initialized as random.",error)
      IF ( (explicit_presample .AND. presample) .OR. &
           (.NOT. explicit_presample) ) THEN
        ! Perform initial MC sampling to get rid of the overlaps.
        ! This IF catches two cases: a) PRESAMPLE given explicitly and it's
        ! value set to .TRUE. and b) PRESAMPLE not given explicitly (then
        ! it's value is set to .FALSE. by default. It does not catch the
        ! case when PRESAMPLE is given explicitly and it's value is set to
        ! .FALSE. - this might be used to skip presampling even though the
        ! coordinates are not given in the restart file.
        CALL helium_sample( helium, pint_env, error )
        helium%force_avrg(:,:) = 0.0_dp
        helium%energy_avrg(:) = 0.0_dp
        helium%wnumber_avrg(:) = 0.0_dp
        helium%sdensity_avrg = 0.0_dp
        helium%plength_avrg(:) = 0.0_dp
        helium%rdf_avrg(:) = 0.0_dp
        helium%num_accepted(:,:) = 0.0_dp
        coords_presampled = .TRUE.
        CALL helium_write_line("Bead coordinates pre-sampled.",error)
      END IF
    END IF

    IF ( helium%solute_present ) THEN
      ! restore helium forces
      NULLIFY(sec)
      sec => section_vals_get_subs_vals(helium_section,"FORCE",error=error)
      CALL section_vals_get(sec,explicit=explicit,error=error)
      IF ( explicit ) THEN
        IF ( .NOT. coords_presampled ) THEN
          CALL helium_force_restore( helium, error )
        END IF
      ELSE
        IF ( .NOT. coords_presampled ) THEN
          CALL helium_force_init( helium, pint_env, error )
          CALL helium_write_line("Forces on the solute initialized as zero.",error)
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
  !> \author Lukasz Walewski
  !> \date   2009-11-09
  ! ***************************************************************************
  SUBROUTINE helium_coord_init( helium, pint_env, error )
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_coord_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ib, ic
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: r1, r2

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    DO ia = 1, helium%atoms
      DO ic = 1, 3
        r1 = next_random_number(helium%rng_stream_uniform,error=error)
        r1 = r1 * helium%cell_size
        DO ib = 1, helium%beads
!TODO use thermal gaussian instead
          r2 = next_random_number(helium%rng_stream_uniform,error=error)
          helium%pos(ic,ia,ib) = r1+0.1_dp*r2
        END DO
      END DO
      DO ib = 1, helium%beads
        CALL helium_pbc( helium, helium%pos(:,ia,ib) )
      END DO
    END DO

    RETURN
  END SUBROUTINE helium_coord_init

  ! ***************************************************************************
  !> \brief  Restore coordinates from the input structure.
  !> \author Lukasz Walewski
  !> \date   2009-11-09
  !> \par    History
  !>         2010-07-22 accomodate additional cpus in the runtime wrt the
  !>                    restart [lwalewski]
  ! ***************************************************************************
  SUBROUTINE helium_coord_restore( helium, error )
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_coord_restore', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: err_str, stmp
    INTEGER                                  :: actlen, msglen, &
                                                num_env_restart, offset, &
                                                status
    LOGICAL                                  :: failure
    LOGICAL, DIMENSION(:, :, :), POINTER     :: m
    REAL(kind=dp), DIMENSION(:), POINTER     :: message
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: f
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    ! assign the pointer to the memory location of the input structure, where
    ! the coordinates are stored
    NULLIFY(message)
    CALL section_vals_val_get( helium%input, &
         "MOTION%PINT%HELIUM%COORD%_DEFAULT_KEYWORD_", &
         r_vals=message, error=error)

    ! check that the number of values in the input match the current runtime
    actlen = SIZE(message)
    num_env_restart = actlen / helium%atoms / helium%beads / 3

    ! distribute coordinates over processors (no message passing)
    msglen = helium%atoms * helium%beads * 3
    offset = msglen * MOD( logger%para_env%mepos, num_env_restart )
    NULLIFY(m,f)
    ALLOCATE(m(3,helium%atoms,helium%beads),STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(f(3,helium%atoms,helium%beads),STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)
    m(:,:,:) = .TRUE.
    f(:,:,:) = 0.0_dp
    helium%pos(:,:,:) = UNPACK(message(offset+1:offset+msglen), MASK=m, FIELD=f )
    DEALLOCATE(f,m,STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)

    IF ( num_env_restart .NE. logger%para_env%num_pe ) THEN
      err_str = "Reading bead coordinates from the input file."
      CALL helium_write_line(err_str,error)
      err_str = "Number of environments in the restart...: '"
      stmp = ""
      WRITE(stmp,*) num_env_restart
      err_str = TRIM(ADJUSTL(err_str)) // TRIM(ADJUSTL(stmp)) // "'."
      CALL helium_write_line(err_str,error)
      err_str = "Number of current run time environments.: '"
      stmp = ""
      WRITE(stmp,*) logger%para_env%num_pe
      err_str = TRIM(ADJUSTL(err_str)) // TRIM(ADJUSTL(stmp)) // "'."
      CALL helium_write_line(err_str,error)
      IF ( num_env_restart .LT. logger%para_env%num_pe ) THEN
        err_str = "Replicated bead coordinates from the restarted environments."
        CALL helium_write_line(err_str,error)
      END IF
      IF ( num_env_restart .GT. logger%para_env%num_pe ) THEN
        err_str = "Dropped bead coordinates from some restarted environments."
        CALL helium_write_line(err_str,error)
      END IF
      err_str = "Done."
      CALL helium_write_line(err_str,error)
    ELSE
      CALL helium_write_line("Bead coordinates read from the input file.",error)
    END IF

    NULLIFY(message)

    RETURN
  END SUBROUTINE helium_coord_restore

  ! ***************************************************************************
  !> \brief  Initialize forces exerted on the solute.
  !> \author Lukasz Walewski
  !> \date   2009-11-10
  !> \descr  The forces are calculated based on both the helium and the solute
  !>         positions, hence this function should be called AFTER
  !>         helium_coord_init/restore.
  ! ***************************************************************************
  SUBROUTINE helium_force_init( helium, pint_env, error )
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_force_init', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    IF ( helium%solute_present ) THEN
!TODO initial forces are set to 0 due to possible overlapping atoms
!      CALL helium_solute_e_f(pint_env, helium, rtmp)
!    ELSE
      helium%force_avrg(:,:) = 0.0_dp
      helium%force_inst(:,:) = 0.0_dp
    END IF

    RETURN
  END SUBROUTINE helium_force_init

  ! ***************************************************************************
  !> \brief  Restore forces from the input structure to the runtime environment.
  !> \author Lukasz Walewski
  !> \date   2009-11-10
  ! ***************************************************************************
  SUBROUTINE helium_force_restore( helium, error )
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_force_restore', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: err_str, stmp
    INTEGER                                  :: actlen, msglen, &
                                                num_env_restart, status
    LOGICAL                                  :: failure
    LOGICAL, DIMENSION(:, :), POINTER        :: m
    REAL(kind=dp), DIMENSION(:), POINTER     :: message
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: f

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    ! assign the pointer to the memory location of the input structure, where
    ! the forces are stored
    NULLIFY(message)
    CALL section_vals_val_get( helium%input, &
         "MOTION%PINT%HELIUM%FORCE%_DEFAULT_KEYWORD_", &
         r_vals=message, error=error)

    ! check the number of environments presumably stored in the restart
    actlen = SIZE(message)
    num_env_restart = actlen / helium%solute_atoms / helium%solute_beads / 3

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
    err_str = TRIM(ADJUSTL(err_str)) // &
              TRIM(ADJUSTL(stmp)) // "'."
    CALL cp_assert(actlen==msglen,cp_failure_level,&
         cp_assertion_failed, routineP, err_str)

    ! restore forces on all processors (no message passing)
    NULLIFY(m,f)
    ALLOCATE(m(helium%solute_beads,helium%solute_atoms*3),STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(f(helium%solute_beads,helium%solute_atoms*3),STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)
    m(:,:) = .TRUE.
    f(:,:) = 0.0_dp
    helium%force_avrg(:,:) = UNPACK(message(1:msglen), MASK=m, FIELD=f )
    helium%force_inst(:,:) = 0.0_dp
    DEALLOCATE(f,m,STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)

    CALL helium_write_line("Forces on the solute read from the input file.",error)

    NULLIFY(message)

    RETURN
  END SUBROUTINE helium_force_restore

  ! ***************************************************************************
  !> \brief  Initialize the permutation state.
  !> \author Lukasz Walewski
  !> \date   2009-11-05
  !> \descr  Assign the identity permutation at each processor. Inverse
  !>         permutation array gets assigned as well.
  ! ***************************************************************************
  SUBROUTINE helium_perm_init( helium, error )
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_perm_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    DO ia = 1, helium%atoms
      helium%permutation(ia) = ia
      helium%iperm(ia) = ia
    END DO

    RETURN
  END SUBROUTINE helium_perm_init

  ! ***************************************************************************
  !> \brief  Restore permutation state from the input structre.
  !> \author Lukasz Walewski
  !> \date   2009-11-05
  !> \descr  Transfer permutation state from the input tree to the runtime
  !>         data structures on each processor. Inverse permutation array is
  !>         recalculated according to the restored permutation state.
  !> \par    History
  !>         2010-07-22 accomodate additional cpus in the runtime wrt the
  !>                    restart [lwalewski]
  ! ***************************************************************************
  SUBROUTINE helium_perm_restore( helium, error )
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_perm_restore', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: err_str, stmp
    INTEGER                                  :: actlen, ia, ic, msglen, &
                                                num_env_restart, offset
    INTEGER, DIMENSION(:), POINTER           :: message
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    ! assign the pointer to the memory location of the input structure, where
    ! the permutation state is stored
    NULLIFY(message)
    CALL section_vals_val_get( helium%input, &
         "MOTION%PINT%HELIUM%PERM%_DEFAULT_KEYWORD_", &
         i_vals=message, error=error)

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
      CALL helium_write_line(err_str,error)
      err_str = "Number of environments in the restart...: '"
      stmp = ""
      WRITE(stmp,*) num_env_restart
      err_str = TRIM(ADJUSTL(err_str)) // TRIM(ADJUSTL(stmp)) // "'."
      CALL helium_write_line(err_str,error)
      err_str = "Number of current run time environments.: '"
      stmp = ""
      WRITE(stmp,*) logger%para_env%num_pe
      err_str = TRIM(ADJUSTL(err_str)) // TRIM(ADJUSTL(stmp)) // "'."
      CALL helium_write_line(err_str,error)
      IF ( num_env_restart .LT. logger%para_env%num_pe ) THEN
        err_str = "Replicated permutation state from the restarted environments."
        CALL helium_write_line(err_str,error)
      END IF
      IF ( num_env_restart .GT. logger%para_env%num_pe ) THEN
        err_str = "Dropped permutation state from some restarted environments."
        CALL helium_write_line(err_str,error)
      END IF
      err_str = "Done."
      CALL helium_write_line(err_str,error)
    ELSE
      CALL helium_write_line("Permutation state read from the input file.",error)
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
    err_str = TRIM(ADJUSTL(err_str)) // &
              TRIM(ADJUSTL(stmp)) // ")."
    CALL cp_assert(ic==msglen,cp_failure_level,&
         cp_assertion_failed, routineP, err_str)

    NULLIFY(message)

    RETURN
  END SUBROUTINE helium_perm_restore

  ! ***************************************************************************
  !> \brief  Create RNG streams and initialize their state.
  !> \author Lukasz Walewski
  !> \date   2009-11-04
  !> \todo   This function shouldn't create (allocate) objects! Only
  !>         initialization, i.e. setting the seed values etc, should be done
  !>         here, allocation should be moved to helium_create
  ! ***************************************************************************
  SUBROUTINE helium_rng_state_init( helium, error )
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_rng_state_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: rank
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(rng_stream_type), POINTER           :: next_rngs, prev_rngs

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

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
         extended_precision=.TRUE.,&
         error=error)
    IF (logger%para_env%mepos .EQ. 0) THEN
      helium%rng_stream_uniform => prev_rngs
    END IF
    CALL create_rng_stream(next_rngs,&
         name="helium_rns_gaussian",&
         last_rng_stream=prev_rngs,&
         distribution_type=GAUSSIAN,&
         extended_precision=.TRUE.,&
         error=error)
    IF (logger%para_env%mepos .EQ. 0) THEN
      helium%rng_stream_gaussian => next_rngs
      NULLIFY(prev_rngs)
    ELSE
      CALL delete_rng_stream(prev_rngs,error=error)
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
           extended_precision=.TRUE.,&
           error=error)
      IF ( logger%para_env%mepos .EQ. rank ) THEN
        helium%rng_stream_uniform => next_rngs
      END IF
      CALL delete_rng_stream(prev_rngs,error=error)
      prev_rngs => next_rngs
      NULLIFY(next_rngs)
      CALL create_rng_stream(next_rngs,&
           name="helium_rns_gaussian",&
           last_rng_stream=prev_rngs,&
           distribution_type=GAUSSIAN,&
           extended_precision=.TRUE.,&
           error=error)
      IF ( logger%para_env%mepos .EQ. rank ) THEN
        helium%rng_stream_gaussian => next_rngs
        NULLIFY(prev_rngs)
      ELSE
        CALL delete_rng_stream(prev_rngs,error=error)
        prev_rngs => next_rngs
      END IF
      NULLIFY(next_rngs)
    END DO

    RETURN
  END SUBROUTINE helium_rng_state_init

  ! ***************************************************************************
  !> \brief  Restore RNG state from the input structure.
  !> \author Lukasz Walewski
  !> \date   2009-11-04
  !> \par    History
  !>         2010-07-22 Create new rng streams if more cpus available in the
  !>         runtime than in the restart [lwalewski]
  ! ***************************************************************************
  SUBROUTINE helium_rng_state_restore( helium, error )
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_rng_state_restore', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: err_str, stmp
    INTEGER                                  :: actlen, msglen, &
                                                num_env_restart, offset, rank
    LOGICAL                                  :: failure, lbf
    LOGICAL, DIMENSION(3, 2)                 :: m
    REAL(kind=dp)                            :: bf, bu
    REAL(kind=dp), DIMENSION(3, 2)           :: bg, cg, f, ig
    REAL(kind=dp), DIMENSION(:), POINTER     :: message
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(rng_stream_type), POINTER           :: next_rngs, prev_rngs

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    ! assign the pointer to the memory location of the input structure
    ! where the RNG state is stored
    NULLIFY(message)
    CALL section_vals_val_get( helium%input, &
         "MOTION%PINT%HELIUM%RNG_STATE%_DEFAULT_KEYWORD_", &
         r_vals=message, error=error)

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
         buffer=bu,buffer_filled=lbf,error=error)
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
         buffer=bu,buffer_filled=lbf,error=error)

    ELSE
      ! On processors that did not receive rng state from the restart file
      ! delete rng streams created in helium_rng_state_init and create them
      ! anew, as they have been initialized with default initial state. Here
      ! the sequence of rng streams starts from the last stream from the
      ! restart file, each stream is initialized from the previously created
      ! one.

      CALL delete_rng_stream(helium%rng_stream_uniform,error)
      CALL delete_rng_stream(helium%rng_stream_gaussian,error)
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
           extended_precision=.TRUE.,&
           error=error)
      ! set it to the last uniform rng stream from the restart
      CALL set_rng_stream(prev_rngs,bg=bg,cg=cg,ig=ig,&
           buffer=bu,buffer_filled=lbf,error=error)
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
           extended_precision=.TRUE.,&
           error=error)
      ! set it to the last gaussian rng stream from the restart
      CALL set_rng_stream(next_rngs,bg=bg,cg=cg,ig=ig,&
           buffer=bu,buffer_filled=lbf,error=error)
      ! use this on the first non-restarted environment as the new gaussian rng
      IF ( logger%para_env%mepos .EQ. num_env_restart ) THEN
        helium%rng_stream_gaussian => next_rngs
        NULLIFY(prev_rngs)
      ELSE
        CALL delete_rng_stream(prev_rngs,error=error)
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
             extended_precision=.TRUE.,&
             error=error)
        IF ( logger%para_env%mepos .EQ. rank ) THEN
          helium%rng_stream_uniform => next_rngs
        END IF
        CALL delete_rng_stream(prev_rngs,error=error)
        prev_rngs => next_rngs
        NULLIFY(next_rngs)
        CALL create_rng_stream(next_rngs,&
             name="helium_rns_gaussian",&
             last_rng_stream=prev_rngs,&
             distribution_type=GAUSSIAN,&
             extended_precision=.TRUE.,&
             error=error)
        IF ( logger%para_env%mepos .EQ. rank ) THEN
          helium%rng_stream_gaussian => next_rngs
          NULLIFY(prev_rngs)
        ELSE
          CALL delete_rng_stream(prev_rngs,error=error)
          prev_rngs => next_rngs
        END IF
        NULLIFY(next_rngs)
      END DO

    END IF

    ! say what has been done
    IF ( num_env_restart .NE. logger%para_env%num_pe ) THEN
      err_str = "Reading RNG state from the input file."
      CALL helium_write_line(err_str,error)
      err_str = "Number of environments in the restart...: '"
      stmp = ""
      WRITE(stmp,*) num_env_restart
      err_str = TRIM(ADJUSTL(err_str)) // TRIM(ADJUSTL(stmp)) // "'."
      CALL helium_write_line(err_str,error)
      err_str = "Number of current run time environments.: '"
      stmp = ""
      WRITE(stmp,*) logger%para_env%num_pe
      err_str = TRIM(ADJUSTL(err_str)) // TRIM(ADJUSTL(stmp)) // "'."
      CALL helium_write_line(err_str,error)
      IF ( num_env_restart .LT. logger%para_env%num_pe ) THEN
        err_str = "Created some new RNGs from the restarted environments."
        CALL helium_write_line(err_str,error)
      END IF
      IF ( num_env_restart .GT. logger%para_env%num_pe ) THEN
        err_str = "Dropped RNG state from some restarted environments."
        CALL helium_write_line(err_str,error)
      END IF
      err_str = "Done."
      CALL helium_write_line(err_str,error)
    ELSE
      CALL helium_write_line("RNG state read from the input file.",error)
    END IF

    NULLIFY(message)

   RETURN
  END SUBROUTINE helium_rng_state_restore


  ! ***************************************************************************
  !> \brief  Restore the densities from the input structure.
  !> \author Lukasz Walewski
  !> \date   2011-06-22
  ! ***************************************************************************
  SUBROUTINE helium_densities_restore( helium, error )

    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_densities_restore', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: err_str, stmp1, stmp2
    INTEGER                                  :: istat, itmp, msglen
    LOGICAL                                  :: failure, ltmp
    LOGICAL, DIMENSION(:, :, :, :), POINTER  :: m
    REAL(kind=dp), DIMENSION(:), POINTER     :: message
    REAL(kind=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: f

    itmp = helium%rho_nbin
    ALLOCATE(helium%rho_rstr(helium%rho_num,itmp,itmp,itmp),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ! assign pointers to the memory locations in the input structure, where
    ! the old densities are stored and check if the sizes match
    NULLIFY(message)
    CALL section_vals_val_get( helium%input, &
         "MOTION%PINT%HELIUM%RHO%CUBE_DATA%_DEFAULT_KEYWORD_", &
         r_vals=message, error=error)
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
      CALL cp_assert(.FALSE.,cp_failure_level,&
           cp_assertion_failed, routineP, err_str)
    END IF

    itmp = helium%rho_nbin
    NULLIFY(m,f)
    ALLOCATE(m(helium%rho_num,itmp,itmp,itmp),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(f(helium%rho_num,itmp,itmp,itmp),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    m(:,:,:,:) = .TRUE.
    f(:,:,:,:) = 0.0_dp

    helium%rho_rstr(:,:,:,:) = UNPACK(message(1:msglen), MASK=m, FIELD=f)

    DEALLOCATE(f,m,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(message)

    RETURN
  END SUBROUTINE helium_densities_restore


  ! ***************************************************************************
  !> \brief Count atoms of different types and store their global indices.
  !> \note  Arrays ALLOCATEd here are (should be) DEALLOCATEd in
  !>        helium_release.
  !> \author Lukasz Walewski
  ! ***************************************************************************
  SUBROUTINE helium_set_solute_indices(helium, pint_env)
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env

    CHARACTER(len=*), PARAMETER :: routineN = 'helium_set_solute_indices', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: my_element_symbol
    CHARACTER(LEN=2), DIMENSION(:), POINTER  :: element
    INTEGER                                  :: i, istat, j, mnum, natoms, &
                                                nelements
    LOGICAL                                  :: found, my_failure
    REAL(KIND=dp)                            :: mass
    TYPE(cp_error_type)                      :: my_error
    TYPE(cp_subsys_type), POINTER            :: my_subsys
    TYPE(f_env_type), POINTER                :: my_f_env
    TYPE(particle_list_type), POINTER        :: my_particles

    my_failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,my_error,my_failure)

    ! set up my_particles structure
    NULLIFY(my_f_env, my_subsys, my_particles)
    CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id, &
      f_env=my_f_env, new_error=my_error, failure=my_failure)
    CALL force_env_get(force_env=my_f_env%force_env, subsys=my_subsys, &
      error=my_error)
    CALL cp_subsys_get(my_subsys, particles=my_particles, error=my_error)
    CALL f_env_rm_defaults(my_f_env, my_error, istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,my_error,my_failure)

    natoms  = helium%solute_atoms
    NULLIFY(helium%solute_element)
    ALLOCATE(helium%solute_element(natoms), STAT=istat)
    CPPostcondition(istat==0, cp_fatal_level, routineP, my_error, my_failure)

    ! in the worst case there will be as many atomic types as atoms
    NULLIFY(element)
    ALLOCATE(element(natoms), STAT=istat)
    CPPostcondition(istat==0, cp_fatal_level, routineP, my_error, my_failure)

    ! find out how many different atomic types are there
    nelements = 0
    DO i=1, natoms
      CALL get_atomic_kind( my_particles%els(i)%atomic_kind, mass=mass)
      mnum = NINT(cp_unit_from_cp2k(mass, "amu", error=my_error))
      SELECT CASE (mnum)
        CASE (1)
          my_element_symbol = "H "
        CASE (16)
          my_element_symbol = "O "
        CASE (35)
          my_element_symbol = "CL"
        CASE DEFAULT
          WRITE(*,*) "Atomic mass number'", mnum, &
            "' not supported by the HELIUM code."
          CPAssert(.FALSE., cp_failure_level, routineP, my_error, my_failure)
      END SELECT
      helium%solute_element(i) = my_element_symbol
      ! check if this element symbol is already present in element
      found = .FALSE.
      DO j=1, nelements
        IF ( element(j) == my_element_symbol ) THEN
          found = .TRUE.
          EXIT
        END IF
      END DO
      ! increase the nelements counter if not
      IF (.NOT. found) THEN
        nelements = nelements + 1
        element(nelements) = my_element_symbol
      END IF
    END DO
    CPPostcondition(nelements.LE.3,cp_fatal_level,routineP,my_error,my_failure)

    ! allocate the arrays, DEALLOCATE them in helium_release
    ! (solute_index a bit superfluous at the moment)
    NULLIFY(helium%solute_number,helium%solute_index)
    ALLOCATE(helium%solute_number(3), STAT=istat)
    ALLOCATE(helium%solute_index(3,natoms), STAT=istat)
    CPPostcondition(istat==0,cp_fatal_level,routineP,my_error,my_failure)

    ! collect atomic indices
    helium%solute_number(:) = 0
    DO i=1, natoms
      SELECT CASE (helium%solute_element(i))
      CASE ("CL")
        helium%solute_number(1) = helium%solute_number(1) + 1
        helium%solute_index(1,helium%solute_number(1)) = i
      CASE ("O ")
        helium%solute_number(2) = helium%solute_number(2) + 1
        helium%solute_index(2,helium%solute_number(2)) = i
      CASE ("H ")
        helium%solute_number(3) = helium%solute_number(3) + 1
        helium%solute_index(3,helium%solute_number(3)) = i
      CASE DEFAULT
        WRITE(*,*) "Atom type '", helium%solute_element(i), &
          "' not supported by the HELIUM code."
        CPAssert(.FALSE., cp_failure_level, routineP, my_error, my_failure)
      END SELECT
    END DO

    DEALLOCATE(element, STAT=istat)
    CPPostconditionNoFail(istat==0, cp_warning_level, routineP, my_error)

    RETURN
  END SUBROUTINE helium_set_solute_indices

  ! ***************************************************************************
  !> \brief Sets helium%solute_cell based on the solute's force_env.
  !> \author Lukasz Walewski
  !> \descr The simulation cell for the solvated molecule is taken from force_env
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

    INTEGER                                  :: status
    LOGICAL                                  :: my_failure, my_orthorhombic
    TYPE(cell_type), POINTER                 :: my_cell
    TYPE(cp_error_type)                      :: my_error
    TYPE(f_env_type), POINTER                :: my_f_env

    my_failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,my_error,my_failure)

    ! get the cell structure from pint_env
    NULLIFY(my_f_env, my_cell)
    CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id, &
      f_env=my_f_env, new_error=my_error, failure=my_failure)
    CALL force_env_get(force_env=my_f_env%force_env, cell=my_cell, &
      error=my_error)
    CALL f_env_rm_defaults(my_f_env, my_error, status)
    CPPostcondition(status==0,cp_failure_level,routineP,my_error,my_failure)

    CALL get_cell(my_cell, orthorhombic=my_orthorhombic)
    IF (.NOT. my_orthorhombic) THEN
      PRINT *, "Helium solvent not implemented for non-orthorhombic cells."
      CPAssert(.FALSE., cp_failure_level, routineP, my_error, my_failure)
    ELSE
      helium%solute_cell => my_cell
    END IF

    RETURN
  END SUBROUTINE helium_set_solute_cell

!TODO headers/comments, beautify

FUNCTION helium_atom_action(helium,n,i) RESULT(res)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: n, i
    REAL(KIND=dp)                            :: res

    INTEGER                                  :: c, j
    REAL(KIND=dp)                            :: r(3), rp(3), s, t

   s = 0.0_dp
   t = 0.0_dp
   IF (n < helium%beads) THEN
      DO c = 1, 3
         r(c) = helium%pos(c,i,n) - helium%pos(c,i,n+1)
      END DO
      CALL helium_pbc( helium, r )
      t = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
      DO j = 1, i-1
         DO c = 1, 3
            r(c)  = helium%pos(c,i,n) - helium%pos(c,j,n)
            rp(c) = helium%pos(c,i,n+1) - helium%pos(c,j,n+1)
         END DO
         CALL helium_pbc( helium, r )
         CALL helium_pbc( helium, rp )
         s = s + helium_eval_expansion(helium,r,rp,helium%uij,1)
      END DO
      DO j = i+1, helium%atoms
         DO c = 1, 3
            r(c)  = helium%pos(c,i,n) - helium%pos(c,j,n)
            rp(c) = helium%pos(c,i,n+1) - helium%pos(c,j,n+1)
         END DO
         CALL helium_pbc( helium, r )
         CALL helium_pbc( helium, rp )
         s = s + helium_eval_expansion(helium,r,rp,helium%uij,1)
      END DO
   ELSE
      DO c = 1, 3
         r(c) = helium%pos(c,i,n) - helium%pos(c,helium%permutation(i),1)
      END DO
      CALL helium_pbc( helium, r )
      t = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
      DO j = 1, i-1
         DO c = 1, 3
            r(c)  = helium%pos(c,i,n) - helium%pos(c,j,n)
            rp(c) = helium%pos(c,helium%permutation(i),1) - helium%pos(c,helium%permutation(j),1)
         END DO
         CALL helium_pbc( helium, r )
         CALL helium_pbc( helium, rp )
         s = s + helium_eval_expansion(helium,r,rp,helium%uij,1)
      END DO
      DO j = i+1, helium%atoms
         DO c = 1, 3
            r(c)  = helium%pos(c,i,n) - helium%pos(c,j,n)
            rp(c) = helium%pos(c,helium%permutation(i),1) - helium%pos(c,helium%permutation(j),1)
         END DO
         CALL helium_pbc( helium, r )
         CALL helium_pbc( helium, rp )
         s = s + helium_eval_expansion(helium,r,rp,helium%uij,1)
      END DO
   END IF
   t = t / (2.0_dp*helium%tau*helium%hb2m)
   s = s * 0.5_dp
   res = s+t
   RETURN

END FUNCTION helium_atom_action

FUNCTION helium_link_action(helium,n) RESULT(res)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: n
    REAL(KIND=dp)                            :: res

    INTEGER                                  :: c, i, j
    REAL(KIND=dp)                            :: r(3), rp(3), s, t

   s = 0.0_dp
   t = 0.0_dp
   IF (n < helium%beads) THEN
      DO i = 1, helium%atoms
         DO c = 1, 3
            r(c) = helium%pos(c,i,n) - helium%pos(c,i,n+1)
         END DO
         CALL helium_pbc( helium, r )
         t = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
         DO j = 1, i-1
            DO c = 1, 3
               r(c)  = helium%pos(c,i,n) - helium%pos(c,j,n)
               rp(c) = helium%pos(c,i,n+1) - helium%pos(c,j,n+1)
            END DO
            CALL helium_pbc( helium, r )
            CALL helium_pbc( helium, rp )
            s = s + helium_eval_expansion(helium,r,rp,helium%uij,1)
         END DO
      END DO
   ELSE
      DO i = 1, helium%atoms
         DO c = 1, 3
            r(c) = helium%pos(c,i,n) - helium%pos(c,helium%permutation(i),1)
         END DO
         CALL helium_pbc( helium, r )
         t = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
         DO j = 1, i-1
            DO c = 1, 3
               r(c)  = helium%pos(c,i,n) - helium%pos(c,j,n)
               rp(c) = helium%pos(c,helium%permutation(i),1) - helium%pos(c,helium%permutation(j),1)
            END DO
            CALL helium_pbc( helium, r )
            CALL helium_pbc( helium, rp )
            s = s + helium_eval_expansion(helium,r,rp,helium%uij,1)
         END DO
      END DO
   END IF
   t = t / (2.0_dp*helium%tau*helium%hb2m)
   res = s+t
   RETURN

END FUNCTION helium_link_action

FUNCTION helium_total_action(helium) RESULT(res)

    TYPE(helium_solvent_type), POINTER       :: helium
    REAL(KIND=dp)                            :: res

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: s

   s = 0.0_dp
   DO i = 1, helium%beads
      s = s + helium_link_action(helium,i)
   END DO
   res = s
   RETURN

END FUNCTION helium_total_action

SUBROUTINE helium_delta_pos(helium,part,ref_bead,delta_bead,d)

    TYPE(helium_solvent_type), POINTER       :: helium
    INTEGER, INTENT(IN)                      :: part, ref_bead, delta_bead
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: d

    INTEGER                                  :: b, bead, db, n, nbead, np, p
    REAL(KIND=dp), DIMENSION(3)              :: r

   b = helium%beads
   n = helium%atoms

   d(:) = 0.0_dp
   IF (delta_bead > 0) THEN
      bead = ref_bead
      p = part
      db = delta_bead
      DO
         IF (db < 1) EXIT
         nbead = bead + 1
         np = p
         IF (nbead > b) THEN
            nbead = nbead - b
            np = helium%permutation(np)
         END IF
         r(:) = helium%pos(:,p,bead) - helium%pos(:,np,nbead)
         CALL helium_pbc( helium, r )
         d(:) = d(:) + r(:)
         bead = nbead
         p = np
         db = db-1
      END DO
   ELSEIF ( delta_bead < 0) THEN
      bead = ref_bead
      p = part
      db = delta_bead
      DO
         IF (db >= 0) EXIT
         nbead = bead - 1
         np = p
         IF (nbead < 1) THEN
            nbead = nbead + b
            np = helium%iperm(np)
         END IF
         r(:) = helium%pos(:,p,bead) - helium%pos(:,np,nbead)
         CALL helium_pbc( helium, r )
         d(:) = d(:) + r(:)
         bead = nbead
         p = np
         db = db + 1
      END DO
   END IF
   RETURN
END SUBROUTINE helium_delta_pos

END MODULE helium_methods
