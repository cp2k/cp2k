!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Module containing a MiMiC communicator class
!> \par History
!>      03.2023 Created [AA]
!> \author Andrej Antalik
! **************************************************************************************************

MODULE mimic_communicator

   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get, &
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE force_env_types,                 ONLY: force_env_get, &
                                              force_env_type
   USE kinds,                           ONLY: default_path_length, &
                                              dp
   USE machine,                         ONLY: m_getcwd
   USE message_passing,                 ONLY: mp_para_env_type
   USE mimic_utils,                     ONLY: mimic_get_input_val
   USE particle_list_types,             ONLY: particle_list_type
   USE pw_env_types,                    ONLY: pw_env_type, &
                                              pw_env_get
   USE pw_pool_types,                   ONLY: pw_pool_create, &
                                              pw_pool_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_environment_types,            ONLY: qs_environment_type, &
                                              get_qs_env, &
                                              set_qs_env
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_kind_types,                   ONLY: get_qs_kind, &
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_rho_types,                    ONLY: qs_rho_type, &
                                              qs_rho_get

   USE mimic_wrapper,                   ONLY: mimic_finalize, &
                                              mimic_receive, &
                                              mimic_send
#if defined(__MIMIC)
   USE mcl
#endif

#include "../base/base_uses.f90"

   IMPLICIT NONE
   
   PRIVATE

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = 'mimic_communicator'
  
   REAL(dp), DIMENSION(:), ALLOCATABLE, SAVE, TARGET  :: grid_buffer
   REAL(dp), DIMENSION(:), ALLOCATABLE, SAVE, TARGET  :: grid_buffer_loc
   REAL(dp), DIMENSION(:), ALLOCATABLE, SAVE, TARGET  :: coord_buffer

! **************************************************************************************************
!> \brief MiMiC communicator class that facilitates MiMiC client-server data exchange 
!> \par History
!>      03.2023 Created [AA]
! **************************************************************************************************
   TYPE, PUBLIC :: mimic_communicator_type
      PRIVATE
         !> communication
         TYPE(mp_para_env_type), POINTER                 :: para_env => Null()
         LOGICAL                                         :: is_ionode
         INTEGER                                         :: command = -1
         INTEGER                                         :: MCL_DATA = -1, &
                                                            MCL_REQUEST = -1, &
                                                            mcl_server = 0, &
                                                            client_id = -1
         !> CP2K data
         TYPE(force_env_type), POINTER                   :: force_env => Null()
         TYPE(pw_pool_type), POINTER                     :: pw_info => Null()
         TYPE(particle_list_type), POINTER               :: atoms => Null()
         TYPE(atomic_kind_list_type), POINTER            :: kinds => Null()
         TYPE(qs_energy_type), POINTER                   :: energy => Null()
         TYPE(pw_r3d_rs_type), POINTER                   :: potential => Null()
         TYPE(qs_rho_type), POINTER                      :: density => Null()
         !TYPE(pw_type), DIMENSION(:), POINTER            :: density => Null()
         INTEGER                                         :: n_atoms = -1, &
                                                            n_kinds = -1, &
                                                            n_spins = -1
         INTEGER, DIMENSION(:,:), ALLOCATABLE            :: npts_pproc
         !> beginning index of the local buffer in the global buffer diminished by 1 
         INTEGER, DIMENSION(:), ALLOCATABLE              :: lbounds_pproc

      CONTAINS

         PROCEDURE          :: initialize
         PROCEDURE          :: finalize
         PROCEDURE          :: receive_command
         PROCEDURE          :: execute_command
         PROCEDURE          :: is_cmd
         PROCEDURE, PRIVATE :: send_value
         PROCEDURE, PRIVATE :: send_cell_info
         PROCEDURE, PRIVATE :: send_grid_coordinates
         PROCEDURE, PRIVATE :: send_atom_kinds
         PROCEDURE, PRIVATE :: send_multipoles
         PROCEDURE, PRIVATE :: send_kind_masses
         PROCEDURE, PRIVATE :: send_atom_ids
         PROCEDURE, PRIVATE :: send_kind_elements
         PROCEDURE, PRIVATE :: send_nuclear_charges
         PROCEDURE, PRIVATE :: send_coords
         PROCEDURE, PRIVATE :: send_density
         PROCEDURE, PRIVATE :: receive_coords
         PROCEDURE, PRIVATE :: receive_potential
         PROCEDURE, PRIVATE :: send_forces

   END TYPE mimic_communicator_type

CONTAINS

! **************************************************************************************************
!> \brief Initializes MiMiC communicator loading input and by saving pointers to relevant data
!> \par History
!>      03.2023 Created [AA]
! **************************************************************************************************
   SUBROUTINE initialize(this, force_env)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      TYPE(force_env_type), TARGET                       :: force_env
      
      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':initialize'
      
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      INTEGER                                            :: handle
      
      CALL timeset(routineN, handle)

      CALL mcl_get_program_id(this%client_id)

      NULLIFY (subsys, qs_env, ks_env, pw_env)
      this%force_env => force_env
      CALL force_env_get(this%force_env, subsys=subsys, para_env=this%para_env, qs_env=qs_env)
      CALL cp_subsys_get(subsys, natom=this%n_atoms, particles=this%atoms, &
                                 nkind=this%n_kinds, atomic_kinds=this%kinds)
      CALL get_qs_env(qs_env, energy=this%energy, vee=this%potential, rho=this%density, &
                              dft_control=dft_control, ks_env=ks_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=this%pw_info)
      
      this%is_ionode = this%para_env%is_source()


      ALLOCATE(this%npts_pproc(3,0:this%para_env%num_pe-1), source=0)
      this%npts_pproc(: ,this%para_env%mepos) = this%pw_info%pw_grid%npts_local
      CALL this%para_env%sum(this%npts_pproc)

      ALLOCATE(this%lbounds_pproc(0:this%para_env%num_pe-1), source=0)
      this%lbounds_pproc(this%para_env%mepos) = this%pw_info%pw_grid%bounds_local(1,1) &
                                                - this%pw_info%pw_grid%bounds(1,1)
      CALL this%para_env%sum(this%lbounds_pproc)

      this%n_spins = dft_control%nspins

      ! allocate buffers
      ALLOCATE(coord_buffer(3*this%n_atoms))
      ALLOCATE(grid_buffer(PRODUCT(this%pw_info%pw_grid%npts)))
      ALLOCATE(grid_buffer_loc(PRODUCT(this%pw_info%pw_grid%npts_local)))

      CALL set_qs_env(qs_env, mimic=.true.)
      dft_control%apply_external_potential = .true.      
      dft_control%eval_external_potential = .false.
      
      ! cAAx
      ! allocate external electrostatic potential
      IF (ASSOCIATED(this%potential)) THEN
         CALL this%potential%release()
         DEALLOCATE (this%potential)
      END IF
      ALLOCATE (this%potential)
      CALL this%pw_info%create_pw(this%potential)
      CALL set_ks_env(ks_env, vee=this%potential)

#if defined(__MIMIC)
      this%MCL_DATA = MCL_DATA
      this%MCL_REQUEST = MCL_REQUEST
#endif

      CALL timestop(handle)

   END SUBROUTINE initialize


! **************************************************************************************************
!> \brief Destroys the endpoint
! **************************************************************************************************
   SUBROUTINE finalize(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':finalize'

      INTEGER                                            :: handle

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      DEALLOCATE(grid_buffer)
      DEALLOCATE(grid_buffer_loc)
      DEALLOCATE(coord_buffer)
      CALL mimic_finalize()

      CALL timestop(handle)

   END SUBROUTINE finalize


! **************************************************************************************************
!> \brief Receives a MiMiC command
! **************************************************************************************************
   SUBROUTINE receive_command(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':receive_command'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      this%command = -1
      IF (this%is_ionode) THEN
         CALL mimic_receive(this%command, 1, this%MCL_REQUEST, this%mcl_server)
      END IF
      CALL this%para_env%bcast(this%command)

      CALL timestop(handle)

   END SUBROUTINE receive_command


! **************************************************************************************************
!> \brief Checks whether 'command' corresponds to the command received from the server
!> \param command command to compare with the received command
!> \retval is_match true if 'command' matches the received command 'this%command'
! **************************************************************************************************
   LOGICAL FUNCTION is_cmd(this, command) RESULT(is_match)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      INTEGER, INTENT(IN)                                :: command

      is_match = (this%command == command)

   END FUNCTION is_cmd


! **************************************************************************************************
!> \brief Executes an action corresponding to one of the communication commands
!> \param executed [output] returns true if one of the listed actions has been executed
!> \par History
!>      03.2023 Created [AA]
! **************************************************************************************************
   SUBROUTINE execute_command(this, executed)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      LOGICAL, INTENT(INOUT)                             :: executed
          
      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':execute_command'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      executed = .TRUE.

#if defined(__MIMIC)
      IF (this%is_cmd(MCL_SEND_CLIENT_ID)) THEN
         CALL this%send_value("clientId")
      ELSE IF (this%is_cmd(MCL_SEND_NUM_PARTICLES)) THEN
         CALL this%send_value("nAtoms")
      ELSE IF (this%is_cmd(MCL_SEND_NUM_PARTICLE_SPECIES)) THEN
         CALL this%send_value("nKinds")
      ELSE IF (this%is_cmd(MCL_SEND_PARTICLE_SPECIES_IDS)) THEN
         CALL this%send_atom_kinds()
      ELSE IF (this%is_cmd(MCL_SEND_NUM_FRAGMENTS)) THEN
         CALL this%send_value("nFragments")
      ELSE IF (this%is_cmd(MCL_SEND_NUM_PARTICLES_IN_FRAGMENTS)) THEN
         CALL this%send_value("nAtomsInFragments")
      ELSE IF (this%is_cmd(MCL_SEND_NUM_CONSTR_BONDS)) THEN
         CALL this%send_value("nBonds")
      ELSE IF (this%is_cmd(MCL_SEND_NUM_CONSTR_ANGLES)) THEN
         CALL this%send_value("nAngles")
      ELSE IF (this%is_cmd(MCL_SEND_PARTICLE_MULTIPOLES)) THEN
         CALL this%send_multipoles()  ! fake placeholder (these are intended only for MM)
      ELSE IF (this%is_cmd(MCL_SEND_SPECIES_MASSES)) THEN
         CALL this%send_kind_masses()
      ELSE IF (this%is_cmd(MCL_SEND_PARTICLE_IDS_IN_FRAGMENTS)) THEN
         CALL this%send_atom_ids()
      ELSE IF (this%is_cmd(MCL_SEND_SPECIES_ELEMENTS)) THEN
         CALL this%send_kind_elements()
      ELSE IF (this%is_cmd(MCL_SEND_PARTICLE_POSITIONS)) THEN
         CALL this%send_coords()
      ELSE IF (this%is_cmd(MCL_RECV_PARTICLE_POSITIONS)) THEN
         CALL this%receive_coords()
      ELSE IF (this%is_cmd(MCL_SEND_PARTICLE_FORCES)) THEN
         CALL this%send_forces()
      ELSE IF (this%is_cmd(MCL_SEND_ENERGY)) THEN
         CALL this%send_value("Energy")
      ELSE IF (this%is_cmd(MCL_SEND_CLIENT_RUNTYPE)) THEN
         CALL this%send_value("ClientType")
      ELSE IF (this%is_cmd(MCL_SEND_NUCLEAR_CHARGES)) THEN
         CALL this%send_nuclear_charges()
      ELSE IF (this%is_cmd(MCL_SEND_BOX_NUM_GRIDPOINTS)) THEN
         CALL this%send_cell_info("nPoints")
      ELSE IF (this%is_cmd(MCL_SEND_BOX_ORIGIN)) THEN
         CALL this%send_cell_info("Origin")
      ELSE IF (this%is_cmd(MCL_SEND_BOX_VECTORS)) THEN
         CALL this%send_cell_info("Lattice")
      ELSE IF (this%is_cmd(MCL_SEND_BOX_GRIDPOINT_COORDS)) THEN
         CALL this%send_grid_coordinates()
      ELSE IF (this%is_cmd(MCL_SEND_DENSITY)) THEN
         CALL this%send_density()
      ELSE IF (this%is_cmd(MCL_RECV_POTENTIAL_ON_GRIDPOINTS)) THEN
         CALL this%receive_potential()
      ELSE
         executed = .FALSE.
      END IF
#endif

      CALL timestop(handle)

   END SUBROUTINE execute_command


! **************************************************************************************************
!> \brief Sends one of the single value data to MiMiC server
!> \param option word corresponding to available options 
!> \note for now hardcoded num_frag = 1  =>  nAtInFrag = nAt
! **************************************************************************************************
   SUBROUTINE send_value(this, option)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      CHARACTER(LEN=*)                                   :: option

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_value'

      REAL(dp)                                           :: energy
      INTEGER                                            :: handle

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      SELECT CASE (option)
         CASE ("clientId")
            CALL mimic_send(this%client_id, 1, this%MCL_DATA, this%mcl_server)
         CASE ("nAtoms", "nAtomsInFragments") 
            CALL mimic_send(this%n_atoms, 1, this%MCL_DATA, this%mcl_server)
         CASE ("nKinds")
            CALL mimic_send(this%n_kinds, 1, this%MCL_DATA, this%mcl_server)
         CASE ("nFragments")
            CALL mimic_send(1, 1, this%MCL_DATA, this%mcl_server)
         CASE ("nBonds") ! later use to communicate constraints
            CALL mimic_send(0, 1, this%MCL_DATA, this%mcl_server)
         CASE ("nAngles") ! later use to communicate constraints
            CALL mimic_send(0, 1, this%MCL_DATA, this%mcl_server)
         CASE ("Energy")
           ! cAA fix this later
            energy = this%energy%total - this%energy%ee
            CALL mimic_send(energy, 1, this%MCL_DATA, this%mcl_server)
         CASE ("ClientType")
            CALL mimic_send(MCL_RUNTYPE_QM_PW, 1, this%MCL_DATA, this%mcl_server)
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE send_value


! **************************************************************************************************
!> \brief Sends specified cell information to MiMiC server
! **************************************************************************************************
   SUBROUTINE send_cell_info(this, option)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      CHARACTER(LEN=*)                                   :: option

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_cell_info'

      INTEGER                                            :: i
      INTEGER, DIMENSION(3)                              :: num_points
      REAL(dp), DIMENSION(3)                             :: origin
      REAL(dp), DIMENSION(9)                             :: lattice
      INTEGER                                            :: handle

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      num_points = this%pw_info%pw_grid%npts

      SELECT CASE (option)
         CASE ("nPoints") 
            CALL mimic_send(num_points, 3, this%MCL_DATA, this%mcl_server)
         CASE ("Origin")
            origin = 0.0_dp
            CALL mimic_send(origin, 3, this%MCL_DATA, this%mcl_server)
         CASE ("Lattice")
            lattice = [(this%pw_info%pw_grid%dh(:,i) * REAL(num_points(i), dp), i=1,3)]
            CALL mimic_send(lattice, 9, this%MCL_DATA, this%mcl_server)
         CASE DEFAULT
            CPABORT("The value chosen in "//routineN//" is not implemented.")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE send_cell_info


! **************************************************************************************************
!> \brief Sends coordinates of ... grids information to MiMiC server
! **************************************************************************************************
   SUBROUTINE send_grid_coordinates(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_grid_coordinates'

      INTEGER                                            :: i, j, k, offset
      INTEGER, DIMENSION(3)                              :: npts_glob, npts, lbounds_glob
      INTEGER, DIMENSION(2,3)                            :: bounds
      REAL(dp), DIMENSION(3)                             :: origin
      REAL(dp), DIMENSION(3,3)                           :: lattice
      REAL(dp), DIMENSION(:,:), ALLOCATABLE              :: coords
      
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      npts_glob = this%pw_info%pw_grid%npts
      npts = this%pw_info%pw_grid%npts_local
      lbounds_glob = this%pw_info%pw_grid%bounds(1,:)
      bounds = this%pw_info%pw_grid%bounds_local 
      origin = 0.0_dp
      lattice = this%pw_info%pw_grid%dh
     
      ALLOCATE(coords(3, PRODUCT(npts_glob)), source=0.0_dp)
      offset = (bounds(1,1)-lbounds_glob(1)) * PRODUCT(npts(2:))
      DO k = bounds(1,3), bounds(2,3)
         DO j = bounds(1,2), bounds(2,2) 
            DO i = bounds(1,1), bounds(2,1)
              offset = offset + 1
              coords(:, offset) = origin + lattice(:,1) * REAL(i-lbounds_glob(1), dp) &
                                         + lattice(:,2) * REAL(j-lbounds_glob(2), dp) &
                                         + lattice(:,3) * REAL(k-lbounds_glob(3), dp)
            END DO
         END DO
      END DO
      CALL this%para_env%sum(coords)

      IF (this%is_ionode) THEN
         CALL mimic_send(PACK(coords, .true.), SIZE(coords), this%MCL_DATA, this%mcl_server)
      END IF

      CALL timestop(handle)

   END SUBROUTINE send_grid_coordinates


! **************************************************************************************************
!> \brief Receive external potential from the MiMiC server
! **************************************************************************************************
   SUBROUTINE receive_potential(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':receive_potential'

      INTEGER                                            :: i, j, k, i_proc
      INTEGER                                            :: offset, length, tag = 1
      INTEGER, DIMENSION(3)                              :: npts
      INTEGER, DIMENSION(2,3)                            :: bounds
      REAL(dp), DIMENSION(:), POINTER                    :: buffer
      REAL(dp), DIMENSION(:), ALLOCATABLE                :: buffer_loc
      
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      NULLIFY(buffer)
      npts = this%pw_info%pw_grid%npts_local
      bounds = this%pw_info%pw_grid%bounds_local 

      ! cAA
      ALLOCATE(buffer_loc(PRODUCT(npts)))
     
      IF (this%is_ionode) THEN
         ! receive potential at the IO process
         CALL mimic_receive(grid_buffer, SIZE(grid_buffer), this%MCL_DATA, this%mcl_server)
         ! distribute across processes
         DO i_proc = 0, this%para_env%num_pe-1
            length = PRODUCT(this%npts_pproc(:,i_proc))
            offset = this%lbounds_pproc(i_proc) * PRODUCT(npts(2:)) + 1
            buffer => grid_buffer(offset:offset+length-1)
            IF (i_proc /= this%para_env%source) THEN
               i = i_proc
               CALL this%para_env%send(buffer, i, tag)
            ELSE
               DO i = 1, size(buffer_loc)
                  buffer_loc(i) = buffer(i)
               END DO
               !grid_buffer_loc = buffer
            END IF
         END DO
      ELSE
         CALL this%para_env%recv(buffer_loc, this%para_env%source, tag)
         !CALL this%para_env%recv(grid_buffer_loc, this%para_env%source, tag)
      END IF

      ! set the potential
      offset = 0
      DO k = bounds(1,3), bounds(2,3)
         DO j = bounds(1,2), bounds(2,2)
            DO i = bounds(1,1), bounds(2,1)
               offset = offset + 1
               this%potential%array(i,j,k) = -buffer_loc(offset)
               !this%potential%array(i,j,k) = -grid_buffer_loc(offset)
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE receive_potential


! **************************************************************************************************
!> \brief Sends electron density to MiMiC server
! **************************************************************************************************
   SUBROUTINE send_density(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_density'

      INTEGER                                            :: i_spin, i_proc, i, j, k
      INTEGER                                            :: offset, length, tag = 1
      INTEGER, DIMENSION(3)                              :: npts
      INTEGER, DIMENSION(2,3)                            :: bounds
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho
      REAL(dp), DIMENSION(:), POINTER                    :: buffer
      
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      NULLIFY(rho, buffer)
      CALL qs_rho_get(this%density, rho_r=rho)
      npts = this%pw_info%pw_grid%npts_local
      bounds = this%pw_info%pw_grid%bounds_local 
     
      ! gather density values
      grid_buffer_loc = 0.0_dp
      DO i_spin = 1, this%n_spins
         offset = 0
         DO k = bounds(1,3), bounds(2,3)
            DO j = bounds(1,2), bounds(2,2) 
               DO i = bounds(1,1), bounds(2,1) 
                  offset = offset + 1
                  grid_buffer_loc(offset) = grid_buffer_loc(offset) + rho(i_spin)%array(i,j,k)
               END DO
            END DO
         END DO
      END DO

      IF (.NOT. this%is_ionode) THEN
         CALL this%para_env%send(grid_buffer_loc, this%para_env%source, tag)
      ELSE
         ! collect from the processes at the IO process
         DO i_proc = 0, this%para_env%num_pe-1
            length = PRODUCT(this%npts_pproc(:,i_proc))
            offset = this%lbounds_pproc(i_proc) * PRODUCT(npts(2:)) + 1
            buffer => grid_buffer(offset:offset+length-1)
            IF (i_proc /= this%para_env%source) THEN
               i = i_proc
               CALL this%para_env%recv(buffer, i, tag)
            ELSE
               buffer = grid_buffer_loc
            END IF
         END DO
         ! send the density
         CALL mimic_send(grid_buffer, SIZE(grid_buffer), this%MCL_DATA, this%mcl_server)
      END IF

      CALL timestop(handle)

   END SUBROUTINE send_density


! **************************************************************************************************
!> \brief Send a list of kind ids for each atom to the server
! **************************************************************************************************
   SUBROUTINE send_atom_kinds(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_atom_kinds'

      INTEGER, DIMENSION(:), ALLOCATABLE                 :: buffer
      INTEGER                                            :: handle, i

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(this%n_atoms))
      DO i=1,this%n_atoms
         buffer(i) = this%atoms%els(i)%atomic_kind%kind_number
      END DO
      CALL mimic_send(buffer, this%n_atoms, this%MCL_DATA, this%mcl_server)
      DEALLOCATE(buffer)

      CALL timestop(handle)

   END SUBROUTINE send_atom_kinds


! **************************************************************************************************
!> \brief TBD
! **************************************************************************************************
   SUBROUTINE send_multipoles(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_multipoles'

      REAL(dp), DIMENSION(:), ALLOCATABLE                :: buffer
      INTEGER                                            :: handle


      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(this%n_atoms))
      buffer = 0.0_dp ! placeholder
      CALL mimic_send(buffer, this%n_atoms, this%MCL_DATA, this%mcl_server)
      DEALLOCATE(buffer)

      CALL timestop(handle)

   END SUBROUTINE send_multipoles


! **************************************************************************************************
!> \brief Send a list of masses for each kind to the server
! **************************************************************************************************
   SUBROUTINE send_kind_masses(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_kind_masses'

      REAL(dp), DIMENSION(:), ALLOCATABLE                :: buffer
      INTEGER                                            :: handle, i

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(this%n_kinds))
      DO i=1,this%n_kinds
         buffer(i) = cp_unit_from_cp2k(this%kinds%els(i)%mass, "AMU")
      END DO
      CALL mimic_send(buffer, this%n_kinds, this%MCL_DATA, this%mcl_server)
      DEALLOCATE(buffer)

      CALL timestop(handle)

   END SUBROUTINE send_kind_masses


! **************************************************************************************************
!> \brief Send a list of atom ids to the server
! **************************************************************************************************
   SUBROUTINE send_atom_ids(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_atom_ids'

      INTEGER, DIMENSION(:), ALLOCATABLE                 :: buffer
      INTEGER                                            :: handle, i

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(this%n_atoms))
      DO i=1,this%n_atoms
         buffer(i) = this%atoms%els(i)%atom_index
      END DO
      CALL mimic_send(buffer, this%n_atoms, this%MCL_DATA, this%mcl_server)
      DEALLOCATE(buffer)

      CALL timestop(handle)

   END SUBROUTINE send_atom_ids


! **************************************************************************************************
!> \brief Send a list of atomic numbers for each kind to the server
! **************************************************************************************************
   SUBROUTINE send_kind_elements(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_kind_elements'

      INTEGER, DIMENSION(:), ALLOCATABLE                 :: buffer
      INTEGER                                            :: handle, i

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(this%n_kinds))
      DO i=1,this%n_kinds
         CALL get_atomic_kind(this%kinds%els(i), z=buffer(i))
      END DO
      CALL mimic_send(buffer, this%n_kinds, this%MCL_DATA, this%mcl_server)
      DEALLOCATE(buffer)

      CALL timestop(handle)

   END SUBROUTINE send_kind_elements


! **************************************************************************************************
!> \brief Send a list of ionic charges for each kind to the server
! **************************************************************************************************
   SUBROUTINE send_nuclear_charges(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_nuclear_charges'

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kinds
      REAL(dp), DIMENSION(:), ALLOCATABLE                :: buffer
      INTEGER                                            :: handle, i

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      NULLIFY(qs_env, qs_kinds)
      CALL force_env_get(this%force_env, qs_env=qs_env)
      CALL get_qs_env(qs_env, qs_kind_set=qs_kinds)

      ALLOCATE(buffer(this%n_kinds))
      DO i=1,this%n_kinds
         CALL get_qs_kind(qs_kinds(i), zeff=buffer(i))
      END DO
      CALL mimic_send(buffer, this%n_kinds, this%MCL_DATA, this%mcl_server)
      DEALLOCATE(buffer)

      CALL timestop(handle)

   END SUBROUTINE send_nuclear_charges


! **************************************************************************************************
! **************************************************************************************************
!> \brief Send coordinates of all atoms to the server
! **************************************************************************************************
   SUBROUTINE send_coords(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_coords'

      INTEGER                                            :: handle, i

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      DO i=1,this%n_atoms
         coord_buffer((3*i-2):(3*i)) = this%atoms%els(i)%r
      END DO
      CALL mimic_send(coord_buffer, 3*this%n_atoms, this%MCL_DATA, this%mcl_server)

      CALL timestop(handle)

   END SUBROUTINE send_coords


! **************************************************************************************************
!> \brief Receive coordinates of all atoms from the server
! **************************************************************************************************
   SUBROUTINE receive_coords(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':receive_coords'

      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      IF (this%is_ionode) THEN
         CALL mimic_receive(coord_buffer, 3*this%n_atoms, this%MCL_DATA, this%mcl_server)
      END IF
      CALL this%para_env%bcast(coord_buffer)
      DO i=1,this%n_atoms
         this%atoms%els(i)%r = coord_buffer((3*i-2):(3*i))
      END DO

      CALL timestop(handle)

   END SUBROUTINE receive_coords


! **************************************************************************************************
!> \brief Send QM forces of all atoms to the server
! **************************************************************************************************
   SUBROUTINE send_forces(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_forces'

      INTEGER                                            :: handle, i

      IF (.NOT. this%is_ionode) RETURN

      CALL timeset(routineN, handle)

      DO i=1,this%n_atoms
         coord_buffer((3*i-2):(3*i)) = this%atoms%els(i)%f
      END DO
      CALL mimic_send(coord_buffer, 3*this%n_atoms, this%MCL_DATA, this%mcl_server)

      CALL timestop(handle)

   END SUBROUTINE send_forces

END MODULE mimic_communicator

