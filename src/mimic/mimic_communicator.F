!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Module containing a MiMiC communicator class
!> \par History
!>      05.2025 Created [AA]
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
   USE message_passing,                 ONLY: mp_para_env_type
   USE mcl_api,                         ONLY: mcl_finalize, &
                                              mcl_get_program_id, &
                                              mcl_receive, &
                                              mcl_send
   USE mcl_requests,                    ONLY: MCL_DATA, &
                                              MCL_REQUEST, &
                                              MCL_RUNTYPE_QM_PW
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

#include "../base/base_uses.f90"

   IMPLICIT NONE
   
   PRIVATE

   CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = 'mimic'
  
   REAL(dp), DIMENSION(:), ALLOCATABLE, SAVE, TARGET  :: grid_buffer
   REAL(dp), DIMENSION(:), ALLOCATABLE, SAVE, TARGET  :: grid_buffer_loc
   REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET  :: coord_buffer

! **************************************************************************************************
!> \brief MiMiC communicator class that facilitates MiMiC client-server data exchange 
!> \par History
!>      05.2025 Created [AA]
! **************************************************************************************************
   TYPE, PUBLIC :: mimic_communicator_type
      PRIVATE
         !> communication
         TYPE(mp_para_env_type), POINTER                 :: para_env => Null()
         LOGICAL                                         :: is_ionode
         INTEGER                                         :: mcl_server = 0, &
                                                            client_id = -1
         !> CP2K data
         TYPE(force_env_type), POINTER                   :: force_env => Null()
         TYPE(pw_pool_type), POINTER                     :: pw_info => Null()
         TYPE(particle_list_type), POINTER               :: atoms => Null()
         TYPE(atomic_kind_list_type), POINTER            :: kinds => Null()
         TYPE(qs_energy_type), POINTER                   :: energy => Null()
         TYPE(pw_r3d_rs_type), POINTER                   :: potential => Null()
         TYPE(qs_rho_type), POINTER                      :: density => Null()
         INTEGER                                         :: n_atoms = -1, &
                                                            n_kinds = -1, &
                                                            n_spins = -1
         INTEGER, DIMENSION(:,:), ALLOCATABLE            :: npts_pproc
         !> beginning index of the local buffer in the global buffer diminished by 1 
         INTEGER, DIMENSION(:), ALLOCATABLE              :: lbounds_pproc

      CONTAINS

         PROCEDURE :: initialize
         PROCEDURE :: finalize
         PROCEDURE :: receive_request
         PROCEDURE :: send_value
         PROCEDURE :: send_box_info
         PROCEDURE :: send_grid_coordinates
         PROCEDURE :: send_atom_kinds
         PROCEDURE :: send_atom_ids
         PROCEDURE :: send_kind_info
         PROCEDURE :: send_density
         PROCEDURE :: send_forces
         PROCEDURE :: send_positions
         PROCEDURE :: receive_positions
         PROCEDURE :: receive_potential

   END TYPE mimic_communicator_type

CONTAINS

! **************************************************************************************************
!> \brief Initializes MiMiC communicator loading input and by saving pointers to relevant data
!> \par History
!>      05.2025 Created [AA]
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
      ALLOCATE(coord_buffer(3, this%n_atoms))
      ALLOCATE(grid_buffer(PRODUCT(this%pw_info%pw_grid%npts)))
      ALLOCATE(grid_buffer_loc(PRODUCT(this%pw_info%pw_grid%npts_local)))

      CALL set_qs_env(qs_env, mimic=.true.)
      dft_control%apply_external_potential = .true.      
      dft_control%eval_external_potential = .false.
      
      ! allocate external electrostatic potential
      IF (ASSOCIATED(this%potential)) THEN
         CALL this%potential%release()
         DEALLOCATE (this%potential)
      END IF
      ALLOCATE (this%potential)
      CALL this%pw_info%create_pw(this%potential)
      CALL set_ks_env(ks_env, vee=this%potential)

      CALL timestop(handle)

   END SUBROUTINE initialize


! **************************************************************************************************
!> \brief Destroys the endpoint
! **************************************************************************************************
   SUBROUTINE finalize(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':finalize'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      DEALLOCATE(grid_buffer)
      DEALLOCATE(grid_buffer_loc)
      DEALLOCATE(coord_buffer)
      CALL mcl_finalize()

      CALL timestop(handle)

   END SUBROUTINE finalize


! **************************************************************************************************
!> \brief Receives a MiMiC request
! **************************************************************************************************
   FUNCTION receive_request(this) RESULT(request)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      INTEGER                                            :: request

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':receive_request'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      request = -1
      CALL mcl_receive(request, 1, MCL_REQUEST, this%mcl_server)
      CALL this%para_env%bcast(request)

      CALL timestop(handle)

   END FUNCTION receive_request


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

      CALL timeset(routineN, handle)

      SELECT CASE (option)
         CASE ("client_id")
            CALL mcl_send(this%client_id, 1, MCL_DATA, this%mcl_server)
         CASE ("client_type")
            CALL mcl_send(MCL_RUNTYPE_QM_PW, 1, MCL_DATA, this%mcl_server)
         CASE ("num_atoms", "num_atoms_in_fragments") 
            CALL mcl_send(this%n_atoms, 1, MCL_DATA, this%mcl_server)
         CASE ("num_kinds")
            CALL mcl_send(this%n_kinds, 1, MCL_DATA, this%mcl_server)
         CASE ("num_fragments")
            CALL mcl_send(1, 1, MCL_DATA, this%mcl_server)
         CASE ("num_bonds") ! later use to communicate constraints
            CALL mcl_send(0, 1, MCL_DATA, this%mcl_server)
         CASE ("num_angles") ! later use to communicate constraints
            CALL mcl_send(0, 1, MCL_DATA, this%mcl_server)
         CASE ("energy")
            energy = this%energy%total - this%energy%ee
            CALL mcl_send(energy, 1, MCL_DATA, this%mcl_server)
         CASE DEFAULT
            CPABORT("The value chosen in "//routineN//" is not implemented.")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE send_value


! **************************************************************************************************
!> \brief Sends the specified box information to MiMiC server
! **************************************************************************************************
   SUBROUTINE send_box_info(this, option)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      CHARACTER(LEN=*)                                   :: option

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_box_info'

      INTEGER                                            :: i
      INTEGER, DIMENSION(3)                              :: num_points
      REAL(dp), DIMENSION(3)                             :: origin
      REAL(dp), DIMENSION(9)                             :: box_vectors
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      num_points = this%pw_info%pw_grid%npts

      SELECT CASE (option)
         CASE ("num_gridpoints") 
            CALL mcl_send(num_points, 3, MCL_DATA, this%mcl_server)
         CASE ("origin")
            origin = 0.0_dp
            CALL mcl_send(origin, 3, MCL_DATA, this%mcl_server)
         CASE ("box_vectors")
            box_vectors = [(this%pw_info%pw_grid%dh(:,i) * REAL(num_points(i), dp), i=1,3)]
            CALL mcl_send(box_vectors, 9, MCL_DATA, this%mcl_server)
         CASE DEFAULT
            CPABORT("The value chosen in "//routineN//" is not implemented.")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE send_box_info


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

      CALL mcl_send(coords, SIZE(coords), MCL_DATA, this%mcl_server)

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
      
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      NULLIFY(buffer)
      npts = this%pw_info%pw_grid%npts_local
      bounds = this%pw_info%pw_grid%bounds_local 

      IF (this%is_ionode) THEN
         ! receive potential at the IO process
         CALL mcl_receive(grid_buffer, SIZE(grid_buffer), MCL_DATA, this%mcl_server)
         ! distribute across processes
         DO i_proc = 0, this%para_env%num_pe-1
            length = PRODUCT(this%npts_pproc(:,i_proc))
            offset = this%lbounds_pproc(i_proc) * PRODUCT(npts(2:)) + 1
            buffer => grid_buffer(offset:offset+length-1)
            IF (i_proc /= this%para_env%source) THEN
               i = i_proc
               CALL this%para_env%send(buffer, i, tag)
            ELSE
               grid_buffer_loc(:) = buffer
            END IF
         END DO
      ELSE
         CALL this%para_env%recv(grid_buffer_loc, this%para_env%source, tag)
      END IF

      ! set the potential
      offset = 0
      DO k = bounds(1,3), bounds(2,3)
         DO j = bounds(1,2), bounds(2,2)
            DO i = bounds(1,1), bounds(2,1)
               offset = offset + 1
               this%potential%array(i,j,k) = -grid_buffer_loc(offset)
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
         CALL mcl_send(grid_buffer, SIZE(grid_buffer), MCL_DATA, this%mcl_server)
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

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(this%n_atoms))
      DO i=1,this%n_atoms
         buffer(i) = this%atoms%els(i)%atomic_kind%kind_number
      END DO
      CALL mcl_send(buffer, SIZE(buffer), MCL_DATA, this%mcl_server)

      CALL timestop(handle)

   END SUBROUTINE send_atom_kinds


! **************************************************************************************************
!> \brief Send the specified data for each kind to the server
! **************************************************************************************************
   SUBROUTINE send_kind_info(this, option)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      CHARACTER(LEN=*)                                   :: option

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_kind_info'

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kinds
      REAL(dp), DIMENSION(:), ALLOCATABLE                :: buffer_dp
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: buffer_i
      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      
      SELECT CASE (option)
         CASE ("elements") 
            ALLOCATE(buffer_i(this%n_kinds))
            DO i=1,this%n_kinds
               CALL get_atomic_kind(this%kinds%els(i), z=buffer_i(i))
            END DO
            CALL mcl_send(buffer_i, SIZE(buffer_i), MCL_DATA, this%mcl_server)
         CASE ("masses") 
            ALLOCATE(buffer_dp(this%n_kinds))
            DO i=1,this%n_kinds
               buffer_dp(i) = cp_unit_from_cp2k(this%kinds%els(i)%mass, "AMU")
            END DO
            CALL mcl_send(buffer_dp, SIZE(buffer_dp), MCL_DATA, this%mcl_server)
         CASE ("nuclear_charges") 
            NULLIFY(qs_env, qs_kinds)
            CALL force_env_get(this%force_env, qs_env=qs_env)
            CALL get_qs_env(qs_env, qs_kind_set=qs_kinds)
            ALLOCATE(buffer_dp(this%n_kinds))
            DO i=1,this%n_kinds
               CALL get_qs_kind(qs_kinds(i), zeff=buffer_dp(i))
            END DO
            CALL mcl_send(buffer_dp, SIZE(buffer_dp), MCL_DATA, this%mcl_server)
         CASE DEFAULT
            CPABORT("The value chosen in "//routineN//" is not implemented.")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE send_kind_info


! **************************************************************************************************
!> \brief Send a list of atom ids to the server
! **************************************************************************************************
   SUBROUTINE send_atom_ids(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_atom_ids'

      INTEGER, DIMENSION(:), ALLOCATABLE                 :: buffer
      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      ALLOCATE(buffer(this%n_atoms))
      DO i=1,this%n_atoms
         buffer(i) = this%atoms%els(i)%atom_index
      END DO
      CALL mcl_send(buffer, SIZE(buffer), MCL_DATA, this%mcl_server)

      CALL timestop(handle)

   END SUBROUTINE send_atom_ids


! **************************************************************************************************
!> \brief Send positions of all atoms to the server
! **************************************************************************************************
   SUBROUTINE send_positions(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_positions'

      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      DO i=1,this%n_atoms
         coord_buffer(:,i) = this%atoms%els(i)%r
      END DO
      CALL mcl_send(coord_buffer, SIZE(coord_buffer), MCL_DATA, this%mcl_server)

      CALL timestop(handle)

   END SUBROUTINE send_positions


! **************************************************************************************************
!> \brief Receive positions of all atoms from the server
! **************************************************************************************************
   SUBROUTINE receive_positions(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':receive_positions'

      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      CALL mcl_receive(coord_buffer, SIZE(coord_buffer), MCL_DATA, this%mcl_server)
      CALL this%para_env%bcast(coord_buffer)
      DO i=1,this%n_atoms
         this%atoms%els(i)%r = coord_buffer(:,i)
      END DO

      CALL timestop(handle)

   END SUBROUTINE receive_positions


! **************************************************************************************************
!> \brief Send QM forces of all atoms to the server
! **************************************************************************************************
   SUBROUTINE send_forces(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_forces'

      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      DO i=1,this%n_atoms
         coord_buffer(:,i) = this%atoms%els(i)%f
      END DO
      CALL mcl_send(coord_buffer, SIZE(coord_buffer), MCL_DATA, this%mcl_server)

      CALL timestop(handle)

   END SUBROUTINE send_forces

END MODULE mimic_communicator

