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

   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_result_methods,               ONLY: get_results
   USE cp_result_types,                 ONLY: cp_result_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mcl_api,                         ONLY: mcl_finalize,&
                                              mcl_get_api_version,&
                                              mcl_get_program_id,&
                                              mcl_receive,&
                                              mcl_send
   USE mcl_requests,                    ONLY: MCL_DATA,&
                                              MCL_LENGTH,&
                                              MCL_REQUEST,&
                                              MCL_RUNTYPE_QM_RS_GRID
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_list_types,             ONLY: particle_list_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mimic_communicator'

! **************************************************************************************************
!> \brief MiMiC communicator class that facilitates MiMiC client-server data exchange
!> \par History
!>      05.2025 Created [AA]
! **************************************************************************************************
   TYPE, PUBLIC :: mimic_communicator_type
      PRIVATE
      !> communication
      TYPE(mp_para_env_type), POINTER                 :: para_env => Null()
      LOGICAL                                         :: is_ionode = .FALSE.
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
      INTEGER, DIMENSION(:, :), ALLOCATABLE            :: npts_pproc
      !> beginning index of the local buffer in the global buffer diminished by 1
      INTEGER, DIMENSION(:), ALLOCATABLE              :: lb_pproc

   CONTAINS

      PROCEDURE :: initialize
      PROCEDURE :: finalize
      PROCEDURE :: receive_request
      PROCEDURE :: send_value
      PROCEDURE :: send_client_info
      PROCEDURE :: send_atom_info
      PROCEDURE :: send_kind_info
      PROCEDURE :: send_box_info
      PROCEDURE :: send_result
      PROCEDURE :: send_grid_coordinates
      PROCEDURE :: send_density
      PROCEDURE :: send_forces
      PROCEDURE :: send_positions
      PROCEDURE :: receive_positions
      PROCEDURE :: receive_potential

   END TYPE mimic_communicator_type

CONTAINS

! **************************************************************************************************
!> \brief Initialize the communicator by loading data and saving pointers to relevant data
!> \param this ...
!> \param force_env ...
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

      ALLOCATE (this%npts_pproc(3, 0:this%para_env%num_pe - 1), source=0)
      this%npts_pproc(:, this%para_env%mepos) = this%pw_info%pw_grid%npts_local
      CALL this%para_env%sum(this%npts_pproc)

      ALLOCATE (this%lb_pproc(0:this%para_env%num_pe - 1), source=0)
      this%lb_pproc(this%para_env%mepos) = this%pw_info%pw_grid%bounds_local(1, 1) &
                                           - this%pw_info%pw_grid%bounds(1, 1)
      CALL this%para_env%sum(this%lb_pproc)

      this%n_spins = dft_control%nspins

      CALL set_qs_env(qs_env, mimic=.TRUE.)
      dft_control%apply_external_potential = .TRUE.
      dft_control%eval_external_potential = .FALSE.

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
!> \brief Finalize the simulation by deallocating memory
!> \param this ...
! **************************************************************************************************
   SUBROUTINE finalize(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':finalize'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL this%para_env%sync()

      CALL mcl_finalize()

      CALL timestop(handle)

   END SUBROUTINE finalize

! **************************************************************************************************
!> \brief Receive a request from the server
!> \param this ...
!> \return ...
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
!> \brief Send the specified single value data to the server
!> \param this ...
!> \param option word corresponding to available options
!> \note Several values hardcoded for now
! **************************************************************************************************
   SUBROUTINE send_value(this, option)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      CHARACTER(LEN=*)                                   :: option

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_value'

      REAL(dp)                                           :: energy
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (option)
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
!> \brief Send the specified information about the client to the server
!> \param this ...
!> \param option word corresponding to available options
! **************************************************************************************************
   SUBROUTINE send_client_info(this, option)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      CHARACTER(LEN=*)                                   :: option

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_client_info'

      CHARACTER(LEN=*), PARAMETER                        :: client_name = "CP2K"
      INTEGER, DIMENSION(3)                              :: api_version
      INTEGER                                            :: handle, length

      CALL timeset(routineN, handle)

      SELECT CASE (option)
      CASE ("id")
         CALL mcl_send(this%client_id, 1, MCL_DATA, this%mcl_server)
      CASE ("name")
         length = LEN(client_name)
         CALL mcl_send(length, 1, MCL_LENGTH, this%mcl_server)
         CALL mcl_send(client_name, length, MCL_DATA, this%mcl_server)
      CASE ("run_type")
         CALL mcl_send(MCL_RUNTYPE_QM_RS_GRID, 1, MCL_DATA, this%mcl_server)
      CASE ("api_version")
         CALL mcl_get_api_version(api_version)
         CALL mcl_send(api_version, 3, MCL_DATA, this%mcl_server)
      CASE DEFAULT
         CPABORT("The value chosen in "//routineN//" is not implemented.")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE send_client_info

! **************************************************************************************************
!> \brief Send the specified data for each atom to the server
!> \param this ...
!> \param option word corresponding to available options
! **************************************************************************************************
   SUBROUTINE send_atom_info(this, option)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      CHARACTER(LEN=*)                                   :: option

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_atom_info'

      INTEGER, DIMENSION(:), ALLOCATABLE                 :: buffer
      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      ALLOCATE (buffer(this%n_atoms))
      SELECT CASE (option)
      CASE ("kinds")
         DO i = 1, this%n_atoms
            buffer(i) = this%atoms%els(i)%atomic_kind%kind_number
         END DO
      CASE ("ids")
         DO i = 1, this%n_atoms
            buffer(i) = this%atoms%els(i)%atom_index
         END DO
      CASE DEFAULT
         CPABORT("The value chosen in "//routineN//" is not implemented.")
      END SELECT
      CALL mcl_send(buffer, SIZE(buffer), MCL_DATA, this%mcl_server)

      CALL timestop(handle)

   END SUBROUTINE send_atom_info

! **************************************************************************************************
!> \brief Send the specified data for each kind to the server
!> \param this ...
!> \param option word corresponding to available options
! **************************************************************************************************
   SUBROUTINE send_kind_info(this, option)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      CHARACTER(LEN=*)                                   :: option

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_kind_info'

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kinds
      REAL(dp), DIMENSION(:), ALLOCATABLE                :: buffer_dp
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: buffer_i
      CHARACTER(LEN=:), ALLOCATABLE                      :: labels
      CHARACTER(LEN=default_string_length)               :: label
      INTEGER                                            :: handle, length, i

      CALL timeset(routineN, handle)

      SELECT CASE (option)
      CASE ("labels")
         ALLOCATE (CHARACTER(30*this%n_kinds) :: labels)
         labels = ""
         DO i = 1, this%n_kinds
            CALL get_atomic_kind(this%kinds%els(i), name=label)
            labels = TRIM(labels)//TRIM(label)//","
         END DO
         length = LEN(TRIM(labels)) - 1
         CALL mcl_send(length, 1, MCL_LENGTH, this%mcl_server)
         CALL mcl_send(labels, length, MCL_DATA, this%mcl_server)
      CASE ("elements")
         ALLOCATE (buffer_i(this%n_kinds))
         DO i = 1, this%n_kinds
            CALL get_atomic_kind(this%kinds%els(i), z=buffer_i(i))
         END DO
         CALL mcl_send(buffer_i, SIZE(buffer_i), MCL_DATA, this%mcl_server)
      CASE ("masses")
         ALLOCATE (buffer_dp(this%n_kinds))
         DO i = 1, this%n_kinds
            buffer_dp(i) = cp_unit_from_cp2k(this%kinds%els(i)%mass, "AMU")
         END DO
         CALL mcl_send(buffer_dp, SIZE(buffer_dp), MCL_DATA, this%mcl_server)
      CASE ("nuclear_charges")
         NULLIFY (qs_env, qs_kinds)
         CALL force_env_get(this%force_env, qs_env=qs_env)
         CALL get_qs_env(qs_env, qs_kind_set=qs_kinds)
         ALLOCATE (buffer_dp(this%n_kinds))
         DO i = 1, this%n_kinds
            CALL get_qs_kind(qs_kinds(i), zeff=buffer_dp(i))
         END DO
         CALL mcl_send(buffer_dp, SIZE(buffer_dp), MCL_DATA, this%mcl_server)
      CASE DEFAULT
         CPABORT("The value chosen in "//routineN//" is not implemented.")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE send_kind_info

! **************************************************************************************************
!> \brief Send the specified box information to the server
!> \param this ...
!> \param option word corresponding to available options
! **************************************************************************************************
   SUBROUTINE send_box_info(this, option)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      CHARACTER(LEN=*)                                   :: option

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_box_info'

      INTEGER, DIMENSION(3)                              :: npts_glob
      REAL(dp), DIMENSION(3)                             :: origin
      REAL(dp), DIMENSION(9)                             :: box_vectors
      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      npts_glob = this%pw_info%pw_grid%npts

      SELECT CASE (option)
      CASE ("num_gridpoints")
         CALL mcl_send(npts_glob, 3, MCL_DATA, this%mcl_server)
      CASE ("origin")
         origin = 0.0_dp
         CALL mcl_send(origin, 3, MCL_DATA, this%mcl_server)
      CASE ("box_vectors")
         box_vectors = [(this%pw_info%pw_grid%dh(:, i)*REAL(npts_glob(i), dp), i=1, 3)]
         CALL mcl_send(box_vectors, 9, MCL_DATA, this%mcl_server)
      CASE DEFAULT
         CPABORT("The value chosen in "//routineN//" is not implemented.")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE send_box_info

! **************************************************************************************************
!> \brief Send the specified result to the server
!> \param this ...
!> \param option word corresponding to available options
! **************************************************************************************************
   SUBROUTINE send_result(this, option)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this
      CHARACTER(LEN=*)                                   :: option

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_result'

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_result_type), POINTER                      :: results
      CHARACTER(LEN=default_string_length)               :: description
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: buffer
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      NULLIFY (qs_env, results)
      CALL force_env_get(this%force_env, qs_env=qs_env)
      CALL get_qs_env(qs_env, results=results)

      SELECT CASE (option)
      CASE ("hirshfeld_charges")
         description = "[HIRSHFELD-CHARGES]"
         ALLOCATE (buffer(this%n_atoms), source=0.0_dp)
         CALL get_results(results, description, buffer)
         CALL mcl_send(buffer, SIZE(buffer), MCL_DATA, this%mcl_server)
      CASE DEFAULT
         CPABORT("The value chosen in "//routineN//" is not implemented.")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE send_result

! **************************************************************************************************
!> \brief Send grid point coordinates to the server
!> \param this ...
! **************************************************************************************************
   SUBROUTINE send_grid_coordinates(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_grid_coordinates'

      INTEGER, DIMENSION(3)                              :: npts_glob, npts, lb_glob, lb, ub
      REAL(dp), DIMENSION(3)                             :: origin
      REAL(dp), DIMENSION(3, 3)                           :: box_vectors
      REAL(dp), DIMENSION(:, :), ALLOCATABLE              :: coords
      INTEGER                                            :: handle, i, j, k, offset

      CALL timeset(routineN, handle)

      origin = 0.0_dp
      box_vectors = this%pw_info%pw_grid%dh
      ! number of grid points
      npts_glob = this%pw_info%pw_grid%npts
      npts = this%pw_info%pw_grid%npts_local
      ! bounds
      lb_glob = this%pw_info%pw_grid%bounds(1, :)
      lb = this%pw_info%pw_grid%bounds_local(1, :)
      ub = this%pw_info%pw_grid%bounds_local(2, :)

      ALLOCATE (coords(3, PRODUCT(npts_glob)), source=0.0_dp)
      offset = (lb(1) - lb_glob(1))*PRODUCT(npts(2:))
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               offset = offset + 1
               coords(:, offset) = origin + box_vectors(:, 1)*REAL(i - lb_glob(1), dp) &
                                   + box_vectors(:, 2)*REAL(j - lb_glob(2), dp) &
                                   + box_vectors(:, 3)*REAL(k - lb_glob(3), dp)
            END DO
         END DO
      END DO
      CALL this%para_env%sum(coords)

      CALL mcl_send(coords, SIZE(coords), MCL_DATA, this%mcl_server)

      CALL timestop(handle)

   END SUBROUTINE send_grid_coordinates

! **************************************************************************************************
!> \brief Receive external potential from the server
!> \param this ...
! **************************************************************************************************
   SUBROUTINE receive_potential(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':receive_potential'

      INTEGER, DIMENSION(3)                              :: npts, lb, ub
      REAL(dp), DIMENSION(:), ALLOCATABLE, TARGET        :: buffer
      REAL(dp), DIMENSION(:), ALLOCATABLE                :: buffer_loc
      REAL(dp), DIMENSION(:), POINTER                    :: buffer_p
      INTEGER                                            :: i, j, k, i_proc, offset
      INTEGER                                            :: handle, length, tag

      CALL timeset(routineN, handle)

      NULLIFY (buffer_p)
      npts = this%pw_info%pw_grid%npts_local
      lb = this%pw_info%pw_grid%bounds_local(1, :)
      ub = this%pw_info%pw_grid%bounds_local(2, :)
      ALLOCATE (buffer_loc(PRODUCT(npts)))

      tag = 1

      IF (this%is_ionode) THEN
         ALLOCATE (buffer(PRODUCT(this%pw_info%pw_grid%npts)))
         ! receive potential at the IO process
         CALL mcl_receive(buffer, SIZE(buffer), MCL_DATA, this%mcl_server)
         ! distribute across processes
         DO i_proc = 0, this%para_env%num_pe - 1
            length = PRODUCT(this%npts_pproc(:, i_proc))
            offset = this%lb_pproc(i_proc)*PRODUCT(npts(2:)) + 1
            buffer_p => buffer(offset:offset + length - 1)
            IF (i_proc /= this%para_env%source) THEN
               i = i_proc
               CALL this%para_env%send(buffer_p, i, tag)
            ELSE
               buffer_loc(:) = buffer_p
            END IF
         END DO
      ELSE
         CALL this%para_env%recv(buffer_loc, this%para_env%source, tag)
      END IF

      ! set the potential
      offset = 0
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               offset = offset + 1
               this%potential%array(i, j, k) = -buffer_loc(offset)
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE receive_potential

! **************************************************************************************************
!> \brief Send electron density to the server
!> \param this ...
! **************************************************************************************************
   SUBROUTINE send_density(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_density'

      INTEGER, DIMENSION(3)                              :: npts, lb, ub
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho
      REAL(dp), DIMENSION(:), ALLOCATABLE, TARGET        :: buffer
      REAL(dp), DIMENSION(:), ALLOCATABLE                :: buffer_loc
      REAL(dp), DIMENSION(:), POINTER                    :: buffer_p
      INTEGER                                            :: i_spin, i_proc, i, j, k, offset
      INTEGER                                            :: handle, length, tag

      CALL timeset(routineN, handle)

      NULLIFY (rho, buffer_p)
      CALL qs_rho_get(this%density, rho_r=rho)
      npts = this%pw_info%pw_grid%npts_local
      lb = this%pw_info%pw_grid%bounds_local(1, :)
      ub = this%pw_info%pw_grid%bounds_local(2, :)
      ALLOCATE (buffer_loc(PRODUCT(npts)))

      ! gather density values
      buffer_loc = 0.0_dp
      DO i_spin = 1, this%n_spins
         offset = 0
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  offset = offset + 1
                  buffer_loc(offset) = buffer_loc(offset) + rho(i_spin)%array(i, j, k)
               END DO
            END DO
         END DO
      END DO

      tag = 1

      IF (.NOT. this%is_ionode) THEN
         CALL this%para_env%send(buffer_loc, this%para_env%source, tag)
      ELSE
         ALLOCATE (buffer(PRODUCT(this%pw_info%pw_grid%npts)))
         ! collect from the processes at the IO process
         DO i_proc = 0, this%para_env%num_pe - 1
            length = PRODUCT(this%npts_pproc(:, i_proc))
            offset = this%lb_pproc(i_proc)*PRODUCT(npts(2:)) + 1
            buffer_p => buffer(offset:offset + length - 1)
            IF (i_proc /= this%para_env%source) THEN
               i = i_proc
               CALL this%para_env%recv(buffer_p, i, tag)
            ELSE
               buffer_p = buffer_loc
            END IF
         END DO
         ! send the density
         CALL mcl_send(buffer, SIZE(buffer), MCL_DATA, this%mcl_server)
      END IF

      CALL timestop(handle)

   END SUBROUTINE send_density

! **************************************************************************************************
!> \brief Send positions of all atoms to the server
!> \param this ...
! **************************************************************************************************
   SUBROUTINE send_positions(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_positions'

      REAL(dp), DIMENSION(:, :), ALLOCATABLE             :: buffer
      INTEGER                                            :: handle, i_atom

      CALL timeset(routineN, handle)

      ALLOCATE (buffer(3, this%n_atoms))
      DO i_atom = 1, this%n_atoms
         buffer(:, i_atom) = this%atoms%els(i_atom)%r
      END DO
      CALL mcl_send(buffer, SIZE(buffer), MCL_DATA, this%mcl_server)

      CALL timestop(handle)

   END SUBROUTINE send_positions

! **************************************************************************************************
!> \brief Receive positions of all atoms from the server
!> \param this ...
! **************************************************************************************************
   SUBROUTINE receive_positions(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':receive_positions'

      REAL(dp), DIMENSION(:, :), ALLOCATABLE             :: buffer
      INTEGER                                            :: handle, i_atom

      CALL timeset(routineN, handle)

      ALLOCATE (buffer(3, this%n_atoms))
      CALL mcl_receive(buffer, SIZE(buffer), MCL_DATA, this%mcl_server)
      CALL this%para_env%bcast(buffer)
      DO i_atom = 1, this%n_atoms
         this%atoms%els(i_atom)%r = buffer(:, i_atom)
      END DO

      CALL timestop(handle)

   END SUBROUTINE receive_positions

! **************************************************************************************************
!> \brief Send QM forces of all atoms to the server
!> \param this ...
! **************************************************************************************************
   SUBROUTINE send_forces(this)
      CLASS(mimic_communicator_type), INTENT(INOUT)      :: this

      CHARACTER(LEN=*), PARAMETER                        :: routineN = moduleN//':send_forces'

      REAL(dp), DIMENSION(:, :), ALLOCATABLE             :: buffer
      INTEGER                                            :: handle, i_atom

      CALL timeset(routineN, handle)

      ALLOCATE (buffer(3, this%n_atoms))
      DO i_atom = 1, this%n_atoms
         buffer(:, i_atom) = this%atoms%els(i_atom)%f
      END DO
      CALL mcl_send(buffer, SIZE(buffer), MCL_DATA, this%mcl_server)

      CALL timestop(handle)

   END SUBROUTINE send_forces

END MODULE mimic_communicator
