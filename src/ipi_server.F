!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief i–PI server mode: Communication with i–PI clients
!> \par History
!>      03.2024 created
!> \author Sebastian Seidenath (sebastian.seidenath@uni-jena.de)
! **************************************************************************************************
MODULE ipi_server
   USE ISO_C_BINDING, ONLY: C_CHAR, &
                            C_DOUBLE, &
                            C_INT, &
                            C_LOC, &
                            C_NULL_CHAR, &
                            C_PTR
   USE cell_methods, ONLY: cell_create, &
                           init_cell
   USE cell_types, ONLY: cell_release, &
                         cell_type
   USE cp_external_control, ONLY: external_control
   USE cp_log_handling, ONLY: cp_logger_get_default_io_unit
   USE cp_subsys_types, ONLY: cp_subsys_get, &
                              cp_subsys_set, &
                              cp_subsys_type
   USE global_types, ONLY: global_environment_type
   USE input_section_types, ONLY: section_vals_get_subs_vals, &
                                  section_vals_type, &
                                  section_vals_val_get
   USE ipi_environment_types, ONLY: ipi_environment_type, &
                                    ipi_env_set
   USE kinds, ONLY: default_path_length, &
                    default_string_length, &
                    dp, &
                    int_4
   USE message_passing, ONLY: mp_para_env_type, &
                              mp_request_type, &
                              mp_testany
   USE particle_list_types, ONLY: particle_list_type
   USE particle_types, ONLY: particle_type
#ifndef __NO_SOCKETS
   USE sockets_interface, ONLY: writebuffer, &
                                readbuffer, &
                                uwait, &
                                open_bind_socket, &
                                listen_socket, &
                                accept_socket, &
                                close_socket, &
                                remove_socket_file
#endif
   USE virial_types, ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ipi_server'
   INTEGER, PARAMETER                   :: msglength = 12

   PUBLIC                               :: start_server, &
                                           shutdown_server, &
                                           request_forces

CONTAINS

! **************************************************************************************************
!> \brief Starts the i–PI server. Will block until it recieves a connection.
!> \param driver_section The driver section from the input file
!> \param para_env ...
!> \param ipi_env The ipi environment
!> \author Sebastian Seidenath (sebastian.seidenath@uni-jena.de)
! **************************************************************************************************
   SUBROUTINE start_server(driver_section, para_env, ipi_env)
      TYPE(section_vals_type), POINTER         :: driver_section
      TYPE(mp_para_env_type), POINTER          :: para_env
      TYPE(ipi_environment_type), POINTER      :: ipi_env

      CHARACTER(len=*), PARAMETER :: routineN = 'start_server'

#ifdef __NO_SOCKETS
      INTEGER                                  :: handle
      CALL timeset(routineN, handle)
      CPABORT("CP2K was compiled with the __NO_SOCKETS option!")
#else
      CHARACTER(len=default_path_length)       :: c_hostname, drv_hostname
      INTEGER                                  :: drv_port, handle, i_drv_unix, &
                                                  output_unit, socket, comm_socket
      CHARACTER(len=msglength)                 :: msgbuffer
      CHARACTER(len=msglength), PARAMETER      :: initmsg = "INIT"
      LOGICAL                                  :: drv_unix, ionode

      CALL timeset(routineN, handle)
      ionode = para_env%is_source()
      output_unit = cp_logger_get_default_io_unit()

      ! Read connection parameters
      CALL section_vals_val_get(driver_section, "HOST", c_val=drv_hostname)
      CALL section_vals_val_get(driver_section, "PORT", i_val=drv_port)
      CALL section_vals_val_get(driver_section, "UNIX", l_val=drv_unix)
      IF (output_unit > 0) THEN
         WRITE (output_unit, *) "@ i-PI SERVER BEING STARTED"
         WRITE (output_unit, *) "@ HOSTNAME: ", TRIM(drv_hostname)
         WRITE (output_unit, *) "@ PORT: ", drv_port
         WRITE (output_unit, *) "@ UNIX SOCKET: ", drv_unix
      END IF

      ! opens the socket
      socket = 0
      !inet = 1
      i_drv_unix = 1 ! a bit convoluted. socket.c uses a different convention...
      IF (drv_unix) i_drv_unix = 0

      c_hostname = TRIM(drv_hostname)//C_NULL_CHAR
      IF (ionode) THEN
         CALL open_bind_socket(socket, i_drv_unix, drv_port, c_hostname)
         CALL listen_socket(socket, 1_c_int)
         CALL accept_socket(socket, comm_socket)
         CALL close_socket(socket)
         CALL remove_socket_file(c_hostname)
         CALL ipi_env_set(ipi_env=ipi_env, sockfd=comm_socket)
      END IF

      ! Check if the client needs initialization
      ! We only send a meaningless message since we have no general way of
      ! knowing what the client is expecting
      CALL ask_status(comm_socket, msgbuffer)
      IF (TRIM(msgbuffer) == "NEEDINIT") THEN
         CALL writebuffer(comm_socket, initmsg, msglength)
         CALL writebuffer(comm_socket, 1) ! Bead index - just send 1
         CALL writebuffer(comm_socket, 12) ! Bits in the following message
         CALL writebuffer(comm_socket, "Initializing", 12)
      END IF

#endif

      CALL timestop(handle)

   END SUBROUTINE start_server

! **************************************************************************************************
!> \brief Shut down the i–PI server.
!> \param ipi_env The ipi environment in charge of the server
!> \author Sebastian Seidenath (sebastian.seidenath@uni-jena.de)
! **************************************************************************************************
   SUBROUTINE shutdown_server(ipi_env)
      TYPE(ipi_environment_type), POINTER                :: ipi_env

      CHARACTER(len=msglength), PARAMETER                :: msg = "EXIT"

      INTEGER                                            :: output_unit

      output_unit = cp_logger_get_default_io_unit()
      WRITE (output_unit, *) "@ i–PI: Shutting down server."
      CALL writebuffer(ipi_env%sockfd, msg, msglength)
      CALL close_socket(ipi_env%sockfd)
   END SUBROUTINE shutdown_server

! **************************************************************************************************
!> \brief Send atomic positions to a client and retrieve forces
!> \param ipi_env The ipi environment in charge of the connection
!> \author Sebastian Seidenath
! **************************************************************************************************
   SUBROUTINE request_forces(ipi_env)
      TYPE(ipi_environment_type), POINTER                :: ipi_env

      CHARACTER(len=msglength)                           :: msgbuffer
      INTEGER                                            :: comm_socket, i, nAtom, p, xyz
      REAL(kind=dp)                                      :: energy
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: forces

      i = 0
      nAtom = ipi_env%subsys%particles%n_els
      comm_socket = ipi_env%sockfd

      ! Step 1: See if the client is ready
      CALL ask_status(comm_socket, msgbuffer)
      IF (TRIM(msgbuffer) /= "READY") &
         CPABORT("i–PI: Expected READY header but recieved "//TRIM(msgbuffer))

      ! Step 2: Send cell and position data to client
      CALL send_posdata(comm_socket, subsys=ipi_env%subsys)

      ! Step 3: Ask for status, should be done now
      CALL ask_status(comm_socket, msgbuffer)
      IF (TRIM(msgbuffer) /= "HAVEDATA") &
         CPABORT("i–PI: Expected HAVEDATA header but recieved "//TRIM(msgbuffer))

      ! Step 4: Ask for data
      ALLOCATE (forces(3, nAtom))
      CALL ask_getforce(comm_socket, energy=energy, forces=forces)

      ! Step 4.5: Check for sanity
      IF (SIZE(forces) /= (nAtom*3)) THEN
         CPABORT("i–PI: Mismatch in particle number between CP2K and i–PI client")
      END IF

      ! Step 5: Return data
      DO p = 1, nAtom
         DO xyz = 1, 3
            ipi_env%subsys%particles%els(p)%f(xyz) = forces(xyz, p)
         END DO
      END DO
      CALL ipi_env_set(ipi_env=ipi_env, ipi_energy=energy, ipi_forces=forces)
   END SUBROUTINE request_forces

! **************************************************************************************************
!> \brief ...
!> \param sockfd ...
!> \param buffer ...
! **************************************************************************************************
   SUBROUTINE get_header(sockfd, buffer)
      INTEGER, INTENT(IN)                                :: sockfd
      CHARACTER(len=msglength), INTENT(OUT)              :: buffer

      INTEGER                                            :: output_unit

      CALL readbuffer(sockfd, buffer, msglength)
      output_unit = cp_logger_get_default_io_unit()
      IF (output_unit > 0) WRITE (output_unit, *) " @ i–PI Server: recieved ", TRIM(buffer)
   END SUBROUTINE get_header

! **************************************************************************************************
!> \brief ...
!> \param sockfd ...
!> \param buffer ...
! **************************************************************************************************
   SUBROUTINE ask_status(sockfd, buffer)
      INTEGER, INTENT(IN)                                :: sockfd
      CHARACTER(len=msglength), INTENT(OUT)              :: buffer

      CHARACTER(len=msglength), PARAMETER                :: msg = "STATUS"

      CALL writebuffer(sockfd, msg, msglength)
      CALL get_header(sockfd, buffer)
   END SUBROUTINE ask_status

! **************************************************************************************************
!> \brief ...
!> \param sockfd ...
!> \param energy ...
!> \param forces ...
!> \param virial ...
!> \param extra ...
! **************************************************************************************************
   SUBROUTINE ask_getforce(sockfd, energy, forces, virial, extra)
      INTEGER, INTENT(IN)                                :: sockfd
      REAL(kind=dp), INTENT(OUT)                         :: energy
      REAL(kind=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL, POINTER                               :: forces
      REAL(kind=dp), DIMENSION(3, 3), INTENT(OUT), &
         OPTIONAL                                        :: virial
      CHARACTER(len=:), INTENT(OUT), OPTIONAL, POINTER   :: extra

      CHARACTER(len=msglength), PARAMETER                :: msg = "GETFORCE"

      CHARACTER(len=:), ALLOCATABLE                      :: extra_buffer
      CHARACTER(len=msglength)                           :: msgbuffer
      INTEGER                                            :: extraLength, nAtom
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: forces_buffer
      REAL(kind=dp), DIMENSION(9)                        :: virial_buffer

      ! Exchange headers
      CALL writebuffer(sockfd, msg, msglength)
      CALL get_header(sockfd, msgbuffer)
      IF (TRIM(msgbuffer) /= "FORCEREADY") &
         CPABORT("i–PI: Expected FORCEREADY header but recieved "//TRIM(msgbuffer))

      ! Recieve data
      CALL readbuffer(sockfd, energy)
      CALL readbuffer(sockfd, nAtom)
      ALLOCATE (forces_buffer(3*nAtom))
      CALL readbuffer(sockfd, forces_buffer, nAtom*3)
      CALL readbuffer(sockfd, virial_buffer, 9)
      CALL readbuffer(sockfd, extraLength)
      ALLOCATE (CHARACTER(len=extraLength) :: extra_buffer)
      IF (extraLength /= 0) THEN ! readbuffer(x,y,0) is always an error
         CALL readbuffer(sockfd, extra_buffer, extraLength)
      END IF

      IF (PRESENT(forces)) forces = RESHAPE(forces_buffer, shape=[3, nAtom])
      IF (PRESENT(virial)) virial = RESHAPE(virial_buffer, shape=[3, 3])
      IF (PRESENT(extra)) extra = extra_buffer
   END SUBROUTINE ask_getforce

! **************************************************************************************************
!> \brief ...
!> \param sockfd ...
!> \param subsys ...
! **************************************************************************************************
   SUBROUTINE send_posdata(sockfd, subsys)
      INTEGER, INTENT(IN)                                :: sockfd
      TYPE(cp_subsys_type), POINTER                      :: subsys

      CHARACTER(len=msglength), PARAMETER                :: msg = "POSDATA"

      INTEGER                                            :: i, nAtom, p, xyz
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: particle_buffer
      REAL(kind=dp), DIMENSION(9)                        :: cell_data, icell_data

      i = 0

      CALL writebuffer(sockfd, msg, msglength)

      cell_data = RESHAPE(TRANSPOSE(subsys%cell%hmat), (/9/))
      CALL writebuffer(sockfd, cell_data, 9)

      icell_data = RESHAPE(TRANSPOSE(subsys%cell%h_inv), (/9/))
      CALL writebuffer(sockfd, icell_data, 9)

      nAtom = subsys%particles%n_els
      CALL writebuffer(sockfd, nAtom)

      ALLOCATE (particle_buffer(3*nAtom))
      DO p = 1, nAtom
         DO xyz = 1, 3
            i = i + 1
            particle_buffer(i) = subsys%particles%els(p)%r(xyz)
         END DO
      END DO
      CALL writebuffer(sockfd, particle_buffer, nAtom*3)

   END SUBROUTINE send_posdata

END MODULE ipi_server
