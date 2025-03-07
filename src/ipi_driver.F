!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Driver mode - To communicate with i-PI Python wrapper
!> \par History
!>      none
!> \author Michele Ceriotti 03.2012
! **************************************************************************************************
MODULE ipi_driver
   USE ISO_C_BINDING, ONLY: C_CHAR, &
                            C_DOUBLE, &
                            C_INT, &
                            C_LOC, &
                            C_NULL_CHAR, &
                            C_PTR
   USE bibliography, ONLY: Ceriotti2014, &
                           Kapil2016, &
                           cite_reference
   USE cell_methods, ONLY: cell_create, &
                           init_cell
   USE cell_types, ONLY: cell_release, &
                         cell_type
   USE cp_external_control, ONLY: external_control
   USE cp_log_handling, ONLY: cp_logger_get_default_io_unit
   USE cp_subsys_types, ONLY: cp_subsys_get, &
                              cp_subsys_set, &
                              cp_subsys_type
   USE force_env_methods, ONLY: force_env_calc_energy_force
   USE force_env_types, ONLY: force_env_get, &
                              force_env_type
   USE global_types, ONLY: global_environment_type
   USE input_section_types, ONLY: section_vals_get_subs_vals, &
                                  section_vals_type, &
                                  section_vals_val_get
   USE kinds, ONLY: default_path_length, &
                    default_string_length, &
                    dp, &
                    int_4
   USE message_passing, ONLY: mp_para_env_type, &
                              mp_request_type, &
                              mp_testany
#ifndef __NO_SOCKETS
   USE sockets_interface, ONLY: writebuffer, &
                                readbuffer, &
                                open_connect_socket, &
                                uwait
#endif
   USE virial_types, ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ipi_driver'

   PUBLIC :: run_driver

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param force_env ...
!> \param globenv ...
!> \par History
!>       12.2013 included in repository
!> \author Ceriotti
! **************************************************************************************************

   SUBROUTINE run_driver(force_env, globenv)
      TYPE(force_env_type), POINTER            :: force_env
      TYPE(global_environment_type), POINTER   :: globenv

      CHARACTER(len=*), PARAMETER :: routineN = 'run_driver'

#ifdef __NO_SOCKETS
      INTEGER                                  :: handle
      CALL timeset(routineN, handle)
      CPABORT("CP2K was compiled with the __NO_SOCKETS option!")
      MARK_USED(globenv)
      MARK_USED(force_env)
#else
      INTEGER, PARAMETER                       :: MSGLEN = 12

      CHARACTER(len=default_path_length)       :: c_hostname, drv_hostname
      CHARACTER(LEN=default_string_length)     :: header
      INTEGER                                  :: drv_port, handle, i_drv_unix, &
                                                  idir, ii, inet, ip, iwait, &
                                                  nat, output_unit, socket
      TYPE(mp_request_type), DIMENSION(2) ::                                            wait_req
      INTEGER(KIND=int_4), POINTER             :: wait_msg(:)
      LOGICAL                                  :: drv_unix, fwait, hasdata, &
                                                  ionode, should_stop
      REAL(KIND=dp)                            :: cellh(3, 3), cellih(3, 3), &
                                                  mxmat(9), pot, vir(3, 3)
      REAL(KIND=dp), ALLOCATABLE               :: combuf(:)
      TYPE(cell_type), POINTER                 :: cpcell
      TYPE(mp_para_env_type), POINTER          :: para_env
      TYPE(cp_subsys_type), POINTER            :: subsys
      TYPE(section_vals_type), POINTER         :: drv_section, motion_section
      TYPE(virial_type), POINTER               :: virial
      REAL(KIND=dp)                            :: sleeptime

      CALL timeset(routineN, handle)

      CALL cite_reference(Ceriotti2014)
      CALL cite_reference(Kapil2016)

! server address parsing
! buffers and temporaries for communication
! access cp2k structures

      CPASSERT(ASSOCIATED(force_env))
      CALL force_env_get(force_env, para_env=para_env)

      hasdata = .FALSE.
      ionode = para_env%is_source()

      output_unit = cp_logger_get_default_io_unit()

      ! reads driver parameters from input
      motion_section => section_vals_get_subs_vals(force_env%root_section, "MOTION")
      drv_section => section_vals_get_subs_vals(motion_section, "DRIVER")

      CALL section_vals_val_get(drv_section, "HOST", c_val=drv_hostname)
      CALL section_vals_val_get(drv_section, "PORT", i_val=drv_port)
      CALL section_vals_val_get(drv_section, "UNIX", l_val=drv_unix)
      CALL section_vals_val_get(drv_section, "SLEEP_TIME", r_val=sleeptime)
      CPASSERT(sleeptime >= 0)

      ! opens the socket
      socket = 0
      inet = 1
      i_drv_unix = 1 ! a bit convoluted. socket.c uses a different convention...
      IF (drv_unix) i_drv_unix = 0
      IF (output_unit > 0) THEN
         WRITE (output_unit, *) "@ i-PI DRIVER BEING LOADED"
         WRITE (output_unit, *) "@ INPUT DATA: ", TRIM(drv_hostname), drv_port, drv_unix
      END IF

      c_hostname = TRIM(drv_hostname)//C_NULL_CHAR
      IF (ionode) CALL open_connect_socket(socket, i_drv_unix, drv_port, c_hostname)

      NULLIFY (wait_msg)
      ALLOCATE (wait_msg(1))
      !now we have a socket, so we can initialize the CP2K environments.
      NULLIFY (cpcell)
      CALL cell_create(cpcell)
      driver_loop: DO
         ! do communication on master node only...
         header = ""

         CALL para_env%sync()

         ! non-blocking sync to avoid useless CPU consumption
         IF (ionode) THEN
            CALL readbuffer(socket, header, MSGLEN)
            wait_msg = 0
            DO iwait = 0, para_env%num_pe - 1
               IF (iwait /= para_env%source) THEN
                  CALL para_env%send(msg=wait_msg, dest=iwait, tag=666)
               END IF
            END DO
         ELSE
            CALL para_env%irecv(msgout=wait_msg, source=para_env%source, &
                                tag=666, request=wait_req(2))
            CALL mp_testany(wait_req(2:), flag=fwait)
            DO WHILE (.NOT. fwait)
               CALL mp_testany(wait_req(2:), flag=fwait)
               CALL uwait(sleeptime)
            END DO
         END IF

         CALL para_env%sync()

         CALL para_env%bcast(header)

         IF (output_unit > 0) WRITE (output_unit, *) " @ DRIVER MODE: Message from server: ", TRIM(header)
         IF (TRIM(header) == "STATUS") THEN

            CALL para_env%sync()
            IF (ionode) THEN ! does not  need init (well, maybe it should, just to check atom numbers and the like... )
               IF (hasdata) THEN
                  CALL writebuffer(socket, "HAVEDATA    ", MSGLEN)
               ELSE
                  CALL writebuffer(socket, "READY       ", MSGLEN)
               END IF
            END IF
            CALL para_env%sync()
         ELSE IF (TRIM(header) == "POSDATA") THEN
            IF (ionode) THEN
               CALL readbuffer(socket, mxmat, 9)
               cellh = RESHAPE(mxmat, (/3, 3/))
               CALL readbuffer(socket, mxmat, 9)
               cellih = RESHAPE(mxmat, (/3, 3/))
               CALL readbuffer(socket, nat)
               cellh = TRANSPOSE(cellh)
               cellih = TRANSPOSE(cellih)
            END IF
            CALL para_env%bcast(cellh)
            CALL para_env%bcast(cellih)
            CALL para_env%bcast(nat)
            IF (.NOT. ALLOCATED(combuf)) ALLOCATE (combuf(3*nat))
            IF (ionode) CALL readbuffer(socket, combuf, nat*3)
            CALL para_env%bcast(combuf)

            CALL force_env_get(force_env, subsys=subsys)
            IF (nat /= subsys%particles%n_els) &
               CPABORT("@DRIVER MODE: Uh-oh! Particle number mismatch between i-PI and cp2k input!")
            ii = 0
            DO ip = 1, subsys%particles%n_els
               DO idir = 1, 3
                  ii = ii + 1
                  subsys%particles%els(ip)%r(idir) = combuf(ii)
               END DO
            END DO
            CALL init_cell(cpcell, hmat=cellh)
            CALL cp_subsys_set(subsys, cell=cpcell)

            CALL force_env_calc_energy_force(force_env, calc_force=.TRUE.)

            IF (output_unit > 0) WRITE (output_unit, *) " @ DRIVER MODE: Received positions "

            combuf = 0
            ii = 0
            DO ip = 1, subsys%particles%n_els
               DO idir = 1, 3
                  ii = ii + 1
                  combuf(ii) = subsys%particles%els(ip)%f(idir)
               END DO
            END DO
            CALL force_env_get(force_env, potential_energy=pot)
            CALL force_env_get(force_env, cell=cpcell)
            CALL cp_subsys_get(subsys, virial=virial)
            vir = TRANSPOSE(virial%pv_virial)

            CALL external_control(should_stop, "IPI", globenv=globenv)
            IF (should_stop) EXIT

            hasdata = .TRUE.
         ELSE IF (TRIM(header) == "GETFORCE") THEN
            IF (output_unit > 0) WRITE (output_unit, *) " @ DRIVER MODE: Returning v,forces,stress "
            IF (ionode) THEN
               CALL writebuffer(socket, "FORCEREADY  ", MSGLEN)
               CALL writebuffer(socket, pot)
               CALL writebuffer(socket, nat)
               CALL writebuffer(socket, combuf, 3*nat)
               CALL writebuffer(socket, RESHAPE(vir, (/9/)), 9)

               ! i-pi can also receive an arbitrary string, that will be printed out to the "extra"
               ! trajectory file. this is useful if you want to return additional information, e.g.
               ! atomic charges, wannier centres, etc. one must return the number of characters, then
               ! the string. here we just send back zero characters.
               nat = 0
               CALL writebuffer(socket, nat) ! writes out zero for the length of the "extra" field (not implemented yet!)
            END IF
            hasdata = .FALSE.
         ELSE
            IF (output_unit > 0) WRITE (output_unit, *) " @DRIVER MODE:  Socket disconnected, time to exit."
            EXIT
         END IF
      END DO driver_loop

      ! clean up
      CALL cell_release(cpcell)
      DEALLOCATE (wait_msg)
#endif

      CALL timestop(handle)

   END SUBROUTINE run_driver
END MODULE ipi_driver
