!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE fde_methods
   USE cp_control_types,                ONLY: dft_control_type,&
                                              fde_control_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_to_string
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE xc,                              ONLY: xc_vxc_pw_create
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fde_methods'

   PUBLIC :: &
             fde_build_v_emb, &
             fde_get_env_density, &
             fde_spawn, &
             fde_v_es, &
             fde_write_pw, &
             fde_read_pw, &
             fde_read_results

CONTAINS

! **************************************************************************************************
!> \brief Build the electrostatic potential of a QS calculation
!> \param qs_env ...
!> \param v_es ...
! **************************************************************************************************
   SUBROUTINE fde_v_es(qs_env, v_es)
      TYPE(qs_environment_type), POINTER :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(OUT)  :: v_es

      REAL(KIND=dp)                  :: e_es
      TYPE(pw_c1d_gs_type)           :: rho_tot_g, v_es_g
      TYPE(pw_env_type), POINTER     :: pw_env
      TYPE(pw_poisson_type), POINTER :: poisson_env
      TYPE(pw_pool_type), POINTER    :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER     :: rho

      NULLIFY (pw_env, poisson_env, auxbas_pw_pool, rho)
      CALL get_qs_env(qs_env, pw_env=pw_env, rho=rho)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)

      CALL auxbas_pw_pool%create_pw(rho_tot_g)
      CALL pw_zero(rho_tot_g)
      CALL calc_rho_tot_gspace(rho_tot_g, qs_env, rho)

      CALL auxbas_pw_pool%create_pw(v_es_g)
      CALL pw_poisson_solve(poisson_env, rho_tot_g, e_es, v_es_g)

      CALL auxbas_pw_pool%create_pw(v_es)
      CALL pw_transfer(v_es_g, v_es)

      CALL auxbas_pw_pool%give_back_pw(rho_tot_g)
      CALL auxbas_pw_pool%give_back_pw(v_es_g)

   END SUBROUTINE fde_v_es

! **************************************************************************************************
!> \brief Copy the electron density from a QS calculation into a new grid
!> \param qs_env ...
!> \param rho ...
! **************************************************************************************************
   SUBROUTINE fde_get_env_density(qs_env, rho)
      TYPE(qs_environment_type), POINTER :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(OUT)  :: rho

      INTEGER                                     :: ispin, nspins
      TYPE(dft_control_type), POINTER             :: dft_control
      TYPE(pw_env_type), POINTER                  :: pw_env
      TYPE(pw_pool_type), POINTER                 :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER :: rho_r
      TYPE(qs_rho_type), POINTER                  :: rho_qs

      NULLIFY (dft_control, pw_env, auxbas_pw_pool, rho_r, rho_qs)
      CALL get_qs_env(qs_env, dft_control=dft_control, pw_env=pw_env, rho=rho_qs)
      nspins = dft_control%nspins

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL auxbas_pw_pool%create_pw(rho)
      CALL pw_zero(rho)

      CALL qs_rho_get(rho_qs, rho_r=rho_r)
      DO ispin = 1, nspins
         CALL pw_axpy(rho_r(ispin), rho, 1.0_dp)
      END DO

   END SUBROUTINE fde_get_env_density

! **************************************************************************************************
!> \brief Evaluate an XC+KEDF potential and energy for a density
!> \param ks_env ...
!> \param pw_env ...
!> \param xc_section ...
!> \param rho ...
!> \param v_xc ...
!> \param e_xc ...
! **************************************************************************************************
   SUBROUTINE fde_vxc_on_grid(ks_env, pw_env, xc_section, rho, v_xc, e_xc)
      TYPE(qs_ks_env_type), POINTER     :: ks_env
      TYPE(pw_env_type), POINTER        :: pw_env
      TYPE(section_vals_type), POINTER  :: xc_section
      TYPE(pw_r3d_rs_type), INTENT(IN)  :: rho
      TYPE(pw_r3d_rs_type), INTENT(OUT) :: v_xc
      REAL(KIND=dp), INTENT(OUT)        :: e_xc

      REAL(KIND=dp), DIMENSION(3, 3)              :: virial_xc
      TYPE(pw_c1d_gs_type)                        :: rho_g_single
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER :: rho_g
      TYPE(pw_pool_type), POINTER                 :: auxbas_pw_pool, xc_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER :: rho_r, tau, vxc_rho, vxc_tau
      TYPE(pw_r3d_rs_type), POINTER               :: weights

      NULLIFY (auxbas_pw_pool, xc_pw_pool, rho_r, rho_g, tau, vxc_rho, vxc_tau, weights)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, xc_pw_pool=xc_pw_pool)
      CALL get_ks_env(ks_env, xcint_weights=weights)

      ALLOCATE (rho_r(1), rho_g(1))
      rho_r(1) = rho
      CALL auxbas_pw_pool%create_pw(rho_g_single)
      CALL pw_transfer(rho, rho_g_single)
      rho_g(1) = rho_g_single

      CALL xc_vxc_pw_create(vxc_rho=vxc_rho, &
                            vxc_tau=vxc_tau, &
                            exc=e_xc, &
                            rho_r=rho_r, &
                            rho_g=rho_g, &
                            tau=tau, &
                            xc_section=xc_section, &
                            weights=weights, &
                            pw_pool=xc_pw_pool, &
                            compute_virial=.FALSE., &
                            virial_xc=virial_xc)

      CALL auxbas_pw_pool%create_pw(v_xc)
      CALL pw_copy(vxc_rho(1), v_xc)

      CALL xc_pw_pool%give_back_pw(vxc_rho(1))
      DEALLOCATE (vxc_rho)
      IF (ASSOCIATED(vxc_tau)) THEN
         CALL xc_pw_pool%give_back_pw(vxc_tau(1))
         DEALLOCATE (vxc_tau)
      END IF
      CALL auxbas_pw_pool%give_back_pw(rho_g_single)
      DEALLOCATE (rho_r, rho_g)

   END SUBROUTINE fde_vxc_on_grid

! **************************************************************************************************
!> \brief Build the embedding potential and evaluate nonadditive terms of a QS calculation.
!> \param qs_env ...
!> \param rho ...
!> \param v_emb ...
!> \param e_xc_na ...
!> \param e_kin_na ...
! **************************************************************************************************
   SUBROUTINE fde_build_v_emb(qs_env, rho, v_emb, e_xc_na, e_kin_na)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho
      TYPE(pw_r3d_rs_type), INTENT(OUT)                  :: v_emb
      REAL(KIND=dp), INTENT(OUT)                         :: e_xc_na, e_kin_na

      REAL(KIND=dp)                                      :: exc_a, exc_ab, exc_b
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: rho_ab, rho_b_r, v_ab, v_a, v_tmp
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(section_vals_type), POINTER                   :: input, kg_xc_section, xc_section

      NULLIFY (dft_control, pw_env, auxbas_pw_pool, ks_env, input, xc_section, kg_xc_section)
      CALL get_qs_env(qs_env, dft_control=dft_control, pw_env=pw_env, ks_env=ks_env, input=input)
      CPASSERT(dft_control%nspins == 1)
      CPASSERT(.NOT. dft_control%qs_control%gapw)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      CALL fde_get_env_density(qs_env, rho_b_r)
      CALL auxbas_pw_pool%create_pw(rho_ab)
      CALL pw_copy(rho_b_r, rho_ab)
      CALL pw_axpy(rho, rho_ab, 1.0_dp)

      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      kg_xc_section => section_vals_get_subs_vals(input, "DFT%KG_METHOD%XC")

      CALL fde_v_es(qs_env, v_emb)

      CALL fde_vxc_on_grid(ks_env, pw_env, xc_section, rho_ab, v_ab, exc_ab)
      CALL fde_vxc_on_grid(ks_env, pw_env, xc_section, rho, v_a, exc_a)
      CALL fde_vxc_on_grid(ks_env, pw_env, xc_section, rho_b_r, v_tmp, exc_b) ! energy only
      CALL auxbas_pw_pool%give_back_pw(v_tmp)
      CALL pw_axpy(v_ab, v_emb, 1.0_dp)
      CALL pw_axpy(v_a, v_emb, -1.0_dp)
      CALL auxbas_pw_pool%give_back_pw(v_ab)
      CALL auxbas_pw_pool%give_back_pw(v_a)
      e_xc_na = exc_ab - exc_a - exc_b

      CALL fde_vxc_on_grid(ks_env, pw_env, kg_xc_section, rho_ab, v_ab, exc_ab)
      CALL fde_vxc_on_grid(ks_env, pw_env, kg_xc_section, rho, v_a, exc_a)
      CALL fde_vxc_on_grid(ks_env, pw_env, kg_xc_section, rho_b_r, v_tmp, exc_b) ! energy only
      CALL auxbas_pw_pool%give_back_pw(v_tmp)
      CALL pw_axpy(v_ab, v_emb, 1.0_dp)
      CALL pw_axpy(v_a, v_emb, -1.0_dp)
      CALL auxbas_pw_pool%give_back_pw(v_ab)
      CALL auxbas_pw_pool%give_back_pw(v_a)
      e_kin_na = exc_ab - exc_a - exc_b

      CALL auxbas_pw_pool%give_back_pw(rho_ab)
      CALL auxbas_pw_pool%give_back_pw(rho_b_r)

   END SUBROUTINE fde_build_v_emb

! **************************************************************************************************
!> \brief Spawns the subsystem solver. Only the source rank launches the process; the other ranks
!>        wait for it on the broadcasts below
!> \param fde_control ...
!> \param input_file ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE fde_spawn(fde_control, input_file, para_env)
      TYPE(fde_control_type), POINTER    :: fde_control
      CHARACTER(len=*), INTENT(IN)       :: input_file
      TYPE(mp_para_env_type), INTENT(IN) :: para_env

      CHARACTER(len=:), ALLOCATABLE :: command
      INTEGER                       :: cmdstat, exitstat

      IF (LEN_TRIM(fde_control%workdir) > 0) THEN
         command = "cd "//TRIM(fde_control%workdir)//" && "// &
                   TRIM(fde_control%command)//" "//TRIM(input_file)
      ELSE
         command = TRIM(fde_control%command)//" "//TRIM(input_file)
      END IF

      exitstat = 0
      cmdstat = 0
      IF (para_env%is_source()) &
         CALL EXECUTE_COMMAND_LINE(command, wait=.TRUE., exitstat=exitstat, cmdstat=cmdstat)

      CALL para_env%bcast(cmdstat)
      CALL para_env%bcast(exitstat)

      IF (cmdstat /= 0) &
         CPABORT("FDE: failed to launch command: " // command)
      IF (exitstat /= 0) &
         CPABORT("FDE: command returned a non-zero exit status: " // TRIM(cp_to_string(exitstat)))

   END SUBROUTINE fde_spawn

! **************************************************************************************************
!> \brief Total number of grid points
!> \param pw ...
!> \return ...
! **************************************************************************************************
   FUNCTION fde_pw_npoints(pw) RESULT(npoints)
      TYPE(pw_r3d_rs_type), INTENT(IN) :: pw
      INTEGER                          :: npoints

      npoints = pw%pw_grid%npts(1) * pw%pw_grid%npts(2) * pw%pw_grid%npts(3)

   END FUNCTION fde_pw_npoints

! **************************************************************************************************
!> \brief Write grid points to disk in the Molcas format, optionally with weights and values.
!>        write_value=.TRUE.: x y z weight value
!>        write_value=.FALSE.: x y z
!> \param pw ...
!> \param filename ...
!> \param write_value ...
! **************************************************************************************************
   SUBROUTINE fde_write_pw(pw, filename, write_value)
      TYPE(pw_r3d_rs_type), INTENT(IN) :: pw
      CHARACTER(len=*), INTENT(IN)     :: filename
      LOGICAL, INTENT(IN)              :: write_value

      CHARACTER(len=*), PARAMETER :: routineN = 'fde_write_pw'

      INTEGER                                  :: dest, handle, I1, I2, I3, ip, L1, L2, L3, mepos, &
                                                  num_pe, source, tag, U1, U2, U3, unit_nr
      REAL(KIND=dp)                            :: dvol, x, y, z
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: buf
      TYPE(mp_comm_type)                       :: group

      CALL timeset(routineN, handle)

      L1 = pw%pw_grid%bounds(1, 1)
      U1 = pw%pw_grid%bounds(2, 1)
      L2 = pw%pw_grid%bounds(1, 2)
      U2 = pw%pw_grid%bounds(2, 2)
      L3 = pw%pw_grid%bounds(1, 3)
      U3 = pw%pw_grid%bounds(2, 3)
      dvol = pw%pw_grid%dvol

      group = pw%pw_grid%para%group
      mepos = pw%pw_grid%para%group%mepos
      num_pe = pw%pw_grid%para%group%num_pe
      tag = 1
      dest = 0

      unit_nr = -1
      IF (mepos == dest) THEN
         CALL open_file(file_name=TRIM(filename), file_status="REPLACE", file_form="FORMATTED", &
                        file_action="WRITE", unit_number=unit_nr)
         WRITE (unit_nr, '(I10)') fde_pw_npoints(pw)
      END IF

      ALLOCATE (buf(L3:U3))

      DO I1 = L1, U1
         DO I2 = L2, U2
            IF (pw%pw_grid%para%mode /= PW_MODE_LOCAL) THEN
               source = dest
               DO ip = 0, num_pe - 1
                  IF (pw%pw_grid%para%bo(1, 1, ip, 1) <= I1 - L1 + 1 &
                       .AND. pw%pw_grid%para%bo(2, 1, ip, 1) >= I1 - L1 + 1 &
                       .AND. pw%pw_grid%para%bo(1, 2, ip, 1) <= I2 - L2 + 1 &
                       .AND. pw%pw_grid%para%bo(2, 2, ip, 1) >= I2 - L2 + 1) THEN
                     source = ip
                  END IF
               END DO
            ELSE
               source = dest
            END IF

            IF (source == dest) THEN
               IF (mepos == source) buf(:) = pw%array(I1, I2, :)
            ELSE
               IF (mepos == source) THEN
                  buf(:) = pw%array(I1, I2, :)
                  CALL group%send(buf, dest, tag)
               END IF
               IF (mepos == dest) CALL group%recv(buf, source, tag)
            END IF

            IF (mepos == dest) THEN
               DO I3 = L3, U3
                  x = pw%pw_grid%dh(1, 1)*(I1 - L1) &
                       + pw%pw_grid%dh(2, 1)*(I2 - L2) &
                       + pw%pw_grid%dh(3, 1)*(I3 - L3)
                  y = pw%pw_grid%dh(1, 2)*(I1 - L1) &
                       + pw%pw_grid%dh(2, 2)*(I2 - L2) &
                       + pw%pw_grid%dh(3, 2)*(I3 - L3)
                  z = pw%pw_grid%dh(1, 3)*(I1 - L1) &
                       + pw%pw_grid%dh(2, 3)*(I2 - L2) &
                       + pw%pw_grid%dh(3, 3)*(I3 - L3)
                  IF (write_value) THEN
                     WRITE (unit_nr, '(3F16.8,2ES20.10E3)') x, y, z, dvol, buf(I3)
                  ELSE
                     WRITE (unit_nr, '(3F16.8)') x, y, z
                  END IF
               END DO
            END IF

            CALL group%sync()
         END DO
      END DO

      DEALLOCATE (buf)

      IF (mepos == dest) CALL close_file(unit_number=unit_nr)

      CALL timestop(handle)

   END SUBROUTINE fde_write_pw

! **************************************************************************************************
!> \brief Read the density in x y z weight value format
!> \param pw ...
!> \param filename ...
! **************************************************************************************************
   SUBROUTINE fde_read_pw(pw, filename)
      TYPE(pw_r3d_rs_type), INTENT(INOUT) :: pw
      CHARACTER(len=*), INTENT(IN)        :: filename

      CHARACTER(len=*), PARAMETER :: routineN = 'fde_read_pw'

      INTEGER                                  :: dest, handle, I1, I2, I3, ip, L1, L2, L3, mepos, &
                                                  npoints, num_pe, source, tag, U1, U2, U3, unit_nr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: buf
      TYPE(mp_comm_type)                       :: group

      CALL timeset(routineN, handle)

      L1 = pw%pw_grid%bounds(1, 1)
      U1 = pw%pw_grid%bounds(2, 1)
      L2 = pw%pw_grid%bounds(1, 2)
      U2 = pw%pw_grid%bounds(2, 2)
      L3 = pw%pw_grid%bounds(1, 3)
      U3 = pw%pw_grid%bounds(2, 3)

      group = pw%pw_grid%para%group
      mepos = pw%pw_grid%para%group%mepos
      num_pe = pw%pw_grid%para%group%num_pe
      tag = 1
      dest = 0

      unit_nr = -1
      IF (mepos == dest) THEN
         CALL open_file(file_name=TRIM(filename), file_status="OLD", file_form="FORMATTED", &
                        file_action="READ", unit_number=unit_nr)
         READ (unit_nr, *) npoints
         IF (npoints /= fde_pw_npoints(pw)) &
              CPABORT("Number of grid points in file does not match the pw grid")
      END IF

      ALLOCATE (buf(L3:U3))

      DO I1 = L1, U1
         DO I2 = L2, U2
            IF (pw%pw_grid%para%mode /= PW_MODE_LOCAL) THEN
               source = dest
               DO ip = 0, num_pe - 1
                  IF (pw%pw_grid%para%bo(1, 1, ip, 1) <= I1 - L1 + 1 &
                       .AND. pw%pw_grid%para%bo(2, 1, ip, 1) >= I1 - L1 + 1 &
                       .AND. pw%pw_grid%para%bo(1, 2, ip, 1) <= I2 - L2 + 1 &
                       .AND. pw%pw_grid%para%bo(2, 2, ip, 1) >= I2 - L2 + 1) THEN
                     source = ip
                  END IF
               END DO
            ELSE
               source = dest
            END IF

            IF (mepos == dest) THEN
               DO I3 = L3, U3
                  READ (unit_nr, *) buf(I3)
               END DO
            END IF

            IF (source == dest) THEN
               IF (mepos == source) pw%array(I1, I2, :) = buf(:)
            ELSE
               IF (mepos == dest) CALL group%send(buf, source, tag)
               IF (mepos == source) THEN
                  CALL group%recv(buf, dest, tag)
                  pw%array(I1, I2, :) = buf(:)
               END IF
            END IF

            CALL group%sync()
         END DO
      END DO

      DEALLOCATE (buf)

      IF (mepos == dest) CALL close_file(unit_number=unit_nr)

      CALL timestop(handle)

   END SUBROUTINE fde_read_pw

! **************************************************************************************************
!> \brief Read scalar results returned by the subsystem solver, one value per line, in this order:
!>        E_sub : embedded subsystem total energy
!>        E_emb : embedding energy
!>        E_nuc : nuclear repulsion energy
!>        n_e   : number of active-subsystem electrons
!>        n_a   : number of nuclei of subsystem A
!>        Followed by n_a lines of
!>        Z_I x_I y_I z_I (charge and coordinates in the CP2K frame)
!> \param filename ...
!> \param para_env ...
!> \param e_sub ...
!> \param e_emb ...
!> \param e_nuc ...
!> \param n_e ...
!> \param n_a ...
!> \param z_a ...
!> \param r_a ...
! **************************************************************************************************
   SUBROUTINE fde_read_results(filename, para_env, e_sub, e_emb, e_nuc, n_e, n_a, z_a, r_a)
      CHARACTER(len=*), INTENT(IN)                             :: filename
      TYPE(mp_para_env_type), INTENT(IN)                       :: para_env
      REAL(KIND=dp), INTENT(OUT)                               :: e_sub, e_emb, e_nuc, n_e
      INTEGER, INTENT(OUT)                                     :: n_a
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: z_a
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: r_a

      INTEGER                                            :: iatom, unit_nr
      REAL(KIND=dp), DIMENSION(4)                        :: vals

      vals = 0.0_dp
      n_a = 0
      IF (para_env%is_source()) THEN
         CALL open_file(file_name=TRIM(filename), file_status="OLD", file_form="FORMATTED", &
                        file_action="READ", unit_number=unit_nr)
         READ (unit_nr, *) vals(1)
         READ (unit_nr, *) vals(2)
         READ (unit_nr, *) vals(3)
         READ (unit_nr, *) vals(4)
         READ (unit_nr, *) n_a
      END IF
      CALL para_env%bcast(vals)
      CALL para_env%bcast(n_a)

      ALLOCATE (z_a(n_a), r_a(3, n_a))
      z_a = 0.0_dp
      r_a = 0.0_dp
      IF (para_env%is_source()) THEN
         DO iatom = 1, n_a
            READ (unit_nr, *) z_a(iatom), r_a(1, iatom), r_a(2, iatom), r_a(3, iatom)
         END DO
         CALL close_file(unit_number=unit_nr)
      END IF
      CALL para_env%bcast(z_a)
      CALL para_env%bcast(r_a)

      e_sub = vals(1)
      e_emb = vals(2)
      e_nuc = vals(3)
      n_e = vals(4)

   END SUBROUTINE fde_read_results

END MODULE fde_methods
