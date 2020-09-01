!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

!***************************************************************************************************
!> \brief Interface to the SIRIUS Library
!> \par History
!>      07.2018 initial create
!> \author JHU
!***************************************************************************************************
#if defined(__SIRIUS)
MODULE sirius_interface
   USE ISO_C_BINDING,                   ONLY: C_DOUBLE,&
                                              C_INT,&
                                              C_NULL_PTR,&
                                              C_PTR
   USE atom_kind_orbitals,              ONLY: calculate_atomic_orbitals,&
                                              gth_potential_conversion
   USE atom_types,                      ONLY: atom_gthpot_type
   USE atom_upf,                        ONLY: atom_upfpot_type
   USE atom_utils,                      ONLY: atom_local_potential
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              real_to_scaled
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE external_potential_types,        ONLY: gth_potential_type
   USE input_constants,                 ONLY: do_gapw_log
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: fourpi,&
                                              gamma1
   USE particle_types,                  ONLY: particle_type
   USE pwdft_environment_types,         ONLY: pwdft_energy_type,&
                                              pwdft_env_get,&
                                              pwdft_env_set,&
                                              pwdft_environment_type
   USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                              create_grid_atom,&
                                              deallocate_grid_atom,&
                                              grid_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE sirius,                          ONLY: &
        sirius_add_atom, sirius_add_atom_type, sirius_add_atom_type_radial_function, &
        sirius_create_context, sirius_create_ground_state, sirius_create_kset_from_grid, &
        sirius_dump_runtime_setup, sirius_finalize, sirius_find_ground_state, &
        sirius_get_band_energies, sirius_get_band_occupancies, sirius_get_energy, &
        sirius_get_forces, sirius_get_kpoint_properties, sirius_get_num_bands, &
        sirius_get_num_kpoints, sirius_get_num_spin_components, sirius_get_stress_tensor, &
        sirius_import_parameters, sirius_initialize, sirius_initialize_context, &
        sirius_option_add_string_to, sirius_option_get_length, sirius_option_get_name_and_type, &
        sirius_option_set_double, sirius_option_set_int, sirius_option_set_logical, &
        sirius_option_set_string, sirius_set_atom_position, sirius_set_atom_type_dion, &
        sirius_set_atom_type_radial_grid, sirius_set_lattice_vectors, sirius_update_ground_state
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

!     *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sirius_interface'

!     *** Public subroutines ***

   PUBLIC :: cp_sirius_init, cp_sirius_finalize
   PUBLIC :: cp_sirius_create_env, cp_sirius_energy_force, cp_sirius_update_context

CONTAINS

!***************************************************************************************************
!> \brief ...
!> \param
!> \par History
!>      07.2018 start the Sirius library
!> \author JHU
! **************************************************************************************************
   SUBROUTINE cp_sirius_init()
      CALL sirius_initialize(.FALSE.)
   END SUBROUTINE cp_sirius_init

!***************************************************************************************************
!> \brief ...
!> \param
!> \par History
!>      07.2018 stop the Sirius library
!> \author JHU
! **************************************************************************************************
   SUBROUTINE cp_sirius_finalize()
      CALL sirius_finalize(.FALSE., .FALSE., .FALSE.)
   END SUBROUTINE cp_sirius_finalize

!***************************************************************************************************
!> \brief ...
!> \param pwdft_env ...
!> \param
!> \par History
!>      07.2018 Create the Sirius environment
!> \author JHU
! **************************************************************************************************
   SUBROUTINE cp_sirius_create_env(pwdft_env)
      TYPE(pwdft_environment_type), POINTER              :: pwdft_env

      CHARACTER(len=2)                                   :: element_symbol
      CHARACTER(len=64)                                  :: section_name
      CHARACTER(len=default_string_length)               :: label
      INTEGER                                            :: i, iatom, ibeta, ifun, ikind, iwf, j, l, &
                                                            n, natom, nbeta, nkind, nmesh, &
                                                            num_mag_dims, sirius_mpi_comm
      INTEGER(KIND=C_INT), DIMENSION(3)                  :: k_grid, k_shift
      INTEGER, DIMENSION(:), POINTER                     :: kk
      LOGICAL                                            :: up, use_ref_cell
      LOGICAL(4)                                         :: use_symmetry
      REAL(KIND=C_DOUBLE), ALLOCATABLE, DIMENSION(:)     :: fun
      REAL(KIND=C_DOUBLE), ALLOCATABLE, DIMENSION(:, :)  :: dion
      REAL(KIND=C_DOUBLE), DIMENSION(3)                  :: a1, a2, a3, v1, v2
      REAL(KIND=dp)                                      :: al, angle1, angle2, cval, focc, &
                                                            magnetization, mass, pf, rl, zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: beta, corden, ef, fe, locpot, rc, rp
      REAL(KIND=dp), DIMENSION(3)                        :: vr, vs
      REAL(KIND=dp), DIMENSION(:), POINTER               :: density
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: wavefunction, wfninfo
      TYPE(atom_gthpot_type), POINTER                    :: gth_atompot
      TYPE(atom_upfpot_type), POINTER                    :: upf_pot
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(C_PTR)                                        :: gs_handler = C_NULL_PTR, &
                                                            ks_handler = C_NULL_PTR, &
                                                            sctx = C_NULL_PTR
      TYPE(cell_type), POINTER                           :: my_cell, my_cell_ref
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(grid_atom_type), POINTER                      :: atom_grid
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: qs_subsys
      TYPE(section_vals_type), POINTER                   :: pwdft_section, pwdft_sub_section, &
                                                            xc_fun, xc_section

      CPASSERT(ASSOCIATED(pwdft_env))
! create context of simulation
      CALL pwdft_env_get(pwdft_env, para_env=para_env)
      sirius_mpi_comm = para_env%group
      CALL sirius_create_context(sirius_mpi_comm, sctx)

!     the "fun" starts.

      CALL pwdft_env_get(pwdft_env=pwdft_env, pwdft_input=pwdft_section, xc_input=xc_section)

      ! cp2k should *have* a function that return all xc_functionals. Doing
      ! manually is prone to errors

      IF (ASSOCIATED(xc_section)) THEN
         !     LIBXC is the only repeatable functional section - for each we need
         !     NOT the single values, but the whole section_vals_type independently
         ifun = 0
         DO
            ifun = ifun + 1
            xc_fun => section_vals_get_subs_vals2(xc_section, i_section=ifun)
            IF (.NOT. ASSOCIATED(xc_fun)) EXIT
            IF (TRIM(xc_fun%section%name) == "LIBXC") THEN
               CALL section_vals_get(xc_fun, n_repetition=n)
               DO i = 1, n
                  CALL section_vals_val_get(xc_fun, "FUNCTIONAL", i_rep_section=i, c_val=section_name)
                  CALL sirius_option_add_string_to(sctx, 'parameters', 'xc_functionals', section_name)
               END DO
            ENDIF
         END DO
      ENDIF

!     import control section
      pwdft_sub_section => section_vals_get_subs_vals(pwdft_section, "control")
      IF (ASSOCIATED(pwdft_sub_section)) THEN
         CALL cp_sirius_fill_in_section(sctx, pwdft_sub_section, "control")
      ENDIF

!     import parameters section
      pwdft_sub_section => section_vals_get_subs_vals(pwdft_section, "parameters")
      IF (ASSOCIATED(pwdft_sub_section)) THEN
         CALL cp_sirius_fill_in_section(sctx, pwdft_sub_section, "parameters")
         CALL section_vals_val_get(pwdft_sub_section, "ngridk", i_vals=kk)
         k_grid(1) = kk(1)
         k_grid(2) = kk(2)
         k_grid(3) = kk(3)

         CALL section_vals_val_get(pwdft_sub_section, "shiftk", i_vals=kk)
         k_shift(1) = kk(1)
         k_shift(2) = kk(2)
         k_shift(3) = kk(3)
         CALL section_vals_val_get(pwdft_sub_section, "num_mag_dims", i_val=num_mag_dims)
         CALL section_vals_val_get(pwdft_sub_section, "use_symmetry", l_val=use_symmetry)
      ENDIF

!     import mixer section
      pwdft_sub_section => section_vals_get_subs_vals(pwdft_section, "mixer")
      IF (ASSOCIATED(pwdft_sub_section)) THEN
         CALL cp_sirius_fill_in_section(sctx, pwdft_sub_section, "mixer")
      ENDIF

      !     import solver section
      pwdft_sub_section => section_vals_get_subs_vals(pwdft_section, "iterative_solver")
      IF (ASSOCIATED(pwdft_sub_section)) THEN
         CALL cp_sirius_fill_in_section(sctx, pwdft_sub_section, "iterative_solver")
      ENDIF

      CALL sirius_dump_runtime_setup(sctx, "runtime.json")
      CALL sirius_import_parameters(sctx)

! lattice vectors of the unit cell should be in [a.u.] (length is in [a.u.])
      CALL pwdft_env_get(pwdft_env=pwdft_env, qs_subsys=qs_subsys)
      CALL qs_subsys_get(qs_subsys, cell=my_cell, cell_ref=my_cell_ref, use_ref_cell=use_ref_cell)
      a1(:) = my_cell%hmat(:, 1)
      a2(:) = my_cell%hmat(:, 2)
      a3(:) = my_cell%hmat(:, 3)
      CALL sirius_set_lattice_vectors(sctx, a1(1), a2(1), a3(1))

! set up the atomic type definitions
      CALL qs_subsys_get(qs_subsys, &
                         atomic_kind_set=atomic_kind_set, &
                         qs_kind_set=qs_kind_set, &
                         particle_set=particle_set)
      nkind = SIZE(atomic_kind_set)
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), &
                              name=label, element_symbol=element_symbol, mass=mass)
         CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
         NULLIFY (upf_pot, gth_potential)
         CALL get_qs_kind(qs_kind_set(ikind), upf_potential=upf_pot, gth_potential=gth_potential)
         IF (ASSOCIATED(upf_pot)) THEN
            CALL sirius_add_atom_type(sctx, label, fname=upf_pot%filename, &
                                      mass=REAL(mass, KIND=C_DOUBLE))
         ELSEIF (ASSOCIATED(gth_potential)) THEN
!
            NULLIFY (atom_grid)
            CALL allocate_grid_atom(atom_grid)
            nmesh = 929
            atom_grid%nr = nmesh
            CALL create_grid_atom(atom_grid, nmesh, 1, 1, 0, do_gapw_log)
            ALLOCATE (rp(nmesh), fun(nmesh))
            IF (atom_grid%rad(1) < atom_grid%rad(nmesh)) THEN
               up = .TRUE.
            ELSE
               up = .FALSE.
            END IF
            IF (up) THEN
               rp(1:nmesh) = atom_grid%rad(1:nmesh)
            ELSE
               DO i = 1, nmesh
                  rp(i) = atom_grid%rad(nmesh - i + 1)
               END DO
            END IF
! add new atom type
            CALL sirius_add_atom_type(sctx, label, &
                                      zn=NINT(zeff + 0.001d0), &
                                      mass=REAL(mass, KIND=C_DOUBLE), &
                                      spin_orbit=.FALSE.)
!
            ALLOCATE (gth_atompot)
            CALL gth_potential_conversion(gth_potential, gth_atompot)
! set radial grid
            fun(1:nmesh) = rp(1:nmesh)
            CALL sirius_set_atom_type_radial_grid(sctx, label, nmesh, fun(1))
! set beta-projectors
            ALLOCATE (ef(nmesh), beta(nmesh))
            ibeta = 0
            DO l = 0, 3
               IF (gth_atompot%nl(l) == 0) CYCLE
               rl = gth_atompot%rcnl(l)
! we need to multiply by r so that data transferred to sirius are r \beta(r) not beta(r)
               ef(1:nmesh) = EXP(-0.5_dp*rp(1:nmesh)*rp(1:nmesh)/(rl*rl))
               DO i = 1, gth_atompot%nl(l)
                  pf = rl**(l + 0.5_dp*(4._dp*i - 1._dp))
                  j = l + 2*i - 1
                  pf = SQRT(2._dp)/(pf*SQRT(gamma1(j)))
                  beta(:) = pf*rp**(l + 2*i - 2)*ef
                  ibeta = ibeta + 1
                  fun(1:nmesh) = beta(1:nmesh)*rp(1:nmesh)
                  CALL sirius_add_atom_type_radial_function(sctx, label, &
                                                            "beta", fun(1), nmesh, l=l)
               END DO
            END DO
            DEALLOCATE (ef, beta)
            nbeta = ibeta

! nonlocal PP matrix elements
            ALLOCATE (dion(nbeta, nbeta))
            dion = 0.0_dp
            DO l = 0, 3
               IF (gth_atompot%nl(l) == 0) CYCLE
               ibeta = SUM(gth_atompot%nl(0:l - 1)) + 1
               i = ibeta + gth_atompot%nl(l) - 1
               dion(ibeta:i, ibeta:i) = gth_atompot%hnl(1:gth_atompot%nl(l), 1:gth_atompot%nl(l), l)
            END DO
            CALL sirius_set_atom_type_dion(sctx, label, nbeta, dion(1, 1))
            DEALLOCATE (dion)

! set non-linear core correction
            IF (gth_atompot%nlcc) THEN
               ALLOCATE (corden(nmesh), fe(nmesh), rc(nmesh))
               corden(:) = 0.0_dp
               n = gth_atompot%nexp_nlcc
               DO i = 1, n
                  al = gth_atompot%alpha_nlcc(i)
                  rc(:) = rp(:)/al
                  fe(:) = EXP(-0.5_dp*rc(:)*rc(:))
                  DO j = 1, gth_atompot%nct_nlcc(i)
                     cval = gth_atompot%cval_nlcc(j, i)
                     corden(:) = corden(:) + fe(:)*rc(:)**(2*j - 2)*cval
                  END DO
               END DO
               fun(1:nmesh) = corden(1:nmesh)*rp(1:nmesh)
               CALL sirius_add_atom_type_radial_function(sctx, label, "ps_rho_core", &
                                                         fun(1), nmesh)
               DEALLOCATE (corden, fe, rc)
            END IF

! local potential
            ALLOCATE (locpot(nmesh))
            locpot(:) = 0.0_dp
            CALL atom_local_potential(locpot, gth_atompot, rp)
            fun(1:nmesh) = locpot(1:nmesh)
            CALL sirius_add_atom_type_radial_function(sctx, label, "vloc", &
                                                      fun(1), nmesh)
            DEALLOCATE (locpot)
!
            NULLIFY (density, wavefunction, wfninfo)
            CALL calculate_atomic_orbitals(atomic_kind_set(ikind), qs_kind_set(ikind), &
                                           density=density, wavefunction=wavefunction, &
                                           wfninfo=wfninfo, agrid=atom_grid)

! set the atomic radial functions
            DO iwf = 1, SIZE(wavefunction, 2)
               focc = wfninfo(1, iwf)
               l = NINT(wfninfo(2, iwf))
               IF (up) THEN
                  fun(1:nmesh) = wavefunction(1:nmesh, iwf)*rp(i)
               ELSE
                  DO i = 1, nmesh
                     fun(i) = wavefunction(nmesh - i + 1, iwf)*rp(i)
                  END DO
               END IF
               CALL sirius_add_atom_type_radial_function(sctx, &
                                                         label, "ps_atomic_wf", &
                                                         fun(1), nmesh, l=l, occ=REAL(focc, KIND=C_DOUBLE), n=-1)
            END DO

! set total charge density of a free atom (to compute initial rho(r))
            IF (up) THEN
               fun(1:nmesh) = fourpi*density(1:nmesh)*atom_grid%rad(1:nmesh)**2
            ELSE
               DO i = 1, nmesh
                  fun(i) = fourpi*density(nmesh - i + 1)*atom_grid%rad(nmesh - i + 1)**2
               END DO
            END IF
            CALL sirius_add_atom_type_radial_function(sctx, label, "ps_rho_total", &
                                                      fun(1), nmesh)

            IF (ASSOCIATED(density)) DEALLOCATE (density)
            IF (ASSOCIATED(wavefunction)) DEALLOCATE (wavefunction)
            IF (ASSOCIATED(wfninfo)) DEALLOCATE (wfninfo)

            CALL deallocate_grid_atom(atom_grid)
            DEALLOCATE (rp, fun)
            DEALLOCATE (gth_atompot)
!
         ELSE
            CALL cp_abort(__LOCATION__, &
                          'CP2K/SIRIUS: atomic kind needs UPF or GTH potential definition')
         END IF
      END DO

! add atoms to the unit cell
! WARNING: sirius accepts only fractional coordinates;
      natom = SIZE(particle_set)
      DO iatom = 1, natom
         vr(1:3) = particle_set(iatom)%r(1:3)
         CALL real_to_scaled(vs, vr, my_cell)
         atomic_kind => particle_set(iatom)%atomic_kind
         ikind = atomic_kind%kind_number
         CALL get_atomic_kind(atomic_kind, name=label)
         CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff, magnetization=magnetization)
! angle of magnetization might come from input Atom x y z mx my mz
! or as an angle?
! Answer : SIRIUS only accept the magnetization as mx, my, mz
         IF (num_mag_dims .EQ. 3) THEN
            angle1 = 0.0_dp
            angle2 = 0.0_dp
            v1(1) = zeff*magnetization*SIN(angle1)*COS(angle2)
            v1(2) = zeff*magnetization*SIN(angle1)*SIN(angle2)
            v1(3) = zeff*magnetization*COS(angle1)
         ELSE
            v1 = 0._dp
            v1(3) = zeff*magnetization
         ENDIF
         v2(1:3) = vs(1:3)
         CALL sirius_add_atom(sctx, label, v2(1), v1(1))
      ENDDO
! initialize global variables/indices/arrays/etc. of the simulation
      CALL sirius_initialize_context(sctx)

      ! strictly speaking the parameter use_symmetry is initialized at the
      ! beginning but it does no harm to do it that way
      IF (use_symmetry) THEN
         CALL sirius_create_kset_from_grid(sctx, k_grid(1), k_shift(1), use_symmetry=.TRUE., kset_handler=ks_handler)
      ELSE
         CALL sirius_create_kset_from_grid(sctx, k_grid(1), k_shift(1), use_symmetry=.FALSE., kset_handler=ks_handler)
      ENDIF
! create ground-state class
      CALL sirius_create_ground_state(ks_handler, gs_handler)

      CALL pwdft_env_set(pwdft_env, sctx=sctx, gs_handler=gs_handler, ks_handler=ks_handler)
   END SUBROUTINE cp_sirius_create_env

!***************************************************************************************************
!> \brief ...
!> \param pwdft_env ...
!> \param
!> \par History
!>      07.2018 Update the Sirius environment
!> \author JHU
! **************************************************************************************************
   SUBROUTINE cp_sirius_update_context(pwdft_env)
      TYPE(pwdft_environment_type), POINTER              :: pwdft_env

      INTEGER                                            :: iatom, natom
      REAL(KIND=C_DOUBLE), DIMENSION(3)                  :: a1, a2, a3, v2
      REAL(KIND=dp), DIMENSION(3)                        :: vr, vs
      TYPE(C_PTR)                                        :: gs_handler = C_NULL_PTR, &
                                                            sctx = C_NULL_PTR
      TYPE(cell_type), POINTER                           :: my_cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_subsys_type), POINTER                      :: qs_subsys

      CPASSERT(ASSOCIATED(pwdft_env))
      CALL pwdft_env_get(pwdft_env, sctx=sctx, gs_handler=gs_handler)

! get current positions and lattice vectors
      CALL pwdft_env_get(pwdft_env=pwdft_env, qs_subsys=qs_subsys)

! lattice vectors of the unit cell should be in [a.u.] (length is in [a.u.])
      CALL qs_subsys_get(qs_subsys, cell=my_cell)
      a1(:) = my_cell%hmat(:, 1)
      a2(:) = my_cell%hmat(:, 2)
      a3(:) = my_cell%hmat(:, 3)
      CALL sirius_set_lattice_vectors(sctx, a1(1), a2(1), a3(1))

! new atomic positions
      CALL qs_subsys_get(qs_subsys, particle_set=particle_set)
      natom = SIZE(particle_set)
      DO iatom = 1, natom
         vr(1:3) = particle_set(iatom)%r(1:3)
         CALL real_to_scaled(vs, vr, my_cell)
         v2(1:3) = vs(1:3)
         CALL sirius_set_atom_position(sctx, iatom, v2(1))
      ENDDO

! update ground-state class
      CALL sirius_update_ground_state(gs_handler)

      CALL pwdft_env_set(pwdft_env, sctx=sctx, gs_handler=gs_handler)

   END SUBROUTINE cp_sirius_update_context

! **************************************************************************************************
!> \brief ...
!> \param sctx ...
!> \param section ...
!> \param section_name ...
! **************************************************************************************************
   SUBROUTINE cp_sirius_fill_in_section(sctx, section, section_name)
      TYPE(C_PTR), INTENT(INOUT)                         :: sctx
      TYPE(section_vals_type), POINTER                   :: section
      CHARACTER(*), INTENT(in), TARGET                   :: section_name

      CHARACTER(len=256)                                 :: option_name, option_name1
      CHARACTER(len=80)                                  :: str
      CHARACTER(len=80), DIMENSION(:), POINTER           :: tmp
      INTEGER                                            :: ctype, elem, i, ival, j, length, &
                                                            lvalsi(16), number_of_options
      INTEGER, DIMENSION(:), POINTER                     :: ivals
      LOGICAL                                            :: found
      LOGICAL(4)                                         :: lval
      LOGICAL(4), DIMENSION(:), POINTER                  :: lvals
      REAL(kind=dp)                                      :: rval
      REAL(kind=dp), DIMENSION(:), POINTER               :: rvals

      CALL sirius_option_get_length(section_name, number_of_options)

      DO elem = 0, number_of_options - 1
         option_name = CHAR(0)
         CALL sirius_option_get_name_and_type(section_name, elem, option_name, ctype)
         option_name1 = TRIM(ADJUSTL(option_name))
         option_name = TRIM(ADJUSTL(option_name))//CHAR(0)

         IF ((option_name1 /= 'memory_usage') .AND. (option_name1 /= 'xc_functionals')) THEN
            CALL section_vals_val_get(section, option_name1, explicit=found)
            IF (found) THEN
               SELECT CASE (ctype)
               CASE (1)
                  CALL section_vals_val_get(section, option_name1, i_val=ival)
                  CALL sirius_option_set_int(sctx, section_name, option_name, ival, 0)
               CASE (11)
                  CALL section_vals_val_get(section, option_name1, i_vals=ivals)
                  CALL sirius_option_set_int(sctx, section_name, option_name, ivals(1), SIZE(ivals))
               CASE (2)
                  CALL section_vals_val_get(section, option_name1, r_val=rval)
                  CALL sirius_option_set_double(sctx, section_name, option_name, rval, 0)
               CASE (12)
                  CALL section_vals_val_get(section, option_name1, r_vals=rvals)
                  CALL sirius_option_set_double(sctx, section_name, option_name, rvals(1), SIZE(rvals))
               CASE (3)
                  CALL section_vals_val_get(section, option_name1, l_val=lval)
                  IF (lval) THEN
                     CALL sirius_option_set_logical(sctx, section_name, option_name, 1, 0)
                  ELSE
                     CALL sirius_option_set_logical(sctx, section_name, option_name, 0, 0)
                  ENDIF
               CASE (13)
                  CALL section_vals_val_get(section, option_name, l_vals=lvals)
                  length = SIZE(lvals)
                  DO i = 1, length
                     IF (lvals(i)) THEN
                        lvalsi(i) = 1
                     ELSE
                        lvalsi(i) = 0
                     ENDIF
                  END DO
                  CALL sirius_option_set_logical(sctx, section_name, option_name, lvalsi(1), length)
               CASE (4)         ! string nightmare
                  CALL section_vals_val_get(section, option_name1, c_val=str)
                  str = TRIM(ADJUSTL(str))//CHAR(0)
                  CALL sirius_option_set_string(sctx, section_name, option_name, str)
               CASE (14)
                  CALL section_vals_val_get(section, option_name1, n_rep_val=length)
                  DO j = 1, length
                     CALL section_vals_val_get(section, option_name1, i_rep_val=j, c_vals=tmp)
                     CALL sirius_option_add_string_to(sctx, section_name, option_name, str)
                  END DO
               CASE DEFAULT
               END SELECT
            END IF
         END IF
      END DO
   END SUBROUTINE cp_sirius_fill_in_section

!***************************************************************************************************
!> \brief ...
!> \param pwdft_env ...
!> \param calculate_forces ...
!> \param calculate_stress_tensor ...
!> \param
!> \par History
!>      07.2018 start the Sirius library
!> \author JHU
! **************************************************************************************************
   SUBROUTINE cp_sirius_energy_force(pwdft_env, calculate_forces, calculate_stress_tensor)
      TYPE(pwdft_environment_type), INTENT(INOUT), &
         POINTER                                         :: pwdft_env
      LOGICAL, INTENT(IN)                                :: calculate_forces, calculate_stress_tensor

      INTEGER                                            :: n1, n2
      LOGICAL                                            :: do_print
      REAL(KIND=C_DOUBLE)                                :: etotal
      REAL(KIND=C_DOUBLE), ALLOCATABLE, DIMENSION(:, :)  :: cforces
      REAL(KIND=C_DOUBLE), DIMENSION(3, 3)               :: cstress
      REAL(KIND=dp), DIMENSION(3, 3)                     :: stress
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: forces
      TYPE(C_PTR)                                        :: gs_handler
      TYPE(pwdft_energy_type), POINTER                   :: energy
      TYPE(section_vals_type), POINTER                   :: print_section, pwdft_input

      CPASSERT(ASSOCIATED(pwdft_env))

      gs_handler = C_NULL_PTR

      CALL pwdft_env_get(pwdft_env=pwdft_env, gs_handler=gs_handler)
      CALL sirius_find_ground_state(gs_handler)
      CALL pwdft_env_get(pwdft_env=pwdft_env, energy=energy)
      etotal = 0.0_C_DOUBLE

      CALL sirius_get_energy(gs_handler, 'total', etotal)
      energy%etotal = etotal

      IF (calculate_forces) THEN
         CALL pwdft_env_get(pwdft_env=pwdft_env, forces=forces)
         n1 = SIZE(forces, 1)
         n2 = SIZE(forces, 2)

         ALLOCATE (cforces(n2, n1))
         cforces = 0.0_C_DOUBLE
         CALL sirius_get_forces(gs_handler, 'total', cforces(1, 1))
         ! Sirius computes the forces but cp2k use the gradient everywhere
         ! so a minus sign is needed.
         ! note also that sirius and cp2k store the forces transpose to each other
         ! sirius : forces(coordinates, atoms)
         ! cp2k : forces(atoms, coordinates)
         forces = -TRANSPOSE(cforces(:, :))
         DEALLOCATE (cforces)
      ENDIF

      IF (calculate_stress_tensor) THEN
         cstress = 0.0_C_DOUBLE
         CALL sirius_get_stress_tensor(gs_handler, 'total', cstress(1, 1))
         stress(1:3, 1:3) = cstress(1:3, 1:3)
         CALL pwdft_env_set(pwdft_env=pwdft_env, stress=stress)
      ENDIF

      CALL pwdft_env_get(pwdft_env=pwdft_env, pwdft_input=pwdft_input)
      print_section => section_vals_get_subs_vals(pwdft_input, "PRINT")
      CALL section_vals_get(print_section, explicit=do_print)
      IF (do_print) THEN
         CALL cp_sirius_print_results(pwdft_env, print_section)
      END IF

   END SUBROUTINE cp_sirius_energy_force

!***************************************************************************************************
!> \brief ...
!> \param pwdft_env ...
!> \param print_section ...
!> \param
!> \par History
!>      12.2019 init
!> \author JHU
! **************************************************************************************************
   SUBROUTINE cp_sirius_print_results(pwdft_env, print_section)
      TYPE(pwdft_environment_type), INTENT(INOUT), &
         POINTER                                         :: pwdft_env
      TYPE(section_vals_type), POINTER                   :: print_section

      CHARACTER(LEN=default_string_length)               :: my_act, my_pos
      INTEGER                                            :: i, ik, iounit, ispn, iterstep, iv, iw, &
                                                            nbands, nhist, nkpts, nspins
      INTEGER(KIND=C_INT)                                :: cint
      LOGICAL                                            :: append, dos, ionode
      REAL(KIND=C_DOUBLE)                                :: creal
      REAL(KIND=C_DOUBLE), ALLOCATABLE, DIMENSION(:)     :: slist
      REAL(KIND=dp)                                      :: de, e_fermi(2), emax, emin, eval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: wkpt
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ehist, hist, occval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: energies, occupations
      TYPE(C_PTR)                                        :: gs_handler, ks_handler, sctx
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()
      ionode = logger%para_env%ionode
      iounit = cp_logger_get_default_io_unit(logger)

      ! Density of States
      dos = BTEST(cp_print_key_should_output(logger%iter_info, print_section, "DOS"), cp_p_file)
      IF (dos) THEN
         ks_handler = C_NULL_PTR
         CALL pwdft_env_get(pwdft_env, ks_handler=ks_handler)
         gs_handler = C_NULL_PTR
         CALL pwdft_env_get(pwdft_env, gs_handler=gs_handler)
         sctx = C_NULL_PTR
         CALL pwdft_env_get(pwdft_env, sctx=sctx)

         CALL section_vals_val_get(print_section, "DOS%DELTA_E", r_val=de)
         CALL section_vals_val_get(print_section, "DOS%APPEND", l_val=append)

         CALL sirius_get_num_kpoints(ks_handler, cint)
         nkpts = cint
         CALL sirius_get_num_bands(sctx, cint)
         nbands = cint
         CALL sirius_get_num_spin_components(sctx, cint)
         nspins = cint
         e_fermi(:) = 0.0_dp
         ALLOCATE (energies(nbands, nspins, nkpts))
         energies = 0.0_dp
         ALLOCATE (occupations(nbands, nspins, nkpts))
         occupations = 0.0_dp
         ALLOCATE (wkpt(nkpts))
         ALLOCATE (slist(nbands))
         DO ik = 0, nkpts - 1
            CALL sirius_get_kpoint_properties(ks_handler, ik, creal)
            wkpt(ik + 1) = creal
         END DO
         DO ik = 1, nkpts
            DO ispn = 0, nspins - 1
               CALL sirius_get_band_energies(ks_handler, ik, ispn, slist(1))
               energies(1:nbands, ispn + 1, ik) = slist(1:nbands)
               CALL sirius_get_band_occupancies(ks_handler, ik, ispn, slist(1))
               occupations(1:nbands, ispn + 1, ik) = slist(1:nbands)
            END DO
         END DO
         emin = MINVAL(energies)
         emax = MAXVAL(energies)
         nhist = NINT((emax - emin)/de) + 1
         ALLOCATE (hist(nhist, nspins), occval(nhist, nspins), ehist(nhist, nspins))
         hist = 0.0_dp
         occval = 0.0_dp
         ehist = 0.0_dp

         DO ik = 1, nkpts
            DO ispn = 1, nspins
               DO i = 1, nbands
                  eval = energies(i, ispn, ik) - emin
                  iv = NINT(eval/de) + 1
                  CPASSERT((iv > 0) .AND. (iv <= nhist))
                  hist(iv, ispn) = hist(iv, ispn) + wkpt(ik)
                  occval(iv, ispn) = occval(iv, ispn) + wkpt(ik)*occupations(i, ispn, ik)
               END DO
            END DO
         END DO
         hist = hist/REAL(nbands, KIND=dp)
         DO i = 1, nhist
            ehist(i, 1:nspins) = emin + (i - 1)*de
         END DO

         iterstep = logger%iter_info%iteration(logger%iter_info%n_rlevel)
         my_act = "WRITE"
         IF (append .AND. iterstep > 1) THEN
            my_pos = "APPEND"
         ELSE
            my_pos = "REWIND"
         END IF

         iw = cp_print_key_unit_nr(logger, print_section, "DOS", &
                                   extension=".dos", file_position=my_pos, file_action=my_act, &
                                   file_form="FORMATTED")
         IF (iw > 0) THEN
            IF (nspins == 2) THEN
               WRITE (UNIT=iw, FMT="(T2,A,I0,A,2F12.6)") &
                  "# DOS at iteration step i = ", iterstep, ", E_Fermi[a.u.] = ", e_fermi(1:2)
               WRITE (UNIT=iw, FMT="(T2,A, A)") "   Energy[a.u.]  Alpha_Density     Occupation", &
                  "   Beta_Density      Occupation"
            ELSE
               WRITE (UNIT=iw, FMT="(T2,A,I0,A,F12.6)") &
                  "# DOS at iteration step i = ", iterstep, ", E_Fermi[a.u.] = ", e_fermi(1)
               WRITE (UNIT=iw, FMT="(T2,A)") "   Energy[a.u.]       Density     Occupation"
            END IF
            DO i = 1, nhist
               eval = emin + (i - 1)*de
               IF (nspins == 2) THEN
                  WRITE (UNIT=iw, FMT="(F15.8,4F15.4)") eval, hist(i, 1), occval(i, 1), &
                     hist(i, 2), occval(i, 2)
               ELSE
                  WRITE (UNIT=iw, FMT="(F15.8,2F15.4)") eval, hist(i, 1), occval(i, 1)
               END IF
            END DO
         END IF
         CALL cp_print_key_finished_output(iw, logger, print_section, "DOS")

         DEALLOCATE (energies, occupations, wkpt, slist)
         DEALLOCATE (hist, occval, ehist)

      END IF

   END SUBROUTINE cp_sirius_print_results

END MODULE sirius_interface

#else

!***************************************************************************************************
!> \brief Empty implementation in case SIRIUS is not compiled in.
!***************************************************************************************************
MODULE sirius_interface
   USE pwdft_environment_types, ONLY: pwdft_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: cp_sirius_init, cp_sirius_finalize
   PUBLIC :: cp_sirius_create_env, cp_sirius_energy_force, cp_sirius_update_context

CONTAINS

! **************************************************************************************************
!> \brief Empty implementation in case SIRIUS is not compiled in.
! **************************************************************************************************
   SUBROUTINE cp_sirius_init()
   END SUBROUTINE cp_sirius_init

! **************************************************************************************************
!> \brief Empty implementation in case SIRIUS is not compiled in.
! **************************************************************************************************
   SUBROUTINE cp_sirius_finalize()
   END SUBROUTINE cp_sirius_finalize

! **************************************************************************************************
!> \brief Empty implementation in case SIRIUS is not compiled in.
!> \param pwdft_env ...
! **************************************************************************************************
   SUBROUTINE cp_sirius_create_env(pwdft_env)
      TYPE(pwdft_environment_type), POINTER              :: pwdft_env

      MARK_USED(pwdft_env)
      CPABORT("Sirius library is missing")
   END SUBROUTINE cp_sirius_create_env

! **************************************************************************************************
!> \brief Empty implementation in case SIRIUS is not compiled in.
!> \param pwdft_env ...
!> \param calculate_forces ...
!> \param calculate_stress ...
! **************************************************************************************************
   SUBROUTINE cp_sirius_energy_force(pwdft_env, calculate_forces, calculate_stress)
      TYPE(pwdft_environment_type), POINTER              :: pwdft_env
      LOGICAL                                            :: calculate_forces, calculate_stress

      MARK_USED(pwdft_env)
      MARK_USED(calculate_forces)
      MARK_USED(calculate_stress)
      CPABORT("Sirius library is missing")
   END SUBROUTINE cp_sirius_energy_force

! **************************************************************************************************
!> \brief Empty implementation in case SIRIUS is not compiled in.
!> \param pwdft_env ...
! **************************************************************************************************
   SUBROUTINE cp_sirius_update_context(pwdft_env)
      TYPE(pwdft_environment_type), POINTER              :: pwdft_env

      MARK_USED(pwdft_env)
      CPABORT("Sirius library is missing")
   END SUBROUTINE cp_sirius_update_context

END MODULE sirius_interface

#endif
