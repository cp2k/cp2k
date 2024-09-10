!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of dispersion using pair potentials
!> \author JGH
! **************************************************************************************************
MODULE qs_dispersion_pairpot

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_type
   USE bibliography,                    ONLY: Caldeweyher2017,&
                                              Caldeweyher2019,&
                                              Caldeweyher2020,&
                                              Goerigk2017,&
                                              cite_reference,&
                                              grimme2006,&
                                              grimme2010,&
                                              grimme2011
   USE cell_types,                      ONLY: cell_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                              parser_get_object
   USE cp_parser_types,                 ONLY: cp_parser_type,&
                                              parser_create,&
                                              parser_release
   USE input_constants,                 ONLY: vdw_pairpot_dftd2,&
                                              vdw_pairpot_dftd3,&
                                              vdw_pairpot_dftd3bj,&
                                              vdw_pairpot_dftd4,&
                                              xc_vdw_fun_pairpot
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: bohr,&
                                              kcalmol,&
                                              kjmol
   USE qs_dispersion_cnum,              ONLY: setcn,&
                                              setr0ab
   USE qs_dispersion_d2,                ONLY: calculate_dispersion_d2_pairpot,&
                                              dftd2_param
   USE qs_dispersion_d3,                ONLY: calculate_dispersion_d3_pairpot,&
                                              dftd3_c6_param
   USE qs_dispersion_d4,                ONLY: calculate_dispersion_d4_pairpot
   USE qs_dispersion_types,             ONLY: dftd2_pp,&
                                              dftd3_pp,&
                                              dftd4_pp,&
                                              qs_atom_dispersion_type,&
                                              qs_dispersion_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type,&
                                              set_qs_kind
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_pairpot'

   PUBLIC :: qs_dispersion_pairpot_init, calculate_dispersion_pairpot

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param dispersion_env ...
!> \param pp_section ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE qs_dispersion_pairpot_init(atomic_kind_set, qs_kind_set, dispersion_env, pp_section, para_env)
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      TYPE(section_vals_type), OPTIONAL, POINTER         :: pp_section
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_dispersion_pairpot_init'

      CHARACTER(LEN=2)                                   :: symbol
      CHARACTER(LEN=default_string_length)               :: aname, filename
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: tmpstringlist
      INTEGER                                            :: elem, handle, i, ikind, j, max_elem, &
                                                            maxc, n_rep, nkind, nl, vdw_pp_type, &
                                                            vdw_type
      INTEGER, DIMENSION(:), POINTER                     :: exlist
      LOGICAL                                            :: at_end, explicit, found, is_available
      REAL(KIND=dp)                                      :: dum
      TYPE(qs_atom_dispersion_type), POINTER             :: disp

      CALL timeset(routineN, handle)

      nkind = SIZE(atomic_kind_set)

      vdw_type = dispersion_env%type
      SELECT CASE (vdw_type)
      CASE DEFAULT
         ! do nothing
      CASE (xc_vdw_fun_pairpot)
         ! setup information on pair potentials
         vdw_pp_type = dispersion_env%type
         SELECT CASE (dispersion_env%pp_type)
         CASE DEFAULT
            ! do nothing
         CASE (vdw_pairpot_dftd2)
            CALL cite_reference(Grimme2006)
            DO ikind = 1, nkind
               CALL get_atomic_kind(atomic_kind_set(ikind), element_symbol=symbol, z=elem)
               ALLOCATE (disp)
               disp%type = dftd2_pp
               ! get filename of parameter file
               filename = dispersion_env%parameter_file_name
               ! check for local parameters
               found = .FALSE.
               IF (PRESENT(pp_section)) THEN
                  CALL section_vals_val_get(pp_section, "ATOMPARM", n_rep_val=n_rep)
                  DO i = 1, n_rep
                     CALL section_vals_val_get(pp_section, "ATOMPARM", i_rep_val=i, &
                                               c_vals=tmpstringlist)
                     IF (TRIM(tmpstringlist(1)) == TRIM(symbol)) THEN
                        ! we assume the parameters are in atomic units!
                        READ (tmpstringlist(2), *) disp%c6
                        READ (tmpstringlist(3), *) disp%vdw_radii
                        found = .TRUE.
                        EXIT
                     END IF
                  END DO
               END IF
               IF (.NOT. found) THEN
                  ! check for internal parameters
                  CALL dftd2_param(elem, disp%c6, disp%vdw_radii, found)
               END IF
               IF (.NOT. found) THEN
                  ! check on file
                  INQUIRE (FILE=filename, EXIST=is_available)
                  IF (is_available) THEN
                     BLOCK
                        TYPE(cp_parser_type) :: parser
                        CALL parser_create(parser, filename, para_env=para_env)
                        DO
                           at_end = .FALSE.
                           CALL parser_get_next_line(parser, 1, at_end)
                           IF (at_end) EXIT
                           CALL parser_get_object(parser, aname)
                           IF (TRIM(aname) == TRIM(symbol)) THEN
                              CALL parser_get_object(parser, disp%c6)
                              ! we have to change the units J*nm^6*mol^-1 -> Hartree*Bohr^6
                              disp%c6 = disp%c6*1000._dp*bohr**6/kjmol
                              CALL parser_get_object(parser, disp%vdw_radii)
                              disp%vdw_radii = disp%vdw_radii*bohr
                              found = .TRUE.
                              EXIT
                           END IF
                        END DO
                        CALL parser_release(parser)
                     END BLOCK
                  END IF
               END IF
               IF (found) THEN
                  disp%defined = .TRUE.
               ELSE
                  disp%defined = .FALSE.
               END IF
               ! Check if the parameter is defined
               IF (.NOT. disp%defined) &
                  CALL cp_abort(__LOCATION__, &
                                "Dispersion parameters for element ("//TRIM(symbol)//") are not defined! "// &
                                "Please provide a valid set of parameters through the input section or "// &
                                "through an external file! ")
               CALL set_qs_kind(qs_kind_set(ikind), dispersion=disp)
            END DO
         CASE (vdw_pairpot_dftd3, vdw_pairpot_dftd3bj)
            !DFT-D3 Method initial setup
            CALL cite_reference(Grimme2010)
            CALL cite_reference(Grimme2011)
            CALL cite_reference(Goerigk2017)
            max_elem = 94
            maxc = 5
            dispersion_env%max_elem = max_elem
            dispersion_env%maxc = maxc
            ALLOCATE (dispersion_env%maxci(max_elem))
            ALLOCATE (dispersion_env%c6ab(max_elem, max_elem, maxc, maxc, 3))
            ALLOCATE (dispersion_env%r0ab(max_elem, max_elem))
            ALLOCATE (dispersion_env%rcov(max_elem))
            ALLOCATE (dispersion_env%r2r4(max_elem))
            ALLOCATE (dispersion_env%cn(max_elem))

            ! get filename of parameter file
            filename = dispersion_env%parameter_file_name
            CALL dftd3_c6_param(dispersion_env%c6ab, dispersion_env%maxci, filename, para_env)
            CALL setr0ab(dispersion_env%r0ab, dispersion_env%rcov, dispersion_env%r2r4)
            ! the default coordination numbers
            CALL setcn(dispersion_env%cn)
            ! scale r4/r2 values of the atoms by sqrt(Z)
            ! sqrt is also globally close to optimum
            ! together with the factor 1/2 this yield reasonable
            ! c8 for he, ne and ar. for larger Z, C8 becomes too large
            ! which effectively mimics higher R^n terms neglected due
            ! to stability reasons
            DO i = 1, max_elem
               dum = 0.5_dp*dispersion_env%r2r4(i)*REAL(i, dp)**0.5_dp
               ! store it as sqrt because the geom. av. is taken
               dispersion_env%r2r4(i) = SQRT(dum)
            END DO
            ! parameters
            dispersion_env%k1 = 16.0_dp
            dispersion_env%k2 = 4._dp/3._dp
            ! reasonable choices are between 3 and 5
            ! this gives smoth curves with maxima around the integer values
            ! k3=3 give for CN=0 a slightly smaller value than computed
            ! for the free atom. This also yields to larger CN for atoms
            ! in larger molecules but with the same chem. environment
            ! which is physically not right
            ! values >5 might lead to bumps in the potential
            dispersion_env%k3 = -4._dp
            dispersion_env%rcov = dispersion_env%k2*dispersion_env%rcov*bohr
            ! alpha default parameter
            dispersion_env%alp = 14._dp
            !
            DO ikind = 1, nkind
               CALL get_atomic_kind(atomic_kind_set(ikind), element_symbol=symbol, z=elem)
               ALLOCATE (disp)
               disp%type = dftd3_pp
               IF (elem <= 94) THEN
                  disp%defined = .TRUE.
               ELSE
                  disp%defined = .FALSE.
               END IF
               IF (.NOT. disp%defined) &
                  CALL cp_abort(__LOCATION__, &
                                "Dispersion parameters for element ("//TRIM(symbol)//") are not defined! "// &
                                "Please provide a valid set of parameters through the input section or "// &
                                "through an external file! ")
               CALL set_qs_kind(qs_kind_set(ikind), dispersion=disp)
            END DO

            IF (PRESENT(pp_section)) THEN
               ! Check for coordination numbers
               CALL section_vals_val_get(pp_section, "KIND_COORDINATION_NUMBERS", n_rep_val=n_rep)
               IF (n_rep > 0) THEN
                  ALLOCATE (dispersion_env%cnkind(n_rep))
                  DO i = 1, n_rep
                     CALL section_vals_val_get(pp_section, "KIND_COORDINATION_NUMBERS", i_rep_val=i, &
                                               c_vals=tmpstringlist)
                     READ (tmpstringlist(1), *) dispersion_env%cnkind(i)%cnum
                     READ (tmpstringlist(2), *) dispersion_env%cnkind(i)%kind
                  END DO
               END IF
               CALL section_vals_val_get(pp_section, "ATOM_COORDINATION_NUMBERS", n_rep_val=n_rep)
               IF (n_rep > 0) THEN
                  ALLOCATE (dispersion_env%cnlist(n_rep))
                  DO i = 1, n_rep
                     CALL section_vals_val_get(pp_section, "ATOM_COORDINATION_NUMBERS", i_rep_val=i, &
                                               c_vals=tmpstringlist)
                     nl = SIZE(tmpstringlist)
                     ALLOCATE (dispersion_env%cnlist(i)%atom(nl - 1))
                     dispersion_env%cnlist(i)%natom = nl - 1
                     READ (tmpstringlist(1), *) dispersion_env%cnlist(i)%cnum
                     DO j = 1, nl - 1
                        READ (tmpstringlist(j + 1), *) dispersion_env%cnlist(i)%atom(j)
                     END DO
                  END DO
               END IF
               ! Check for exclusion lists
               CALL section_vals_val_get(pp_section, "D3_EXCLUDE_KIND", explicit=explicit)
               IF (explicit) THEN
                  CALL section_vals_val_get(pp_section, "D3_EXCLUDE_KIND", i_vals=exlist)
                  DO j = 1, SIZE(exlist)
                     ikind = exlist(j)
                     CALL get_qs_kind(qs_kind_set(ikind), dispersion=disp)
                     disp%defined = .FALSE.
                  END DO
               END IF
               CALL section_vals_val_get(pp_section, "D3_EXCLUDE_KIND_PAIR", n_rep_val=n_rep)
               dispersion_env%nd3_exclude_pair = n_rep
               IF (n_rep > 0) THEN
                  ALLOCATE (dispersion_env%d3_exclude_pair(n_rep, 2))
                  DO i = 1, n_rep
                     CALL section_vals_val_get(pp_section, "D3_EXCLUDE_KIND_PAIR", i_rep_val=i, &
                                               i_vals=exlist)
                     dispersion_env%d3_exclude_pair(i, :) = exlist
                  END DO
               END IF
            END IF
         CASE (vdw_pairpot_dftd4)
            !most checks are done by the library
            CALL cite_reference(Caldeweyher2017)
            CALL cite_reference(Caldeweyher2019)
            CALL cite_reference(Caldeweyher2020)
            DO ikind = 1, nkind
               CALL get_atomic_kind(atomic_kind_set(ikind), element_symbol=symbol, z=elem)
               ALLOCATE (disp)
               disp%type = dftd4_pp
               disp%defined = .TRUE.
               CALL set_qs_kind(qs_kind_set(ikind), dispersion=disp)
            END DO
         END SELECT
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE qs_dispersion_pairpot_init

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dispersion_env ...
!> \param energy ...
!> \param calculate_forces ...
!> \param atevdw ...
! **************************************************************************************************
   SUBROUTINE calculate_dispersion_pairpot(qs_env, dispersion_env, energy, calculate_forces, atevdw)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      REAL(KIND=dp), INTENT(INOUT)                       :: energy
      LOGICAL, INTENT(IN)                                :: calculate_forces
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: atevdw

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dispersion_pairpot'

      INTEGER                                            :: atom_a, handle, iatom, ikind, natom, &
                                                            nkind, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      LOGICAL                                            :: atenergy, atex, atstress, debugall, &
                                                            use_virial
      REAL(KIND=dp)                                      :: evdw, gnorm
      REAL(KIND=dp), DIMENSION(3)                        :: fdij
      REAL(KIND=dp), DIMENSION(3, 3)                     :: dvirial, pv_loc, pv_virial_thread
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: atstr
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(virial_type), POINTER                         :: virial

      energy = 0._dp
      ! make valgrind happy
      use_virial = .FALSE.

      IF (dispersion_env%type .NE. xc_vdw_fun_pairpot) THEN
         RETURN
      END IF

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set)

      CALL get_qs_env(qs_env=qs_env, nkind=nkind, natom=natom, atomic_kind_set=atomic_kind_set, &
                      cell=cell, virial=virial, para_env=para_env, atprop=atprop)

      debugall = dispersion_env%verbose

      NULLIFY (logger)
      logger => cp_get_default_logger()
      IF (ASSOCIATED(dispersion_env%dftd_section)) THEN
         unit_nr = cp_print_key_unit_nr(logger, dispersion_env%dftd_section, "PRINT_DFTD", &
                                        extension=".dftd")
      ELSE
         unit_nr = -1
      END IF

      ! atomic energy and stress arrays
      atenergy = atprop%energy
      atstress = atprop%stress
      atstr => atprop%atstress
      ! external atomic energy
      atex = .FALSE.
      IF (PRESENT(atevdw)) THEN
         atex = .TRUE.
      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         WRITE (unit_nr, *) " Pair potential vdW calculation"
         IF (dispersion_env%pp_type == vdw_pairpot_dftd2) THEN
            WRITE (unit_nr, *) " Dispersion potential type: DFT-D2"
            WRITE (unit_nr, *) " Scaling parameter (s6) ", dispersion_env%scaling
            WRITE (unit_nr, *) " Exponential prefactor  ", dispersion_env%exp_pre
         ELSE IF (dispersion_env%pp_type == vdw_pairpot_dftd3) THEN
            WRITE (unit_nr, *) " Dispersion potential type: DFT-D3"
         ELSE IF (dispersion_env%pp_type == vdw_pairpot_dftd3bj) THEN
            WRITE (unit_nr, *) " Dispersion potential type: DFT-D3(BJ)"
         ELSE IF (dispersion_env%pp_type == vdw_pairpot_dftd4) THEN
            WRITE (unit_nr, *) " Dispersion potential type: DFT-D4"
         END IF
      END IF

      CALL get_qs_env(qs_env=qs_env, force=force)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      IF (use_virial .AND. debugall) THEN
         dvirial = virial%pv_virial
      END IF
      IF (use_virial) THEN
         pv_loc = virial%pv_virial
      END IF

      evdw = 0._dp
      pv_virial_thread(:, :) = 0._dp

      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)

      IF (dispersion_env%pp_type == vdw_pairpot_dftd2) THEN
         CALL calculate_dispersion_d2_pairpot(qs_env, dispersion_env, evdw, calculate_forces, atevdw)
      ELSEIF (dispersion_env%pp_type == vdw_pairpot_dftd3 .OR. &
              dispersion_env%pp_type == vdw_pairpot_dftd3bj) THEN
         CALL calculate_dispersion_d3_pairpot(qs_env, dispersion_env, evdw, calculate_forces, &
                                              unit_nr, atevdw)
      ELSEIF (dispersion_env%pp_type == vdw_pairpot_dftd4) THEN
         IF (atenergy) THEN
            CPABORT("Atomic energy with DFTD4 not implemented")
         END IF
         IF (atex) THEN
            CPABORT("Atomic energy with DFTD4 not implemented")
         END IF
         IF (atstress) THEN
            CPABORT("Atomic stress with DFTD4 not implemented")
         END IF
         CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, atomic_kind_set=atomic_kind_set, &
                         cell=cell, force=force, para_env=para_env)
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
         IF (calculate_forces .OR. debugall) THEN
            IF (use_virial) THEN
               CALL calculate_dispersion_d4_pairpot(dispersion_env, particle_set, cell, para_env, &
                                                    evdw, force=force, atom_of_kind=atom_of_kind, &
                                                    virial=pv_virial_thread)
               virial%pv_virial = virial%pv_virial + pv_virial_thread
            ELSE
               CALL calculate_dispersion_d4_pairpot(dispersion_env, particle_set, cell, para_env, &
                                                    evdw, force=force, atom_of_kind=atom_of_kind)
            END IF
         ELSE
            CALL calculate_dispersion_d4_pairpot(dispersion_env, particle_set, cell, para_env, evdw)
         END IF
      END IF

      ! set dispersion energy
      CALL para_env%sum(evdw)
      energy = evdw
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) " Total vdW energy [au]     :", evdw
         WRITE (unit_nr, *) " Total vdW energy [kcal]   :", evdw*kcalmol
         WRITE (unit_nr, *)
      END IF
      IF (calculate_forces .AND. debugall) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) " Dispersion Forces         "
            WRITE (unit_nr, *) " Atom   Kind                            Forces    "
         END IF
         gnorm = 0._dp
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atom_a = atom_of_kind(iatom)
            fdij(1:3) = force(ikind)%dispersion(:, atom_a)
            CALL para_env%sum(fdij)
            gnorm = gnorm + SUM(ABS(fdij))
            IF (unit_nr > 0) WRITE (unit_nr, "(i5,i7,3F20.14)") iatom, ikind, fdij
         END DO
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *)
            WRITE (unit_nr, *) "|G| = ", gnorm
            WRITE (unit_nr, *)
         END IF
         IF (use_virial) THEN
            dvirial = virial%pv_virial - dvirial
            CALL para_env%sum(dvirial)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *) "Stress Tensor (dispersion)"
               WRITE (unit_nr, "(3G20.12)") dvirial
               WRITE (unit_nr, *) "  Tr(P)/3 :  ", (dvirial(1, 1) + dvirial(2, 2) + dvirial(3, 3))/3._dp
               WRITE (unit_nr, *)
            END IF
         END IF
      END IF

      IF (calculate_forces .AND. use_virial) THEN
         virial%pv_vdw = virial%pv_vdw + (virial%pv_virial - pv_loc)
      END IF

      IF (ASSOCIATED(dispersion_env%dftd_section)) THEN
         CALL cp_print_key_finished_output(unit_nr, logger, dispersion_env%dftd_section, "PRINT_DFTD")
      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_dispersion_pairpot

END MODULE qs_dispersion_pairpot
