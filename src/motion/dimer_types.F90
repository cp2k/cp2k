!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Contains types used for a Dimer Method calculations
!> \par History
!>      Luca Bellucci 11.2017 added kdimer and beta
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! **************************************************************************************************
MODULE dimer_types

   USE cell_types,                      ONLY: &
        cell_transform_input_cartesian, cell_type, use_perd_x, use_perd_xy, use_perd_xyz, &
        use_perd_xz, use_perd_y, use_perd_yz, use_perd_z
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit,&
                                              cp_to_string
   USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                              parser_search_string
   USE cp_parser_types,                 ONLY: cp_parser_type,&
                                              parser_create,&
                                              parser_release
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE force_env_types,                 ONLY: force_env_type
   USE global_types,                    ONLY: global_environment_type
   USE input_constants,                 ONLY: dimer_init_molden,&
                                              dimer_init_random,&
                                              do_first_rotation_step
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE molecule_kind_list_types,        ONLY: molecule_kind_list_type
   USE molecule_kind_types,             ONLY: fixd_constraint_type,&
                                              get_molecule_kind,&
                                              molecule_kind_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dimer_types'

   PUBLIC :: dimer_env_type, &
             dimer_env_create, &
             dimer_env_retain, &
             dimer_env_release, &
             dimer_fixed_atom_control, &
             dimer_init_vector, &
             vib_get_index_weight

! **************************************************************************************************
!> \brief Type containing all informations abour the rotation of the Dimer
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! **************************************************************************************************
   TYPE dimer_rotational_type
      ! Rotational parameters
      INTEGER                                    :: rotation_step = 0
      LOGICAL                                    :: interpolate_gradient = .FALSE.
      REAL(KIND=dp)                              :: angle_tol = 0.0_dp, angle1 = 0.0_dp, angle2 = 0.0_dp, &
                                                    dCdp = 0.0_dp, curvature = 0.0_dp
      REAL(KIND=dp), POINTER, DIMENSION(:)       :: g0 => NULL(), g1 => NULL(), g1p => NULL()
   END TYPE dimer_rotational_type

! **************************************************************************************************
!> \brief Type containing all informations abour the translation of the Dimer
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! **************************************************************************************************
   TYPE dimer_translational_type
      ! Translational parameters
      REAL(KIND=dp), POINTER, DIMENSION(:)       :: tls_vec => NULL()
   END TYPE dimer_translational_type

! **************************************************************************************************
!> \brief Conjugate Directions type
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! **************************************************************************************************
   TYPE dimer_cg_rot_type
      REAL(KIND=dp)                              :: norm_theta = 0.0_dp, norm_theta_old = 0.0_dp, norm_h = 0.0_dp
      REAL(KIND=dp), POINTER, DIMENSION(:)       :: nvec_old => NULL()
   END TYPE dimer_cg_rot_type

! **************************************************************************************************
!> \brief Defines the environment for a Dimer Method calculation
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! **************************************************************************************************
   TYPE dimer_env_type
      INTEGER                                    :: ref_count = 0
      REAL(KIND=dp)                              :: dr = 0.0_dp
      REAL(KIND=dp), POINTER, DIMENSION(:)       :: nvec => NULL()
      REAL(KIND=dp)                              :: beta = 0.0_dp
      TYPE(dimer_rotational_type)                :: rot = dimer_rotational_type()
      TYPE(dimer_translational_type)             :: tsl = dimer_translational_type()
      TYPE(dimer_cg_rot_type)                    :: cg_rot = dimer_cg_rot_type()
      LOGICAL                                    :: kdimer = .FALSE.
   END TYPE dimer_env_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param dimer_env ...
!> \param subsys ...
!> \param globenv ...
!> \param dimer_section ...
!> \param force_env ...
!> \par History
!>      Luca Bellucci 11.2017 added K-DIMER and BETA
!>      2016/03/03 [LTong] changed input natom to subsys
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! **************************************************************************************************
   SUBROUTINE dimer_env_create(dimer_env, subsys, globenv, dimer_section, force_env)
      TYPE(dimer_env_type), POINTER                      :: dimer_env
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: dimer_section
      TYPE(force_env_type), POINTER                      :: force_env

      INTEGER                                            :: i, j, k, natom, unit_nr
      REAL(KIND=dp)                                      :: norm, xval(3)
      TYPE(cell_type), POINTER                           :: cell

      NULLIFY (cell)
      unit_nr = cp_logger_get_default_io_unit()
      CPASSERT(.NOT. ASSOCIATED(dimer_env))
      ALLOCATE (dimer_env)
      dimer_env%ref_count = 1
      ! Setup NVEC
      ! get natom
      CALL cp_subsys_get(subsys=subsys, cell=cell, natom=natom)
      ! Allocate the working arrays
      ALLOCATE (dimer_env%rot%g0(natom*3))
      ALLOCATE (dimer_env%rot%g1(natom*3))
      ALLOCATE (dimer_env%rot%g1p(natom*3))
      ! Read dimer vector from input
      CALL dimer_init_vector(dimer_env, dimer_section, natom, &
                             unit_nr, globenv, force_env)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, "(/,T2,A,T9,A,T71,I10)") &
            "DIMER|", "Dimension of dimer vector (natom * 3)", natom*3
         DO j = 1, natom
            WRITE (unit_nr, "(T2,A,T9,3(F12.6,3X))") &
               "DIMER|", dimer_env%nvec(3*j - 2:3*j)
         END DO
      END IF
      ! Transform input cell
      DO i = 1, natom
         xval(:) = dimer_env%nvec(3*i - 2:3*i)
         CALL cell_transform_input_cartesian(cell, xval)
         dimer_env%nvec(3*i - 2:3*i) = xval(:)
      END DO
      ! Check for translation in the dimer vector and remove them
      IF (natom > 1) THEN
         xval = 0.0_dp
         DO j = 1, natom
            DO k = 1, 3
               i = (j - 1)*3 + k
               xval(k) = xval(k) + dimer_env%nvec(i)
            END DO
         END DO
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, "(/,T2,A,T9,A)") &
               "DIMER|", "Overall translation to be removed from the initial dimer vector"
            WRITE (unit_nr, "(T2,A,T9,3(F12.6,3X))") &
               "DIMER|", xval(1:3)
         END IF
         ! Subtract net translations
         xval = xval/REAL(natom*3, KIND=dp)
         DO j = 1, natom
            DO k = 1, 3
               i = (j - 1)*3 + k
               dimer_env%nvec(i) = dimer_env%nvec(i) - xval(k)
            END DO
         END DO
      END IF
      ! set nvec components to zero for the corresponding constraints
      CALL dimer_fixed_atom_control(dimer_env%nvec, subsys, unit_nr)
      ! Normalize dimer vector
      norm = SQRT(SUM(dimer_env%nvec**2))
      IF (norm <= EPSILON(0.0_dp)) THEN
         CPABORT("The norm of the dimer vector is 0! Calculation cannot proceed further.")
      END IF
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, "(T2,A,T9,A,T69,F12.6)") &
            "DIMER|", "Norm of dimer vector to be normalized by rescaling", norm
      END IF
      dimer_env%nvec = dimer_env%nvec/norm
      dimer_env%rot%rotation_step = do_first_rotation_step
      CALL section_vals_val_get(dimer_section, "DR", r_val=dimer_env%dr)
      CALL section_vals_val_get(dimer_section, "INTERPOLATE_GRADIENT", &
                                l_val=dimer_env%rot%interpolate_gradient)
      CALL section_vals_val_get(dimer_section, "ANGLE_TOLERANCE", &
                                r_val=dimer_env%rot%angle_tol)
      CALL section_vals_val_get(dimer_section, "K-DIMER", &
                                l_val=dimer_env%kdimer)
      CALL section_vals_val_get(dimer_section, "BETA", &
                                r_val=dimer_env%beta)
      ! initialise values
      dimer_env%cg_rot%norm_h = 1.0_dp
      dimer_env%rot%g0 = 0.0_dp
      dimer_env%rot%g1 = 0.0_dp
      dimer_env%rot%g1p = 0.0_dp
      ALLOCATE (dimer_env%cg_rot%nvec_old(natom*3))
   END SUBROUTINE dimer_env_create

! **************************************************************************************************
!> \brief ...
!> \param dimer_env ...
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! **************************************************************************************************
   SUBROUTINE dimer_env_retain(dimer_env)
      TYPE(dimer_env_type), POINTER                      :: dimer_env

      CPASSERT(ASSOCIATED(dimer_env))
      CPASSERT(dimer_env%ref_count > 0)
      dimer_env%ref_count = dimer_env%ref_count + 1
   END SUBROUTINE dimer_env_retain

! **************************************************************************************************
!> \brief ...
!> \param dimer_env ...
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! **************************************************************************************************
   SUBROUTINE dimer_env_release(dimer_env)
      TYPE(dimer_env_type), POINTER                      :: dimer_env

      IF (ASSOCIATED(dimer_env)) THEN
         CPASSERT(dimer_env%ref_count > 0)
         dimer_env%ref_count = dimer_env%ref_count - 1
         IF (dimer_env%ref_count == 0) THEN
            IF (ASSOCIATED(dimer_env%nvec)) THEN
               DEALLOCATE (dimer_env%nvec)
            END IF
            IF (ASSOCIATED(dimer_env%rot%g0)) THEN
               DEALLOCATE (dimer_env%rot%g0)
            END IF
            IF (ASSOCIATED(dimer_env%rot%g1)) THEN
               DEALLOCATE (dimer_env%rot%g1)
            END IF
            IF (ASSOCIATED(dimer_env%rot%g1p)) THEN
               DEALLOCATE (dimer_env%rot%g1p)
            END IF
            IF (ASSOCIATED(dimer_env%cg_rot%nvec_old)) THEN
               DEALLOCATE (dimer_env%cg_rot%nvec_old)
            END IF
            ! No need to deallocate tls_vec (just a pointer to aother local array)
            NULLIFY (dimer_env%tsl%tls_vec)
            DEALLOCATE (dimer_env)
         END IF
      END IF
   END SUBROUTINE dimer_env_release

! **************************************************************************************************
!> \brief Set parts of a given array vec to zero according to fixed atom constraints.
!>        When atoms are (partially) fixed then the relevant components of
!>        nvec should be set to zero.  Furthermore, the relevant components
!>        of the gradient in CG should also be set to zero.
!> \param vec : vector to be modified
!> \param subsys : subsys type object used by CP2k
!> \param unit : unit to write message to
!> \par History
!>      2016/03/03 [LTong] created
!> \author Lianheng Tong [LTong]
! **************************************************************************************************
   SUBROUTINE dimer_fixed_atom_control(vec, subsys, unit)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: vec
      TYPE(cp_subsys_type), POINTER                      :: subsys
      INTEGER, INTENT(IN), OPTIONAL                      :: unit

      INTEGER                                            :: ii, ikind, ind, iparticle, nfixed_atoms, &
                                                            nkinds
      TYPE(fixd_constraint_type), DIMENSION(:), POINTER  :: fixd_list
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind

      NULLIFY (molecule_kinds, molecule_kind, fixd_list)

      ! need to get constraint information from molecule information
      CALL cp_subsys_get(subsys=subsys, &
                         molecule_kinds=molecule_kinds)
      molecule_kind_set => molecule_kinds%els

      ! get total number of fixed atoms
      ! nkinds is the kinds of molecules, not atoms
      nkinds = molecule_kinds%n_els
      DO ikind = 1, nkinds
         molecule_kind => molecule_kind_set(ikind)
         CALL get_molecule_kind(molecule_kind, &
                                nfixd=nfixed_atoms, &
                                fixd_list=fixd_list)
         IF (ASSOCIATED(fixd_list)) THEN
            IF (PRESENT(unit)) THEN
               IF (unit > 0) THEN
                  WRITE (unit, "(/,T2,A,T9,A,T71,I10)") &
                     "DIMER|", "Number of fixed atoms to adjust dimer vector", nfixed_atoms
               END IF
            END IF
            DO ii = 1, nfixed_atoms
               IF (.NOT. fixd_list(ii)%restraint%active) THEN
                  iparticle = fixd_list(ii)%fixd
                  ind = (iparticle - 1)*3
                  ! apply constraint to nvec
                  SELECT CASE (fixd_list(ii)%itype)
                  CASE (use_perd_x)
                     vec(ind + 1) = 0.0_dp
                  CASE (use_perd_y)
                     vec(ind + 2) = 0.0_dp
                  CASE (use_perd_z)
                     vec(ind + 3) = 0.0_dp
                  CASE (use_perd_xy)
                     vec(ind + 1) = 0.0_dp
                     vec(ind + 2) = 0.0_dp
                  CASE (use_perd_xz)
                     vec(ind + 1) = 0.0_dp
                     vec(ind + 3) = 0.0_dp
                  CASE (use_perd_yz)
                     vec(ind + 2) = 0.0_dp
                     vec(ind + 3) = 0.0_dp
                  CASE (use_perd_xyz)
                     vec(ind + 1) = 0.0_dp
                     vec(ind + 2) = 0.0_dp
                     vec(ind + 3) = 0.0_dp
                  END SELECT
               END IF ! .NOT.fixd_list(ii)%restraint%active
            END DO ! ii
         END IF ! ASSOCIATED(fixd_list)
      END DO ! ikind
   END SUBROUTINE dimer_fixed_atom_control

! **************************************************************************************************
!> \brief Read dimer vector from input and store into dimer_env%nvec
!> \param dimer_env ...
!> \param dimer_section ...
!> \param natom ...
!> \param unit ...
!> \param globenv ...
!> \param force_env ...
!> \par History
!>      2026/05  created
!> \author HE Zilong
! **************************************************************************************************
   SUBROUTINE dimer_init_vector(dimer_env, dimer_section, natom, unit, globenv, force_env)
      TYPE(dimer_env_type), INTENT(INOUT), POINTER       :: dimer_env
      TYPE(section_vals_type), INTENT(INOUT), POINTER    :: dimer_section
      INTEGER, INTENT(IN)                                :: natom, unit
      TYPE(global_environment_type), INTENT(IN), POINTER :: globenv
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(LEN=17)                                  :: vib_name
      CHARACTER(LEN=default_path_length)                 :: molden_name
      INTEGER                                            :: dimer_init_method, i, ierr, isize, j, &
                                                            n_rep_val, nvec_size, vib_list_size
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: vib_id_list
      LOGICAL                                            :: explicit, found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: array_r, vib_wt_list
      REAL(KIND=dp), DIMENSION(:), POINTER               :: array
      TYPE(cp_parser_type)                               :: parser
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: nvec_section

      para_env => force_env%para_env
      nvec_section => section_vals_get_subs_vals(dimer_section, "DIMER_VECTOR")
      NULLIFY (array)
      nvec_size = natom*3

      ALLOCATE (dimer_env%nvec(nvec_size))
      dimer_env%nvec(:) = 0.0_dp
      CALL section_vals_get(nvec_section, explicit=explicit)
      IF (explicit) THEN
         IF (unit > 0) WRITE (unit, "(/,T2,A,T9,A)") &
            "DIMER|", "Initial dimer vector read from input DIMER_VECTOR section"
         isize = 0
         CALL section_vals_val_get(nvec_section, "_DEFAULT_KEYWORD_", n_rep_val=n_rep_val)
         DO i = 1, n_rep_val
            CALL section_vals_val_get(nvec_section, "_DEFAULT_KEYWORD_", r_vals=array, i_rep_val=i)
            DO j = 1, SIZE(array)
               isize = isize + 1
               IF (isize <= nvec_size) THEN
                  dimer_env%nvec(isize) = array(j)
               ELSE
                  CPABORT("Size of input DIMER_VECTOR more than natom * 3")
               END IF
            END DO
         END DO
         IF (isize /= nvec_size) THEN
            CPABORT("Size of input DIMER_VECTOR inconsistent with natom * 3")
         END IF
      ELSE
         CALL section_vals_val_get(dimer_section, "INITIALIZATION_METHOD", i_val=dimer_init_method)
         SELECT CASE (dimer_init_method)
         CASE (dimer_init_random)
            IF (unit > 0) WRITE (unit, "(/,T2,A,T9,A)") &
               "DIMER|", "Initial dimer vector generated randomly"
            CALL globenv%gaussian_rng_stream%fill(dimer_env%nvec)
         CASE (dimer_init_molden)
            CALL section_vals_val_get(dimer_section, "VIB_MOLDEN_NAME", explicit=explicit)
            IF (.NOT. explicit) THEN
               CPABORT("INITIALIZATION_METHOD MOLDEN requires VIB_MOLDEN_NAME")
            END IF
            CALL section_vals_val_get(dimer_section, "VIB_MOLDEN_NAME", c_val=molden_name)
            IF (unit > 0) THEN
               WRITE (unit, "(/,T2,A,T9,A)") &
                  "DIMER|", "Initial dimer vector by linear combination of normal modes from:"
               WRITE (unit, "(T2,A,T9,A)") &
                  "DIMER|", TRIM(ADJUSTL(molden_name))
            END IF
            ! Build lists of indices and weights for linear combination
            CALL vib_get_index_weight(dimer_section, vib_id_list, vib_wt_list, &
                                      vib_list_size, unit)
            ! Find vibrational modes in molden and do linear combination
            ALLOCATE (array_r(nvec_size))
            CALL parser_create(parser, molden_name, para_env=para_env, apply_preprocessing=.FALSE.)
            Vib_modes: DO i = 1, vib_list_size
               ! Format below is from molden_utils.F
               WRITE (UNIT=vib_name, FMT='(T2,A,1X,I6)') "vibration", vib_id_list(i)
               CALL parser_search_string(parser, TRIM(vib_name), &
                                         ignore_case=.FALSE., found=found, &
                                         begin_line=.TRUE., search_from_begin_of_file=.TRUE.)
               IF (.NOT. found) CALL cp_abort(__LOCATION__, &
                                              "Could not found <"//vib_name//"> from molden file")
               DO j = 1, natom
                  CALL parser_get_next_line(parser, 1)
                  READ (UNIT=parser%input_line, FMT=*, IOSTAT=ierr) &
                     array_r(3*j - 2), array_r(3*j - 1), array_r(3*j)
                  IF (ierr /= 0) THEN
                     CALL cp_abort(__LOCATION__, &
                                   "Error while reading MOLDEN file: cannot parse the line "// &
                                   TRIM(ADJUSTL(cp_to_string(j)))//" of <"//vib_name//"> "// &
                                   "for components of the normal mode")
                  END IF
               END DO
               dimer_env%nvec(:) = dimer_env%nvec(:) + array_r(:)*vib_wt_list(i)
            END DO Vib_modes
            CALL parser_release(parser)
            IF (ALLOCATED(array_r)) DEALLOCATE (array_r)
            IF (ALLOCATED(vib_id_list)) DEALLOCATE (vib_id_list)
            IF (ALLOCATED(vib_wt_list)) DEALLOCATE (vib_wt_list)
         CASE DEFAULT
            CPABORT("Invalid or not yet implemented dimer initialization method")
         END SELECT
      END IF

   END SUBROUTINE dimer_init_vector

! **************************************************************************************************
!> \brief Read VIB_INDEX and VIB_WEIGHT under DIMER section
!> \param dimer_section ...
!> \param vib_id_list ...
!> \param vib_wt_list ...
!> \param vib_list_size ...
!> \param unit ...
!> \par History
!>      2026/05  created
!> \author HE Zilong
! **************************************************************************************************
   SUBROUTINE vib_get_index_weight(dimer_section, vib_id_list, vib_wt_list, vib_list_size, unit)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: dimer_section
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: vib_id_list
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: vib_wt_list
      INTEGER, INTENT(OUT)                               :: vib_list_size
      INTEGER, INTENT(IN)                                :: unit

      INTEGER                                            :: i, isize, n_rep_val, vib_id_size, &
                                                            vib_wt_size
      INTEGER, DIMENSION(:), POINTER                     :: tmpilist
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tmprlist

      vib_id_size = 0
      CALL section_vals_val_get(dimer_section, "VIB_INDEX", n_rep_val=n_rep_val)
      DO i = 1, n_rep_val
         CALL section_vals_val_get(dimer_section, "VIB_INDEX", &
                                   i_vals=tmpilist, i_rep_val=i)
         vib_id_size = vib_id_size + SIZE(tmpilist)
      END DO
      ALLOCATE (vib_id_list(vib_id_size))
      isize = 0
      DO i = 1, n_rep_val
         CALL section_vals_val_get(dimer_section, "VIB_INDEX", &
                                   i_vals=tmpilist, i_rep_val=i)
         vib_id_list(isize + 1:isize + SIZE(tmpilist)) = tmpilist
         isize = isize + SIZE(tmpilist)
      END DO

      vib_wt_size = 0
      CALL section_vals_val_get(dimer_section, "VIB_WEIGHT", n_rep_val=n_rep_val)
      DO i = 1, n_rep_val
         CALL section_vals_val_get(dimer_section, "VIB_WEIGHT", &
                                   r_vals=tmprlist, i_rep_val=i)
         vib_wt_size = vib_wt_size + SIZE(tmprlist)
      END DO
      ALLOCATE (vib_wt_list(vib_wt_size))
      isize = 0
      DO i = 1, n_rep_val
         CALL section_vals_val_get(dimer_section, "VIB_WEIGHT", &
                                   r_vals=tmprlist, i_rep_val=i)
         vib_wt_list(isize + 1:isize + SIZE(tmprlist)) = tmprlist
         isize = isize + SIZE(tmprlist)
      END DO

      IF (vib_id_size /= vib_wt_size) THEN
         CALL cp_abort(__LOCATION__, &
                       "Inconsistent count of values speficied in input between "// &
                       "VIB_INDEX ("//TRIM(ADJUSTL(cp_to_string(vib_id_size)))//") "// &
                       "and VIB_WEIGHT ("//TRIM(ADJUSTL(cp_to_string(vib_wt_size)))//")")
      ELSE
         vib_list_size = vib_id_size
      END IF
      IF (unit > 0) THEN
         WRITE (unit, "(/,T2,A,T9,A)") &
            "DIMER|", "Indices and weights of vibrational modes in linear combination"
         DO i = 1, vib_list_size
            WRITE (unit, "(T2,A,T10,I6,T18,F12.6)") &
               "DIMER|", vib_id_list(i), vib_wt_list(i)
         END DO
      END IF

   END SUBROUTINE vib_get_index_weight

END MODULE dimer_types
