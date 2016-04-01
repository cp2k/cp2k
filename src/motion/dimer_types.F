!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Contains types used for a Dimer Method calculations
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! **************************************************************************************************
MODULE dimer_types

   USE cell_types,                      ONLY: use_perd_x,&
                                              use_perd_xy,&
                                              use_perd_xyz,&
                                              use_perd_xz,&
                                              use_perd_y,&
                                              use_perd_yz,&
                                              use_perd_z
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE global_types,                    ONLY: global_environment_type
   USE input_constants,                 ONLY: do_first_rotation_step
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
   USE molecule_kind_types,             ONLY: fixd_constraint_type,&
                                              get_molecule_kind,&
                                              molecule_kind_type
   USE parallel_rng_types,              ONLY: random_numbers
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dimer_types'
   INTEGER, PRIVATE, SAVE :: last_dimer_id = 0

   PUBLIC :: dimer_env_type, &
             dimer_env_create, &
             dimer_env_retain, &
             dimer_env_release, &
             dimer_fixed_atom_control

! **************************************************************************************************
!> \brief Type containing all informations abour the rotation of the Dimer
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! **************************************************************************************************
   TYPE dimer_rotational_type
      ! Rotational parameters
      INTEGER                                    :: rotation_step
      LOGICAL                                    :: interpolate_gradient
      REAL(KIND=dp)                              :: angle_tol, angle1, angle2, dCdp, curvature
      REAL(KIND=dp), POINTER, DIMENSION(:)       :: g0, g1, g1p
   END TYPE dimer_rotational_type

! **************************************************************************************************
!> \brief Type containing all informations abour the translation of the Dimer
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! **************************************************************************************************
   TYPE dimer_translational_type
      ! Translational parameters
      REAL(KIND=dp), POINTER, DIMENSION(:)       :: tls_vec
   END TYPE dimer_translational_type

! **************************************************************************************************
!> \brief Conjugate Directions type
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! **************************************************************************************************
   TYPE dimer_cg_rot_type
      REAL(KIND=dp)                              :: norm_theta, norm_theta_old, norm_h
      REAL(KIND=dp), POINTER, DIMENSION(:)       :: nvec_old
   END TYPE dimer_cg_rot_type

! **************************************************************************************************
!> \brief Defines the environment for a Dimer Method calculation
!> \par History
!>      none
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! **************************************************************************************************
   TYPE dimer_env_type
      INTEGER                                    :: ref_count, id_nr
      REAL(KIND=dp)                              :: dr
      REAL(KIND=dp), POINTER, DIMENSION(:)       :: nvec
      TYPE(dimer_rotational_type)                :: rot
      TYPE(dimer_translational_type)             :: tsl
      TYPE(dimer_cg_rot_type)                    :: cg_rot
   END TYPE dimer_env_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param dimer_env ...
!> \param subsys ...
!> \param globenv ...
!> \param dimer_section ...
!> \par History
!>      2016/03/03 [LTong] changed input natom to subsys
!> \author Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
! **************************************************************************************************
   SUBROUTINE dimer_env_create(dimer_env, subsys, globenv, dimer_section)
      TYPE(dimer_env_type), POINTER                      :: dimer_env
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: dimer_section

      CHARACTER(len=*), PARAMETER :: routineN = 'dimer_env_create', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, isize, j, k, n_rep_val, natom, unit_nr
      LOGICAL                                            :: explicit
      REAL(KIND=dp)                                      :: norm, xval(3)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: array
      TYPE(section_vals_type), POINTER                   :: nvec_section

      unit_nr = cp_logger_get_default_io_unit()
      CPASSERT(.NOT. ASSOCIATED(dimer_env))
      ALLOCATE (dimer_env)
      dimer_env%ref_count = 1
      last_dimer_id = last_dimer_id+1
      dimer_env%id_nr = last_dimer_id
      ! Setup NVEC
      NULLIFY (dimer_env%nvec, dimer_env%rot%g0, dimer_env%rot%g1, dimer_env%rot%g1p, &
               dimer_env%tsl%tls_vec)
      ! get natom
      CALL cp_subsys_get(subsys=subsys, natom=natom)
      ! Allocate the working arrays
      ALLOCATE (dimer_env%nvec(natom*3))
      ALLOCATE (dimer_env%rot%g0(natom*3))
      ALLOCATE (dimer_env%rot%g1(natom*3))
      ALLOCATE (dimer_env%rot%g1p(natom*3))
      ! Check if the dimer vector is available in the input or not..
      nvec_section => section_vals_get_subs_vals(dimer_section, "DIMER_VECTOR")
      CALL section_vals_get(nvec_section, explicit=explicit)
      IF (explicit) THEN
         IF (unit_nr > 0) WRITE (unit_nr, *) "Reading Dimer Vector from file!"
         NULLIFY (array)
         CALL section_vals_val_get(nvec_section, "_DEFAULT_KEYWORD_", n_rep_val=n_rep_val)
         isize = 0
         DO i = 1, n_rep_val
            CALL section_vals_val_get(nvec_section, "_DEFAULT_KEYWORD_", r_vals=array, i_rep_val=i)
            DO j = 1, SIZE(array)
               isize = isize+1
               dimer_env%nvec(isize) = array(j)
            END DO
         END DO
         CPASSERT(isize == SIZE(dimer_env%nvec))
      ELSE
         CALL random_numbers(dimer_env%nvec, globenv%gaussian_rng_stream)
      END IF
      ! Check for translation in the dimer vector and remove them
      IF (natom > 1) THEN
         xval = 0.0_dp
         DO j = 1, natom
            DO k = 1, 3
               i = (j-1)*3+k
               xval(k) = xval(k)+dimer_env%nvec(i)
            END DO
         END DO
         ! Subtract net translations
         xval = xval/REAL(natom*3, KIND=dp)
         DO j = 1, natom
            DO k = 1, 3
               i = (j-1)*3+k
               dimer_env%nvec(i) = dimer_env%nvec(i)-xval(k)
            END DO
         END DO
      END IF
      ! set nvec components to zero for the corresponding constraints
      CALL dimer_fixed_atom_control(dimer_env%nvec, subsys)
      norm = SQRT(SUM(dimer_env%nvec**2))
      IF (norm <= EPSILON(0.0_dp)) &
         CPABORT("The norm of the dimer vector is 0! Calculation cannot proceed further.")
      dimer_env%nvec = dimer_env%nvec/norm
      dimer_env%rot%rotation_step = do_first_rotation_step
      CALL section_vals_val_get(dimer_section, "DR", r_val=dimer_env%dr)
      CALL section_vals_val_get(dimer_section, "INTERPOLATE_GRADIENT", &
                                l_val=dimer_env%rot%interpolate_gradient)
      CALL section_vals_val_get(dimer_section, "ANGLE_TOLERANCE", &
                                r_val=dimer_env%rot%angle_tol)
      ! initialise values
      dimer_env%cg_rot%norm_h = 1.0_dp
      dimer_env%cg_rot%norm_theta = 0.0_dp
      dimer_env%cg_rot%norm_theta_old = 0.0_dp
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

      CHARACTER(len=*), PARAMETER :: routineN = 'dimer_env_retain', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(dimer_env))
      CPASSERT(dimer_env%ref_count > 0)
      dimer_env%ref_count = dimer_env%ref_count+1
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

      CHARACTER(len=*), PARAMETER :: routineN = 'dimer_env_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(dimer_env)) THEN
         CPASSERT(dimer_env%ref_count > 0)
         dimer_env%ref_count = dimer_env%ref_count-1
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
!>        nvec should be set to zero.  Furthermore, the relevent components
!>        of the gradient in CG should also be set to zero.
!> \param vec : vector to be modified
!> \param subsys : subsys type object used by CP2k
!> \par History
!>      2016/03/03 [LTong] created
!> \author Lianheng Tong [LTong]
! **************************************************************************************************
   SUBROUTINE dimer_fixed_atom_control(vec, subsys)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: vec
      TYPE(cp_subsys_type), POINTER                      :: subsys

      CHARACTER(len=*), PARAMETER :: routineN = 'dimer_fixed_atom_control', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ii, ikind, ind, iparticle, nfixed_atoms, &
                                                            nkinds
      TYPE(fixd_constraint_type), DIMENSION(:), POINTER  :: fixd_list
      TYPE(mol_kind_new_list_type), POINTER              :: molecule_kinds
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind

      NULLIFY (molecule_kinds, molecule_kind, fixd_list)

      ! need to get constraint information from molecule information
      CALL cp_subsys_get(subsys=subsys, &
                         molecule_kinds_new=molecule_kinds)
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
            DO ii = 1, nfixed_atoms
               IF (.NOT. fixd_list(ii)%restraint%active) THEN
                  iparticle = fixd_list(ii)%fixd
                  ind = (iparticle-1)*3
                  ! apply constraint to nvec
                  SELECT CASE (fixd_list (ii)%itype)
                  CASE (use_perd_x)
                     vec(ind+1) = 0.0_dp
                  CASE (use_perd_y)
                     vec(ind+2) = 0.0_dp
                  CASE (use_perd_z)
                     vec(ind+3) = 0.0_dp
                  CASE (use_perd_xy)
                     vec(ind+1) = 0.0_dp
                     vec(ind+2) = 0.0_dp
                  CASE (use_perd_xz)
                     vec(ind+1) = 0.0_dp
                     vec(ind+3) = 0.0_dp
                  CASE (use_perd_yz)
                     vec(ind+2) = 0.0_dp
                     vec(ind+3) = 0.0_dp
                  CASE (use_perd_xyz)
                     vec(ind+1) = 0.0_dp
                     vec(ind+2) = 0.0_dp
                     vec(ind+3) = 0.0_dp
                  END SELECT
               END IF ! .NOT.fixd_list(ii)%restraint%active
            END DO ! ii
         END IF ! ASSOCIATED(fixd_list)
      END DO ! ikind
   END SUBROUTINE dimer_fixed_atom_control

END MODULE dimer_types
