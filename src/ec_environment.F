!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Energy correction environment setup and handling
!> \par History
!>       2019.09 created
!> \author JGH
! **************************************************************************************************
MODULE ec_environment
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_container_types,       ONLY: add_basis_set_to_container,&
                                              remove_basis_from_container
   USE basis_set_types,                 ONLY: copy_gto_basis_set,&
                                              create_primitive_basis_set,&
                                              gto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE ec_env_types,                    ONLY: energy_correction_type
   USE input_constants,                 ONLY: ec_functional_harris,&
                                              xc_vdw_fun_nonloc,&
                                              xc_vdw_fun_pairpot
   USE input_cp2k_check,                ONLY: xc_functionals_expand
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: init_orbital_pointers
   USE qs_dispersion_nonloc,            ONLY: qs_dispersion_nonloc_init
   USE qs_dispersion_pairpot,           ONLY: qs_dispersion_pairpot_init
   USE qs_dispersion_types,             ONLY: qs_dispersion_type
   USE qs_dispersion_utils,             ONLY: qs_dispersion_env_set
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE string_utilities,                ONLY: uppercase
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_environment'

   PUBLIC :: ec_env_create

CONTAINS

! **************************************************************************************************
!> \brief Allocates and intitializes ec_env
!> \param qs_env ...
!> \param ec_env the object to create
!> \param dft_section ...
!> \par History
!>       2019.09 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ec_env_create(qs_env, ec_env, dft_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_env_create', routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(ec_env))
      ALLOCATE (ec_env)
      CALL init_ec_env(qs_env, ec_env, dft_section)

   END SUBROUTINE ec_env_create

! **************************************************************************************************
!> \brief Initializes ec_env
!> \param qs_env ...
!> \param ec_env ...
!> \param dft_section ...
!> \par History
!>       2019.09 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE init_ec_env(qs_env, ec_env, dft_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(section_vals_type), OPTIONAL, POINTER         :: dft_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'init_ec_env', routineP = moduleN//':'//routineN

      INTEGER                                            :: ikind, maxlgto, nkind
      LOGICAL                                            :: explicit
      REAL(KIND=dp)                                      :: eps_pgf_orb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set, harris_basis
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(section_vals_type), POINTER                   :: nl_section, pp_section, section1, &
                                                            section2, xc_section

      NULLIFY (atomic_kind_set, dispersion_env, para_env)
      NULLIFY (ec_env%sab_orb, ec_env%sac_ppl, ec_env%sap_ppnl)
      NULLIFY (ec_env%matrix_ks, ec_env%matrix_h, ec_env%matrix_s)
      NULLIFY (ec_env%matrix_t, ec_env%matrix_p, ec_env%matrix_w)
      NULLIFY (ec_env%task_list)
      NULLIFY (ec_env%force)
      NULLIFY (ec_env%mao_coef)
      NULLIFY (ec_env%dispersion_env)
      NULLIFY (ec_env%xc_section)
      NULLIFY (ec_env%cpmos)
      NULLIFY (ec_env%matrix_hz)
      NULLIFY (ec_env%p_env)
      NULLIFY (ec_env%vh_rspace)
      NULLIFY (ec_env%vxc_rspace)
      NULLIFY (ec_env%vtau_rspace)
      ec_env%should_update = .TRUE.
      ec_env%mao = .FALSE.

      IF (qs_env%energy_correction) THEN
         CALL section_vals_val_get(dft_section, "ENERGY_CORRECTION%ALGORITHM", &
                                   i_val=ec_env%ks_solver)
         CALL section_vals_val_get(dft_section, "ENERGY_CORRECTION%ENERGY_FUNCTIONAL", &
                                   i_val=ec_env%energy_functional)
         CALL section_vals_val_get(dft_section, "ENERGY_CORRECTION%FACTORIZATION", &
                                   i_val=ec_env%factorization)
         CALL section_vals_val_get(dft_section, "ENERGY_CORRECTION%EPS_DEFAULT", &
                                   r_val=ec_env%eps_default)
         CALL section_vals_val_get(dft_section, "ENERGY_CORRECTION%HARRIS_BASIS", &
                                   c_val=ec_env%basis)
         CALL section_vals_val_get(dft_section, "ENERGY_CORRECTION%MAO", &
                                   l_val=ec_env%mao)
         CALL section_vals_val_get(dft_section, "ENERGY_CORRECTION%MAO_MAX_ITER", &
                                   i_val=ec_env%mao_max_iter)
         CALL section_vals_val_get(dft_section, "ENERGY_CORRECTION%MAO_EPS_GRAD", &
                                   r_val=ec_env%mao_eps_grad)

         ! set basis
         CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, nkind=nkind)
         CALL uppercase(ec_env%basis)
         SELECT CASE (ec_env%basis)
         CASE ("ORBITAL")
            DO ikind = 1, nkind
               qs_kind => qs_kind_set(ikind)
               CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="ORB")
               IF (ASSOCIATED(basis_set)) THEN
                  NULLIFY (harris_basis)
                  CALL get_qs_kind(qs_kind=qs_kind, basis_set=harris_basis, basis_type="HARRIS")
                  IF (ASSOCIATED(harris_basis)) THEN
                     CALL remove_basis_from_container(qs_kind%basis_sets, basis_type="HARRIS")
                  END IF
                  NULLIFY (harris_basis)
                  CALL copy_gto_basis_set(basis_set, harris_basis)
                  CALL add_basis_set_to_container(qs_kind%basis_sets, harris_basis, "HARRIS")
               END IF
            END DO
         CASE ("PRIMITIVE")
            DO ikind = 1, nkind
               qs_kind => qs_kind_set(ikind)
               CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="ORB")
               IF (ASSOCIATED(basis_set)) THEN
                  NULLIFY (harris_basis)
                  CALL get_qs_kind(qs_kind=qs_kind, basis_set=harris_basis, basis_type="HARRIS")
                  IF (ASSOCIATED(harris_basis)) THEN
                     CALL remove_basis_from_container(qs_kind%basis_sets, basis_type="HARRIS")
                  END IF
                  NULLIFY (harris_basis)
                  CALL create_primitive_basis_set(basis_set, harris_basis)
                  CALL get_qs_env(qs_env, dft_control=dft_control)
                  eps_pgf_orb = dft_control%qs_control%eps_pgf_orb
                  CALL init_interaction_radii_orb_basis(harris_basis, eps_pgf_orb)
                  harris_basis%kind_radius = basis_set%kind_radius
                  CALL add_basis_set_to_container(qs_kind%basis_sets, harris_basis, "HARRIS")
               END IF
            END DO
         CASE ("HARRIS")
            DO ikind = 1, nkind
               qs_kind => qs_kind_set(ikind)
               NULLIFY (harris_basis)
               CALL get_qs_kind(qs_kind=qs_kind, basis_set=harris_basis, basis_type="HARRIS")
               IF (.NOT. ASSOCIATED(harris_basis)) THEN
                  CPWARN("Harris Basis not defined for all types of atoms.")
               END IF
            END DO
         CASE DEFAULT
            CPABORT("Unknown basis set for energy correction (Harris functional)")
         END SELECT
         !
         CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto, basis_type="HARRIS")
         CALL init_orbital_pointers(maxlgto + 1)
         ! set functional
         SELECT CASE (ec_env%energy_functional)
         CASE (ec_functional_harris)
            ec_env%ec_name = "Harris"
         CASE DEFAULT
            CPABORT("unknown energy correction")
         END SELECT
         ! select the XC section
         NULLIFY (xc_section)
         xc_section => section_vals_get_subs_vals(dft_section, "XC")
         section1 => section_vals_get_subs_vals(dft_section, "ENERGY_CORRECTION%XC")
         section2 => section_vals_get_subs_vals(dft_section, "ENERGY_CORRECTION%XC%XC_FUNCTIONAL")
         CALL section_vals_get(section2, explicit=explicit)
         IF (explicit) THEN
            CALL xc_functionals_expand(section2, section1)
            ec_env%xc_section => section1
         ELSE
            ec_env%xc_section => xc_section
         END IF
         ! dispersion
         ALLOCATE (dispersion_env)
         NULLIFY (xc_section)
         xc_section => ec_env%xc_section
         CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, para_env=para_env)
         CALL qs_dispersion_env_set(dispersion_env, xc_section)
         IF (dispersion_env%type == xc_vdw_fun_pairpot) THEN
            NULLIFY (pp_section)
            pp_section => section_vals_get_subs_vals(xc_section, "VDW_POTENTIAL%PAIR_POTENTIAL")
            CALL qs_dispersion_pairpot_init(atomic_kind_set, qs_kind_set, dispersion_env, pp_section, para_env)
         ELSE IF (dispersion_env%type == xc_vdw_fun_nonloc) THEN
            NULLIFY (nl_section)
            nl_section => section_vals_get_subs_vals(xc_section, "VDW_POTENTIAL%NON_LOCAL")
            CALL qs_dispersion_nonloc_init(dispersion_env, para_env)
         END IF
         ec_env%dispersion_env => dispersion_env
      END IF

   END SUBROUTINE init_ec_env

END MODULE ec_environment
