!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines used for Harris functional
!>        Kohn-Sham calculation
!> \par History
!>       10.2020 created
!> \author Fabian Belleflamme
! **************************************************************************************************
MODULE ec_methods
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_m_by_n_from_row_template
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_to_string
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_init_p,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              qs_kind_type
   USE qs_matrix_pools,                 ONLY: mpools_create,&
                                              mpools_rebuild_fm_pools,&
                                              mpools_release,&
                                              qs_matrix_pools_type
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                              xc_prep_2nd_deriv
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_release
   USE xc_rho_set_types,                ONLY: xc_rho_set_release,&
                                              xc_rho_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_methods'

   PUBLIC :: create_kernel
   PUBLIC :: ec_mos_init

CONTAINS

! **************************************************************************************************
!> \brief Creation of second derivative xc-potential
!> \param qs_env ...
!> \param vxc will contain the partially integrated second derivative
!>        taken with respect to rho, evaluated in rho and folded with rho1
!>        vxc is allocated here and needs to be deallocated by the caller.
!> \param rho density at which derivatives were calculated
!> \param rho1_r density in r-space with which to fold
!> \param rho1_g density in g-space with which to fold
!> \param xc_section XC parameters of this potential
!> \param compute_virial Enable stress-tensor calculation
!> \param virial_xc Will contain GGA contribution of XC-kernel to stress-tensor
!> \date    11.2019
!> \author  fbelle
! **************************************************************************************************
   SUBROUTINE create_kernel(qs_env, vxc, rho, rho1_r, rho1_g, xc_section, compute_virial, virial_xc)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), DIMENSION(:), INTENT(OUT), &
         POINTER                                         :: vxc
      TYPE(qs_rho_type), INTENT(IN), POINTER             :: rho
      TYPE(pw_p_type), DIMENSION(:), INTENT(IN), POINTER :: rho1_r, rho1_g
      TYPE(section_vals_type), INTENT(IN), POINTER       :: xc_section
      LOGICAL, INTENT(IN), OPTIONAL                      :: compute_virial
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(INOUT), &
         OPTIONAL                                        :: virial_xc

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_kernel'

      INTEGER                                            :: handle
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(xc_rho_set_type), POINTER                     :: rho_set

      CALL timeset(routineN, handle)

      NULLIFY (auxbas_pw_pool, deriv_set, pw_env, rho_r, rho_set)

      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool)
      ! Get density
      CALL qs_rho_get(rho, rho_r=rho_r)

      CALL xc_prep_2nd_deriv(deriv_set=deriv_set, &    ! containing potentials
                             rho_set=rho_set, &        ! density at which derivs are calculated
                             rho_r=rho_r, &            ! place where derivative is evaluated
                             pw_pool=auxbas_pw_pool, & ! pool for grids
                             xc_section=xc_section)

      ! folding of second deriv with density in rho1_set
      CALL xc_calc_2nd_deriv(v_xc=vxc, &               ! XC-potential
                             deriv_set=deriv_set, &    ! deriv of xc-potential
                             rho_set=rho_set, &        ! density at which deriv are calculated
                             rho1_r=rho1_r, &          ! density with which to fold
                             rho1_g=rho1_g, &          ! density with which to fold
                             pw_pool=auxbas_pw_pool, & ! pool for grids
                             xc_section=xc_section, &
                             gapw=.FALSE., &
                             compute_virial=compute_virial, &
                             virial_xc=virial_xc)

      ! Release second deriv stuff
      CALL xc_dset_release(deriv_set)
      CALL xc_rho_set_release(rho_set=rho_set, pw_pool=auxbas_pw_pool)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Allocate and initiate molecular orbitals environment
!>
!> \param qs_env ...
!> \param matrix_s Used as template
!> \param
!>
!> \par History
!>       2020.10 created [Fabian Belleflamme]
!> \author Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_mos_init(qs_env, matrix_s)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type)                                   :: matrix_s

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_mos_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, multiplicity, n_ao, &
                                                            nelectron, nmo, nspins
      INTEGER, DIMENSION(2)                              :: n_mo, nelectron_spin
      REAL(dp)                                           :: maxocc
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_b
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_matrix_pools_type), POINTER                :: my_mpools

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, dft_control, mo_coeff, mo_coeff_b, mos, my_mpools, qs_kind_set)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      blacs_env=blacs_env, &
                      qs_kind_set=qs_kind_set, &
                      nelectron_spin=nelectron_spin, &
                      para_env=para_env)
      nspins = dft_control%nspins

      ! Start setup
      CALL get_qs_kind_set(qs_kind_set, nsgf=n_ao, nelectron=nelectron)

      ! the total number of electrons
      nelectron = nelectron - dft_control%charge
      multiplicity = dft_control%multiplicity

      ! setting maxocc and n_mo
      IF (dft_control%nspins == 1) THEN
         maxocc = 2.0_dp
         nelectron_spin(1) = nelectron
         nelectron_spin(2) = 0
         IF (MODULO(nelectron, 2) == 0) THEN
            n_mo(1) = nelectron/2
         ELSE
            n_mo(1) = INT(nelectron/2._dp) + 1
         END IF
         n_mo(2) = 0
      ELSE
         maxocc = 1.0_dp

         ! The simplist spin distribution is written here. Special cases will
         ! need additional user input
         IF (MODULO(nelectron + multiplicity - 1, 2) /= 0) THEN
            CPABORT("LSD: try to use a different multiplicity")
         END IF

         nelectron_spin(1) = (nelectron + multiplicity - 1)/2
         nelectron_spin(2) = (nelectron - multiplicity + 1)/2

         IF (nelectron_spin(2) < 0) THEN
            CPABORT("LSD: too few electrons for this multiplicity")
         END IF

         n_mo(1) = nelectron_spin(1)
         n_mo(2) = nelectron_spin(2)

      END IF

      ! Allocate MO set
      ALLOCATE (mos(nspins))
      DO ispin = 1, nspins
         NULLIFY (mos(ispin)%mo_set)
         CALL allocate_mo_set(mo_set=mos(ispin)%mo_set, &
                              nao=n_ao, &
                              nmo=n_mo(ispin), &
                              nelectron=nelectron_spin(ispin), &
                              n_el_f=REAL(nelectron_spin(ispin), dp), &
                              maxocc=maxocc, &
                              flexible_electron_count=dft_control%relax_multiplicity)
      END DO

      CALL set_qs_env(qs_env, mos=mos)

      ! init pools
      NULLIFY (my_mpools)
      CALL mpools_create(mpools=my_mpools)
      CALL mpools_rebuild_fm_pools(mpools=my_mpools, &
                                   mos=mos, &
                                   blacs_env=blacs_env, &
                                   para_env=para_env)

      ! finish initialization of the MOs
      CPASSERT(ASSOCIATED(mos))
      NULLIFY (mo_coeff, mo_coeff_b)
      DO ispin = 1, SIZE(mos)
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, mo_coeff_b=mo_coeff_b)
         IF (.NOT. ASSOCIATED(mo_coeff)) THEN
            CALL init_mo_set(mos(ispin)%mo_set, &
                             fm_pool=my_mpools%ao_mo_fm_pools(ispin)%pool, &
                             name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))// &
                             "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))))
         END IF

         IF (.NOT. ASSOCIATED(mo_coeff_b)) THEN
            CALL cp_fm_get_info(mos(ispin)%mo_set%mo_coeff, ncol_global=nmo)
            CALL dbcsr_init_p(mos(ispin)%mo_set%mo_coeff_b)
            CALL cp_dbcsr_m_by_n_from_row_template(mos(ispin)%mo_set%mo_coeff_b, &
                                                   template=matrix_s, &
                                                   n=nmo, &
                                                   sym=dbcsr_type_no_symmetry)
         END IF
      END DO

      CALL mpools_release(mpools=my_mpools)

      CALL timestop(handle)

   END SUBROUTINE ec_mos_init

END MODULE ec_methods
