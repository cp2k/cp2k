!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Auxiliary routines needed for RPA-AXK
!>        given blacs_env to another
!> \par History
!>      09.2016 created [Vladimir Rybkin]
!> \author Vladimir Rybkin
! **************************************************************************************************
MODULE rpa_axk_util
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE hfx_types,                       ONLY: hfx_create,&
                                              hfx_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_axk_util'

   PUBLIC :: hfx_create_subgroup

CONTAINS

! **************************************************************************************************
!> \brief ... Initializes x_data on a subgroup
!> \param qs_env ...
!> \param para_env_sub ...
!> \param hfx_section ...
!> \param x_data ...
!> \param n_rep_hf ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE hfx_create_subgroup(qs_env, para_env_sub, hfx_section, x_data, n_rep_hf)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      TYPE(section_vals_type), POINTER                   :: hfx_section
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      INTEGER                                            :: n_rep_hf

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_create_subgroup', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, natom
      LOGICAL                                            :: do_hfx
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: my_cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      NULLIFY (my_cell, atomic_kind_set, particle_set, dft_control, x_data, qs_kind_set)

      CALL get_qs_env(qs_env, &
                      subsys=subsys, &
                      input=input)

      CALL qs_subsys_get(subsys, &
                         cell=my_cell, &
                         atomic_kind_set=atomic_kind_set, &
                         qs_kind_set=qs_kind_set, &
                         particle_set=particle_set)

      do_hfx = .TRUE.
      hfx_section => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")
      !hfx_section => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(hfx_section, explicit=do_hfx, n_repetition=n_rep_hf)
      CALL get_qs_env(qs_env, dft_control=dft_control)

      IF (do_hfx) THEN
         ! Retrieve particle_set and atomic_kind_set
         natom = SIZE(particle_set)
         CALL hfx_create(x_data, para_env_sub, hfx_section, natom, atomic_kind_set, &
                         qs_kind_set, dft_control, my_cell, do_exx=.TRUE.)
      END IF

      CALL timestop(handle)

   END SUBROUTINE hfx_create_subgroup

END MODULE rpa_axk_util
