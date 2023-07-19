!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for X-ray absorption spectroscopy using TDDFPT
!> \author AB (01.2018)
!> This subroutine was directly copied from xas_tdp_utils to prevent a
!> circular dependency during the debug-Phase
! **************************************************************************************************

MODULE qs_tddfpt2_soc_utils
   USE cp_control_types,                ONLY: tddfpt2_control_type
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE input_constants,                 ONLY: tddfpt_dipole_berry,&
                                              tddfpt_dipole_length,&
                                              tddfpt_dipole_velocity
   USE kinds,                           ONLY: dp
   USE moments_utils,                   ONLY: get_reference_point
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_operators_ao,                 ONLY: p_xyz_ao,&
                                              rRc_xyz_ao

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_soc_utils'

   PUBLIC :: soc_dipole_operator, soc_contract_evect

   !A helper type for SOC
   TYPE dbcsr_soc_package_type
      TYPE(dbcsr_type), POINTER     :: dbcsr_sg => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_tp => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_sc => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_sf => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_prod => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_ovlp => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_tmp => Null()
      TYPE(dbcsr_type), POINTER     :: dbcsr_work => Null()
   END TYPE dbcsr_soc_package_type

CONTAINS

! **************************************************************************************************
!> \brief Build the atomic dipole operator
!> \param dipmat dipole operator
!> \param tddfpt_control informations on how to build the operaot
!> \param qs_env Qucikstep environment
! **************************************************************************************************
   SUBROUTINE soc_dipole_operator(dipmat, tddfpt_control, qs_env)
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(OUT), &
         POINTER                                         :: dipmat
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'soc_dipole_operator'

      INTEGER                                            :: dim_op, handle, i_dim, nao, nspin
      REAL(kind=dp), DIMENSION(3)                        :: reference_point
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s

      CALL timeset(routineN, handle)

      NULLIFY (matrix_s)

      IF (tddfpt_control%dipole_form == tddfpt_dipole_berry) THEN
         CPABORT("BERRY DIPOLE FORM NOT IMPLEMENTED FOR SOC")
      END IF
           !! ONLY RCS have been implemented, Therefore, nspin sould always be 1!
      nspin = 1
           !! Number of dimensions should be 3, unless multipole is implemented in the future
      dim_op = 3

           !! Initzilize the dipmat structure
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)

      !NULLIFY (dipmat)
      ALLOCATE (dipmat(dim_op))
      DO i_dim = 1, dim_op
         ALLOCATE (dipmat(i_dim)%matrix)
         CALL dbcsr_copy(dipmat(i_dim)%matrix, &
                         matrix_s(1)%matrix, &
                         name="dipole operator matrix")
         !CALL dbcsr_set(dipmat(i_dim)%matrix, 0.0_dp) !!This line needs to go for qs_tddfpt_properties like ao_op
      END DO

      SELECT CASE (tddfpt_control%dipole_form)
      CASE (tddfpt_dipole_length)
                   !!This routine is analog to qs_tddfpt_prperties but only until the rRc_xyz_ao routine
                   !! This will lead to an operator within the nao x nao basis
                   !! qs_tddpft_properies uses nvirt x nocc
         CALL get_reference_point(reference_point, qs_env=qs_env, &
                                  reference=tddfpt_control%dipole_reference, &
                                  ref_point=tddfpt_control%dipole_ref_point)

         CALL rRc_xyz_ao(op=dipmat, qs_env=qs_env, rc=reference_point, order=1, &
                         minimum_image=.FALSE., soft=.FALSE.)
      CASE (tddfpt_dipole_velocity)
                  !!This Routine calcluates the dipole OPerator within the velocity-form within the ao basis
                  !!This Operation is only used in xas_tdp and qs_tddfpt_soc, lines uses rmc_x_p_xyz_ao
         CALL p_xyz_ao(dipmat, qs_env, minimum_image=.FALSE.)
      CASE DEFAULT
         CPABORT("Unimplemented form of the dipole operator")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE soc_dipole_operator

!******************************************************************************
! \brief Will write the tddft Excitation vectors into a single fm
!          It was optaint as a nspin x nstat matrix of fms
! \param fm_start a nspin x nstat Matrix of matries
! \param fm_res one output matrix with the dimensions nao x nocc*nspin*nstat
!******************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param fm_start ...
!> \param fm_res ...
! **************************************************************************************************
   SUBROUTINE soc_contract_evect(fm_start, fm_res)

      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(in)      :: fm_start
      TYPE(cp_fm_type), INTENT(inout)                    :: fm_res

      CHARACTER(len=*), PARAMETER :: routineN = 'soc_contract_evect'

      INTEGER                                            :: handle, ii, jj, nactive, nao, nspins, &
                                                            nstates, ntmp1, ntmp2

      CALL timeset(routineN, handle)

      nstates = SIZE(fm_start, 2)
      nspins = SIZE(fm_start, 1)

      CALL cp_fm_set_all(fm_res, 0.0_dp)
         !! Evects are written into one matrix.
      DO ii = 1, nstates
         DO jj = 1, nspins
            CALL cp_fm_get_info(fm_start(jj, ii), nrow_global=nao, ncol_global=nactive)
            CALL cp_fm_get_info(fm_res, nrow_global=ntmp1, ncol_global=ntmp2)
            CALL cp_fm_to_fm_submat(fm_start(jj, ii), &
                                    fm_res, &
                                    nao, nactive, &
                                    1, 1, 1, &
                                    1 + nactive*(ii - 1) + (jj - 1)*nao*nstates)
         END DO !nspins
      END DO !nsstates

      CALL timestop(handle)

   END SUBROUTINE soc_contract_evect
END MODULE qs_tddfpt2_soc_utils
