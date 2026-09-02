!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shared assembly of complex k-point operators: real-space DBCSR matrices to
!>        a complex full matrix on the k-point group distribution.
!>        The work is split so that bulk MPI collectives can be started for all operators
!>        and finished later, as done by the STANDARD k-point driver.
!> \author MI
! **************************************************************************************************
MODULE qs_kp_evp_methods
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale_and_add_fm
   USE cp_cfm_types,                    ONLY: cp_cfm_to_fm,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_desymmetrize,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_types,                     ONLY: copy_info_type,&
                                              cp_fm_cleanup_copy_general,&
                                              cp_fm_finish_copy_general,&
                                              cp_fm_start_copy_general,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: rskp_grid_type,&
                                              rskp_transform,&
                                              rskp_transform_grid_extract
   USE mathconstants,                   ONLY: gaussi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kp_evp_methods'

   PUBLIC :: kp_evp_op_start, kp_evp_op_finish_cfm, kp_evp_cfm_to_mo

CONTAINS

! **************************************************************************************************
!> \brief Transform one real-space operator to a k point and start its redistribution.
!>        The real and imaginary parts are densified into fm_re and fm_im, and their
!>        transfer to the k-point group distribution is started. Both transfers receive
!>        into private buffers, so any number of starts may target the same fmlocal.
!>        The results are collected with kp_evp_op_finish_cfm.
!> \param rsmat real-space image matrices of the operator
!> \param ispin spin component of rsmat (1 for spin-free matrices such as S)
!> \param ik global k-point index
!> \param xkp coordinates of this k point
!> \param cell_to_index mapping of cell indices to real-space index
!> \param sab_nl neighbor lists defining the real-space sparsity
!> \param grid prepared reciprocal-grid cache, enables the FFT branch with use_grid
!> \param use_grid extract from grid instead of a direct phase sum
!> \param rmatrix workspace holding the real part of the k-point matrix
!> \param cmatrix workspace holding the imaginary part of the k-point matrix
!> \param tmpmat unsymmetric DBCSR workspace
!> \param fm_re full-matrix workspace of the global communicator holding the real part
!> \param fm_im full-matrix workspace of the global communicator holding the imaginary part
!> \param target_fm transfer target on the k-point group (fmlocal, or fmdummy to keep
!>        collective calls balanced on groups one does not own)
!> \param para_env communicator of source and target distributions
!> \param info_re transfer state of the real part
!> \param info_im transfer state of the imaginary part
! **************************************************************************************************
   SUBROUTINE kp_evp_op_start(rsmat, ispin, ik, xkp, cell_to_index, sab_nl, grid, use_grid, &
                              rmatrix, cmatrix, tmpmat, fm_re, fm_im, target_fm, para_env, info_re, info_im)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat
      INTEGER, INTENT(IN)                                :: ispin, ik
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(rskp_grid_type), INTENT(IN), OPTIONAL         :: grid
      LOGICAL, INTENT(IN)                                :: use_grid
      TYPE(dbcsr_type)                                   :: rmatrix, cmatrix, tmpmat
      TYPE(cp_fm_type)                                   :: fm_re, fm_im, target_fm
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(copy_info_type)                               :: info_re, info_im

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'kp_evp_op_start'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (use_grid) THEN
         CPASSERT(PRESENT(grid))
         CALL rskp_transform_grid_extract(grid, ik, rmatrix, cmatrix)
      ELSE
         CALL dbcsr_set(rmatrix, 0.0_dp)
         CALL dbcsr_set(cmatrix, 0.0_dp)
         CALL rskp_transform(rmatrix=rmatrix, cmatrix=cmatrix, rsmat=rsmat, ispin=ispin, &
                             xkp=xkp, cell_to_index=cell_to_index, sab_nl=sab_nl)
      END IF
      CALL dbcsr_desymmetrize(rmatrix, tmpmat)
      CALL copy_dbcsr_to_fm(tmpmat, fm_re)
      CALL dbcsr_desymmetrize(cmatrix, tmpmat)
      CALL copy_dbcsr_to_fm(tmpmat, fm_im)

      CALL cp_fm_start_copy_general(fm_re, target_fm, para_env, info_re)
      CALL cp_fm_start_copy_general(fm_im, target_fm, para_env, info_im)

      CALL timestop(handle)

   END SUBROUTINE kp_evp_op_start

! **************************************************************************************************
!> \brief Finish one transfer started by kp_evp_op_start and add it to a complex matrix.
!> \param fmlocal transfer target filled by the finished copy
!> \param info transfer state of one part
!> \param cmat complex matrix accumulating the operator
!> \param imaginary_part add fmlocal as the imaginary part instead of the real part
! **************************************************************************************************
   SUBROUTINE kp_evp_op_finish_cfm(fmlocal, info, cmat, imaginary_part)

      TYPE(cp_fm_type)                                   :: fmlocal
      TYPE(copy_info_type)                               :: info
      TYPE(cp_cfm_type)                                  :: cmat
      LOGICAL, INTENT(IN)                                :: imaginary_part

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kp_evp_op_finish_cfm'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL cp_fm_finish_copy_general(fmlocal, info)
      IF (imaginary_part) THEN
         CALL cp_cfm_scale_and_add_fm(z_one, cmat, gaussi, fmlocal)
      ELSE
         CALL cp_cfm_scale_and_add_fm(z_zero, cmat, z_one, fmlocal)
      END IF

      CALL timestop(handle)

   END SUBROUTINE kp_evp_op_finish_cfm

!> \brief Split a complex matrix into the real and imaginary MO sets and keep the
!>        eigenvalues of both spin-imaginary sets in sync.
!> \param cmat complex MO coefficients
!> \param mo_re real part MO set
!> \param mo_im imaginary part MO set
!> \param eigenvalues eigenvalues of mo_re, copied to mo_im
! **************************************************************************************************
   SUBROUTINE kp_evp_cfm_to_mo(cmat, mo_re, mo_im, eigenvalues)

      TYPE(cp_cfm_type)                                  :: cmat
      TYPE(mo_set_type)                                  :: mo_re, mo_im
      REAL(KIND=dp), DIMENSION(:)                        :: eigenvalues

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'kp_evp_cfm_to_mo'

      INTEGER                                            :: handle
      TYPE(cp_fm_type), POINTER                          :: imos, rmos

      CALL timeset(routineN, handle)

      CALL get_mo_set(mo_re, mo_coeff=rmos)
      CALL get_mo_set(mo_im, mo_coeff=imos)
      mo_im%eigenvalues = eigenvalues
      CALL cp_cfm_to_fm(cmat, rmos, imos)

      CALL timestop(handle)

   END SUBROUTINE kp_evp_cfm_to_mo

END MODULE qs_kp_evp_methods
