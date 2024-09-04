!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_assign
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_geadd,&
                                              cp_fm_scale
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE exstates_types,                  ONLY: wfn_history_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_assign'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

   PUBLIC :: assign_state

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_s ...
!> \param evects ...
!> \param psi0 ...
!> \param wfn_history ...
!> \param my_state ...
! **************************************************************************************************
   SUBROUTINE assign_state(qs_env, matrix_s, evects, psi0, wfn_history, my_state)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(cp_fm_type), DIMENSION(:, :)                  :: evects
      TYPE(cp_fm_type), DIMENSION(:)                     :: psi0
      TYPE(wfn_history_type)                             :: wfn_history
      INTEGER, INTENT(INOUT)                             :: my_state

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'assign_state'

      INTEGER                                            :: handle, is, ispin, natom, ncol, nspins, &
                                                            nstate
      REAL(KIND=dp)                                      :: xsum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dv, rdiag
      TYPE(dbcsr_type), POINTER                          :: smat
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, natom=natom, para_env=para_env)
      nspins = SIZE(psi0)
      nstate = SIZE(evects, 2)
      !
      smat => matrix_s(1)%matrix
      !
      IF (ASSOCIATED(wfn_history%evect)) THEN
         ALLOCATE (dv(nstate))
         !
         wfn_history%gsval = 0.0_dp
         wfn_history%gsmin = 1.0_dp
         DO ispin = 1, nspins
            CALL cp_fm_get_info(wfn_history%evect(ispin), ncol_global=ncol)
            CALL lowdin_orthogonalization(wfn_history%cpmos(ispin), wfn_history%evect(ispin), &
                                          ncol, smat)
            ALLOCATE (rdiag(ncol))
            CALL wfn_align(psi0(ispin), wfn_history%cpmos(ispin), wfn_history%evect(ispin), &
                           rdiag, smat)
            wfn_history%gsval = wfn_history%gsval + SUM(rdiag)/REAL(ncol*nspins, KIND=dp)
            wfn_history%gsmin = MIN(wfn_history%gsmin, MINVAL(rdiag))
            DEALLOCATE (rdiag)
         END DO
         DO is = 1, nstate
            xsum = 0.0_dp
            DO ispin = 1, nspins
               CALL cp_fm_get_info(wfn_history%evect(ispin), ncol_global=ncol)
               ALLOCATE (rdiag(ncol))
               CALL xvec_ovlp(evects(ispin, is), wfn_history%evect(ispin), rdiag, smat)
               xsum = xsum + SUM(rdiag)
               DEALLOCATE (rdiag)
            END DO
            dv(is) = ABS(xsum)/SQRT(REAL(nspins, dp))
         END DO
         my_state = MAXVAL(MAXLOC(dv))
         wfn_history%xsval = dv(my_state)
         IF (wfn_history%xsval < 0.75_dp) THEN
            dv(my_state) = 0.0_dp
            IF (wfn_history%xsval/MAXVAL(dv) < 0.5_dp) THEN
               CALL cp_warn(__LOCATION__, "Uncertain assignment for State following."// &
                            " Reduce trust radius in Geometry Optimization or timestep"// &
                            " in MD runs.")
            END IF
         END IF
         DO ispin = 1, nspins
            CALL cp_fm_get_info(wfn_history%evect(ispin), ncol_global=ncol)
            CALL cp_fm_to_fm(evects(ispin, my_state), wfn_history%evect(ispin))
            CALL cp_fm_to_fm(psi0(ispin), wfn_history%cpmos(ispin), ncol, 1, 1)
         END DO
         !
         DEALLOCATE (dv)
      ELSE
         !
         ALLOCATE (wfn_history%evect(nspins))
         ALLOCATE (wfn_history%cpmos(nspins))
         DO ispin = 1, nspins
            CALL cp_fm_create(wfn_history%evect(ispin), evects(ispin, 1)%matrix_struct, "Xvec")
            CALL cp_fm_create(wfn_history%cpmos(ispin), evects(ispin, 1)%matrix_struct, "Cvec")
         END DO
         DO ispin = 1, nspins
            CALL cp_fm_get_info(wfn_history%evect(ispin), ncol_global=ncol)
            CALL cp_fm_to_fm(evects(ispin, my_state), wfn_history%evect(ispin))
            CALL cp_fm_to_fm(psi0(ispin), wfn_history%cpmos(ispin), ncol, 1, 1)
         END DO
         wfn_history%xsval = 1.0_dp
         wfn_history%gsval = 1.0_dp
         wfn_history%gsmin = 1.0_dp
      END IF

      CALL timestop(handle)

   END SUBROUTINE assign_state

! **************************************************************************************************
!> \brief return a set of S orthonormal vectors (C^T S C == 1) where
!>      a Lodwin transformation is applied to keep the rotated vectors as close
!>      as possible to the original ones
!> \param vmatrix ...
!> \param xmatrix ...
!> \param ncol ...
!> \param matrix_s ...
!> \param
!> \par History
!>      05.2009 created [MI]
!>      06.2023 adapted to include a second set of vectors [JGH]
!> \note
! **************************************************************************************************
   SUBROUTINE lowdin_orthogonalization(vmatrix, xmatrix, ncol, matrix_s)

      TYPE(cp_fm_type), INTENT(IN)                       :: vmatrix, xmatrix
      INTEGER, INTENT(IN)                                :: ncol
      TYPE(dbcsr_type)                                   :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'lowdin_orthogonalization'
      REAL(KIND=dp), PARAMETER                           :: rone = 1.0_dp, rzero = 0.0_dp

      INTEGER                                            :: handle, n, ncol_global, ndep
      REAL(dp)                                           :: threshold, xsum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: rdiag
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type)                                   :: csc, sc, work

      IF (ncol .EQ. 0) RETURN

      CALL timeset(routineN, handle)

      threshold = 1.0E-7_dp
      CALL cp_fm_get_info(matrix=vmatrix, nrow_global=n, ncol_global=ncol_global)
      IF (ncol .GT. ncol_global) CPABORT("Wrong ncol value")

      CALL cp_fm_create(sc, xmatrix%matrix_struct, "SC")
      CALL cp_dbcsr_sm_fm_multiply(matrix_s, vmatrix, sc, ncol)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol, &
                               para_env=vmatrix%matrix_struct%para_env, &
                               context=vmatrix%matrix_struct%context)
      CALL cp_fm_create(csc, fm_struct_tmp, "csc")
      CALL cp_fm_create(work, fm_struct_tmp, "work")
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL parallel_gemm('T', 'N', ncol, ncol, n, rone, vmatrix, sc, rzero, csc)
      CALL cp_fm_power(csc, work, -0.5_dp, threshold, ndep)
      CALL parallel_gemm('N', 'N', n, ncol, ncol, rone, vmatrix, csc, rzero, sc)
      CALL cp_fm_to_fm(sc, vmatrix, ncol, 1, 1)
      !
      CALL parallel_gemm('N', 'N', n, ncol, ncol, rone, xmatrix, csc, rzero, sc)
      CALL cp_fm_to_fm(sc, xmatrix, ncol, 1, 1)
      ! projecton for xSv = 0
      CALL cp_dbcsr_sm_fm_multiply(matrix_s, xmatrix, sc, ncol)
      CALL parallel_gemm('T', 'N', ncol, ncol, n, rone, vmatrix, sc, rzero, csc)
      CALL parallel_gemm('N', 'N', n, ncol, ncol, rone, vmatrix, csc, rzero, sc)
      CALL cp_fm_geadd(-1.0_dp, 'N', sc, 1.0_dp, xmatrix)
      ! normalisation
      CALL cp_dbcsr_sm_fm_multiply(matrix_s, xmatrix, sc, ncol)
      CALL parallel_gemm('T', 'N', ncol, ncol, n, rone, xmatrix, sc, rzero, csc)
      ALLOCATE (rdiag(ncol))
      CALL cp_fm_get_diag(csc, rdiag)
      xsum = SUM(rdiag)
      DEALLOCATE (rdiag)
      xsum = 1._dp/SQRT(xsum)
      CALL cp_fm_scale(xsum, xmatrix)

      CALL cp_fm_release(csc)
      CALL cp_fm_release(sc)
      CALL cp_fm_release(work)

      CALL timestop(handle)

   END SUBROUTINE lowdin_orthogonalization

! **************************************************************************************************
!> \brief ...
!> \param gmatrix ...
!> \param vmatrix ...
!> \param xmatrix ...
!> \param rdiag ...
!> \param matrix_s ...
! **************************************************************************************************
   SUBROUTINE wfn_align(gmatrix, vmatrix, xmatrix, rdiag, matrix_s)

      TYPE(cp_fm_type), INTENT(IN)                       :: gmatrix, vmatrix, xmatrix
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: rdiag
      TYPE(dbcsr_type)                                   :: matrix_s

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'wfn_align'
      REAL(KIND=dp), PARAMETER                           :: rone = 1.0_dp, rzero = 0.0_dp

      INTEGER                                            :: handle, n, ncol, ncol_global
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type)                                   :: csc, sc

      CALL timeset(routineN, handle)

      ncol = SIZE(rdiag)
      CALL cp_fm_get_info(matrix=vmatrix, nrow_global=n, ncol_global=ncol_global)
      IF (ncol .GT. ncol_global) CPABORT("Wrong ncol value")

      CALL cp_fm_create(sc, vmatrix%matrix_struct, "SC")
      CALL cp_dbcsr_sm_fm_multiply(matrix_s, vmatrix, sc, ncol)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol, &
                               para_env=vmatrix%matrix_struct%para_env, &
                               context=vmatrix%matrix_struct%context)
      CALL cp_fm_create(csc, fm_struct_tmp, "csc")
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL parallel_gemm('T', 'N', ncol, ncol, n, rone, gmatrix, sc, rzero, csc)
      CALL parallel_gemm('N', 'T', n, ncol, ncol, rone, vmatrix, csc, rzero, sc)
      CALL cp_fm_to_fm(sc, vmatrix, ncol, 1, 1)
      CALL parallel_gemm('N', 'T', n, ncol, ncol, rone, xmatrix, csc, rzero, sc)
      CALL cp_fm_to_fm(sc, xmatrix, ncol, 1, 1)
      !
      CALL lowdin_orthogonalization(vmatrix, xmatrix, ncol, matrix_s)
      !
      CALL cp_dbcsr_sm_fm_multiply(matrix_s, vmatrix, sc, ncol)
      CALL parallel_gemm('T', 'N', ncol, ncol, n, rone, gmatrix, sc, rzero, csc)
      CALL cp_fm_get_diag(csc, rdiag)

      CALL cp_fm_release(csc)
      CALL cp_fm_release(sc)

      CALL timestop(handle)

   END SUBROUTINE wfn_align

! **************************************************************************************************
!> \brief ...
!> \param ematrix ...
!> \param xmatrix ...
!> \param rdiag ...
!> \param matrix_s ...
! **************************************************************************************************
   SUBROUTINE xvec_ovlp(ematrix, xmatrix, rdiag, matrix_s)

      TYPE(cp_fm_type), INTENT(IN)                       :: ematrix, xmatrix
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: rdiag
      TYPE(dbcsr_type)                                   :: matrix_s

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'xvec_ovlp'
      REAL(KIND=dp), PARAMETER                           :: rone = 1.0_dp, rzero = 0.0_dp

      INTEGER                                            :: handle, n, ncol, ncol_global
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type)                                   :: csc, sc

      CALL timeset(routineN, handle)

      ncol = SIZE(rdiag)
      CALL cp_fm_get_info(matrix=xmatrix, nrow_global=n, ncol_global=ncol_global)
      IF (ncol .GT. ncol_global) CPABORT("Wrong ncol value")

      CALL cp_fm_create(sc, xmatrix%matrix_struct, "SC")
      CALL cp_dbcsr_sm_fm_multiply(matrix_s, xmatrix, sc, ncol)

      NULLIFY (fm_struct_tmp)
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol, &
                               para_env=xmatrix%matrix_struct%para_env, &
                               context=xmatrix%matrix_struct%context)
      CALL cp_fm_create(csc, fm_struct_tmp, "csc")
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL parallel_gemm('T', 'N', ncol, ncol, n, rone, ematrix, sc, rzero, csc)
      CALL cp_fm_get_diag(csc, rdiag)

      CALL cp_fm_release(csc)
      CALL cp_fm_release(sc)

      CALL timestop(handle)

   END SUBROUTINE xvec_ovlp

END MODULE qs_tddfpt2_assign
