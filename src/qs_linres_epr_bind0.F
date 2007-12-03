!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****m* cp2k/qs_linres_epr_bind0 *
!!
!!   NAME
!!     qs_linres_epr_bind0
!!
!!   FUNCTION
!!     from the response current density calculates
!!     the susceptibility
!!
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     created 02-2006 [MI]
!!
!!   SOURCE
!!****

MODULE qs_linres_epr_bind0

  USE cell_types,                      ONLY: cell_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_trace
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: gaussi
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_transfer,&
                                             pw_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_epr_op,                ONLY: ind_m2,&
                                             set_vecp
  USE qs_linres_types,                 ONLY: epr_env_type,&
                                             get_epr_env
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_operators_ao,                 ONLY: rRc_xyz_der_ao,&
                                             set_up_op_sm
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             deallocate_matrix_set,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: chicorr_analytic, mult_G_ov_G2_grid !R

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_epr_bind0'

!!***
! **************

CONTAINS

  SUBROUTINE chicorr_analytic(epr_env,qs_env,ispin,iB,error)
    !
    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: ispin, iB
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'chicorr_analytic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, idir, idir2, ii, &
                                                iiB, iii, iiiB, istat, &
                                                istate, nao, nstates(2), &
                                                output_unit
    LOGICAL                                  :: failure, gapw, ionode
    REAL(dp)                                 :: chi(3), contrib, contrib2, &
                                                dk(3)
    REAL(dp), DIMENSION(:), POINTER          :: occupation
    REAL(dp), DIMENSION(:, :), POINTER       :: vecbuf_0, vecbuf_1_d, &
                                                vecbuf_1_rxp
    TYPE(cp_2d_r_p_type)                     :: vecbuf_1_p(3)
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: centers_set
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: psi1_D, psi1_p, psi1_rxp
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff, momxpsi_istate, &
                                                psi0_istate, psi_p1_istate, &
                                                psi_p2_istate, psi_rxp_istate
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, op_mom_ao, op_p_ao
    TYPE(real_matrix_p_type), &
      DIMENSION(:, :), POINTER               :: op_mom_der_ao

!

    CALL timeset(routineN,"I"," ",handle)
    
    failure = .FALSE.
    NULLIFY(dft_control,matrix_s,mos,para_env)
    NULLIFY(mo_coeff,occupation)
    NULLIFY(op_mom_ao,op_mom_der_ao,centers_set)
    NULLIFY(psi0_istate,psi_rxp_istate,psi_p1_istate,psi_p2_istate,momxpsi_istate)
    NULLIFY(op_p_ao,psi1_p,psi1_rxp,psi1_D)
    logger => cp_error_get_logger(error) 
    ionode = logger%para_env%mepos==logger%para_env%source 
    output_unit = -1 
    IF (ionode) THEN 
      output_unit= cp_logger_get_default_unit_nr(logger) 
    END IF 

    CALL get_qs_env(qs_env=qs_env,&
         &          dft_control=dft_control,&
         &          matrix_s=matrix_s,&
         &          mos=mos,&
         &          para_env=para_env,&
         &          error=error)

    gapw = dft_control%qs_control%gapw
    IF (output_unit>0) THEN
       WRITE(output_unit,*) 'WARNING WARNING WARNING WARNING WARNING WARNING WARNING'
       WRITE(output_unit,*) 'chi_analytic: gapw=',gapw,' reset to .false.'
       WRITE(output_unit,*) 'WARNING WARNING WARNING WARNING WARNING WARNING WARNING'
    END IF
    gapw=.FALSE.

    CALL get_epr_env(epr_env=epr_env,&
         &           nao=nao,&
         &           nstates=nstates,&
         &           centers_set=centers_set,&
         &           op_p_ao=op_p_ao,&
         &           psi1_p=psi1_p,&
         &           psi1_rxp=psi1_rxp,&
         &           psi1_D=psi1_D,&
         &           error=error)
    !
    ! Allocate sparse matrices for dipole, quadrupole and their derivatives => 9x3
    ! Remember the derivatives are antisymmetric
    CALL allocate_matrix_set(op_mom_ao,9,error=error)
    CALL allocate_matrix_set(op_mom_der_ao,9,3,error=error)
    CALL set_up_op_sm(op_mom_ao(1)%matrix,qs_env,symmetry="none",&
         &            name="op_mom",error=error)
    DO idir2=1,3
       CALL replicate_matrix_structure(op_mom_ao(1)%matrix, &
            &                          op_mom_der_ao(1,idir2)%matrix,&
            &                          "op_mom_der_ao"//"-"//TRIM(&
            &                          ADJUSTL(cp_to_string(idir2))),&
            &                          target_symmetry="none",&
            &                          error=error)
       CALL set_matrix(op_mom_der_ao(1,idir2)%matrix,0.0_dp)
    ENDDO

    DO idir = 2,SIZE(op_mom_ao,1)
       CALL replicate_matrix_structure(op_mom_ao(1)%matrix, &
            &                          op_mom_ao(idir)%matrix,&
            &                          "op_mom_ao"//"-"//TRIM(&
            &                          ADJUSTL(cp_to_string(idir))),&
            &                          target_symmetry="none",&
            &                          error=error)
       CALL set_matrix(op_mom_ao(idir)%matrix,0.0_dp)
       DO idir2=1,3
          CALL replicate_matrix_structure(op_mom_ao(1)%matrix,&
               &                          op_mom_der_ao(idir,idir2)%matrix,&
               &                          "op_mom_der_ao"//"-"//TRIM(&
               &                          ADJUSTL(cp_to_string(idir*idir2))),&
               &                          target_symmetry="none",&
               &                          error=error)
          CALL set_matrix(op_mom_der_ao(idir,idir2)%matrix,0.0_dp)
       ENDDO
    ENDDO
    !
    ! get iiB and iiiB
    CALL set_vecp(iB,iiB,iiiB)
    !
    ! get ground state MOS
    CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,&
         &          occupation_numbers=occupation)
    !
    ! Allocate full matrices for only one vector
    NULLIFY(tmp_fm_struct)
    CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
         &                   ncol_global=1,para_env=para_env,&
         &                   context=mo_coeff%matrix_struct%context,&
         &                   error=error)
    CALL cp_fm_create(psi0_istate, tmp_fm_struct ,error=error)
    CALL cp_fm_create(psi_rxp_istate, tmp_fm_struct ,error=error)
    CALL cp_fm_create(psi_p1_istate, tmp_fm_struct ,error=error)
    CALL cp_fm_create(psi_p2_istate, tmp_fm_struct ,error=error)
    CALL cp_fm_create(momxpsi_istate, tmp_fm_struct ,error=error)
    CALL cp_fm_struct_release( tmp_fm_struct ,error=error)
    !
    ! Allocate buffer vectors
    NULLIFY(vecbuf_0,vecbuf_1_rxp,vecbuf_1_d)
    ALLOCATE(vecbuf_0(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(vecbuf_1_rxp(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(vecbuf_1_d(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO idir=1,3
       NULLIFY(vecbuf_1_p(idir)%array)
       ALLOCATE(vecbuf_1_p(idir)%array(1,nao),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDDO
    !
    ! Initialize the temporary vector chi
    chi = 0.0_dp
    !
    ! Start loop over the occupied  states
    DO istate = 1,nstates(ispin)
       !
       ! Get the Wannier center of the istate-th ground state orbital
       dk(1:3) = centers_set(ispin)%array(1:3,istate)
       !
       ! Compute the multipole integrals for the state istate,
       ! using as reference center the corresponding Wannier center
       DO idir = 1,9
          CALL set_matrix(op_mom_ao(idir)%matrix, 0.0_dp)
          DO idir2 = 1,3
             CALL set_matrix(op_mom_der_ao(idir,idir2)%matrix, 0.0_dp)
          ENDDO
       ENDDO
       !
       !WRITE(*,*) ' state ', istate, dk(1:3)
       CALL rRc_xyz_der_ao(op_mom_ao,op_mom_der_ao,qs_env,dk,order=2,&
            &              minimum_image=.FALSE.,soft=gapw,error=error)
       !
       ! EXTRACT from the ground state full MOS matrix the coefficients of istate
       CALL cp_fm_get_submatrix(mo_coeff,vecbuf_0,1,istate,nao,1,&
            &                   transpose=.TRUE.,error=error)
       vecbuf_0(1,1:nao) = vecbuf_0(1,1:nao) * occupation(istate)
       CALL cp_fm_set_submatrix(psi0_istate,vecbuf_0,1,1,nao,1,&
            &                   transpose=.TRUE.,error=error)
       !
       ! Extract from the full MOS matrices of the response functions the coefficients
       ! referring to istate
       ! psi1_rxp_iB_istate
       CALL cp_fm_get_submatrix(psi1_rxp(ispin,iB)%matrix,vecbuf_1_rxp,1,&
            &                   istate,nao,1,transpose=.TRUE.,error=error)
       
       IF(epr_env%full_epr) THEN
          ! psi1_d_iB_istate
          CALL cp_fm_get_submatrix(psi1_D(ispin,iB)%matrix,vecbuf_1_d,1,istate,&
               &                   nao,1,transpose=.TRUE.,error=error)
          vecbuf_1_rxp(1,1:nao) = vecbuf_1_rxp(1,1:nao) - vecbuf_1_d(1,1:nao)
       ENDIF
       CALL cp_fm_set_submatrix(psi_rxp_istate,vecbuf_1_rxp,1,1,nao,1,&
            &                   transpose=.TRUE.,error=error)
       !
       ! psi1_p_iiB_istate and psi1_p_iiiB_istate
       CALL cp_fm_get_submatrix(psi1_p(ispin,iiB)%matrix,&
            &                   vecbuf_1_p(iiB)%array,1,istate,nao,1,&
            &                   transpose=.TRUE.,error=error)
       CALL cp_fm_set_submatrix(psi_p1_istate,vecbuf_1_p(iiB)%array,1,1,&
            &                   nao,1,transpose=.TRUE.,error=error)
       CALL cp_fm_get_submatrix(psi1_p(ispin,iiiB)%matrix,&
            &                   vecbuf_1_p(iiiB)%array,1,istate,nao,1,&
            &                   transpose=.TRUE.,error=error)
       CALL cp_fm_set_submatrix(psi_p2_istate,vecbuf_1_p(iiiB)%array,1,1,&
            &                   nao,1,transpose=.TRUE.,error=error)
       !
       ! Multuply left and right by the appropriate coefficients and sum into the
       ! correct component of the chi tensor using the appropriate multiplicative factor
       ! (don't forget the occupation number)
       ! Loop over the cartesian components of the tensor
       ! The loop over the components of the external field is external, thereby
       ! only one column of the chi tensor is computed here
       DO idir = 1,3
          ! get ii and iii
          CALL set_vecp(idir,ii,iii)
          !
          ! term: 2[C0| (r-dk)_ii |d_iii(C1(rxp-D))]-2[C0| (r-dk)_iii |d_ii(C1(rxp-D))]
          ! the factor 2 should be already included in the matrix elements
          contrib = 0.0_dp
          CALL cp_sm_fm_multiply(op_mom_der_ao(ii,iii)%matrix,psi_rxp_istate,&
               &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
          chi(idir) = chi(idir) + 2.0_dp * contrib

          contrib = 0.0_dp
          CALL cp_sm_fm_multiply(op_mom_der_ao(iii,ii)%matrix,psi_rxp_istate,&
               &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
          chi(idir) = chi(idir) - 2.0_dp *contrib
          !
          ! correction: dk_ii*2[C0| d_iii(C1(rxp-D))] - dk_iii*2[C0| d_ii(C1(rxp-D))]
          ! factor 2 not included in the matrix elements
          contrib = 0.0_dp
          CALL cp_sm_fm_multiply(op_p_ao(iii)%matrix,psi_rxp_istate,&
               &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
          chi(idir) = chi(idir) + 2.0_dp*dk(ii)*contrib
          !
          contrib2 = 0.0_dp
          CALL cp_sm_fm_multiply(op_p_ao(ii)%matrix,psi_rxp_istate,&
               &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
          chi(idir) = chi(idir) - 2.0_dp*dk(iii)*contrib2
          !
          ! term: -2[C0| (r-dk)_ii  (r-dk)_iiB | d_iii(C1(piiiB))] \
          !       +2[C0| (r-dk)_iii (r-dk)_iiB | d_ii(C1(piiiB))]
          ! the factor 2 should be already included in the matrix elements
          contrib = 0.0_dp
          idir2 = ind_m2(ii,iiB)
          CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,iii)%matrix,psi_p2_istate,&
               &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
          chi(idir) = chi(idir) - 2.0_dp * contrib
          contrib2 = 0.0_dp
          IF(iiB==iii) THEN
             CALL cp_sm_fm_multiply(op_mom_ao(ii)%matrix,psi_p2_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
             chi(idir) = chi(idir) - contrib2
          ENDIF

          contrib = 0.0_dp
          idir2 = ind_m2(iii,iiB)
          CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,ii)%matrix,psi_p2_istate,&
               &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
          chi(idir) = chi(idir) + 2.0_dp * contrib
          contrib2 = 0.0_dp
          IF(iiB==ii) THEN
             CALL cp_sm_fm_multiply(op_mom_ao(iii)%matrix,psi_p2_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
             chi(idir) = chi(idir) + contrib2
          ENDIF
          !
          ! correction: -dk_ii * 2[C0|(r-dk)_iiB | d_iii(C1(piiiB))] \
          !             +dk_iii * 2[C0|(r-dk)_iiB | d_ii(C1(piiiB))]
          ! the factor 2 should be already included in the matrix elements
          ! no additional correction terms because of the orthogonality between C0 and C1
          contrib = 0.0_dp
          CALL cp_sm_fm_multiply(op_mom_der_ao(iiB,iii)%matrix,psi_p2_istate,&
               &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
          chi(idir) = chi(idir) - 2.0_dp*dk(ii)*contrib
          !
          contrib2 = 0.0_dp
          CALL cp_sm_fm_multiply(op_mom_der_ao(iiB,ii)%matrix,psi_p2_istate,&
               &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
          chi(idir) = chi(idir) + 2.0_dp*dk(iii)*contrib2
          !
          ! term: +2[C0| (r-dk)_ii  (r-dk)_iiiB | d_iii(C1(piiB))] \
          !       -2[C0| (r-dk)_iii (r-dk)_iiiB | d_ii(C1(piiB))]
          ! the factor 2 should be already included in the matrix elements
          contrib = 0.0_dp
          idir2 = ind_m2(ii,iiiB)
          CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,iii)%matrix,psi_p1_istate,&
               &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
          chi(idir) = chi(idir) + 2.0_dp*contrib
          contrib2 = 0.0_dp
          IF(iiiB==iii) THEN
             CALL cp_sm_fm_multiply(op_mom_ao(ii)%matrix,psi_p1_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
             chi(idir) = chi(idir) + contrib2
          ENDIF

          contrib = 0.0_dp
          idir2 = ind_m2(iii,iiiB)
          CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,ii)%matrix,psi_p1_istate,&
               &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
          chi(idir) = chi(idir) - 2.0_dp * contrib
          contrib2 = 0.0_dp
          IF(iiiB==ii) THEN
             CALL cp_sm_fm_multiply(op_mom_ao(iii)%matrix,psi_p1_istate,&
                  &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
             CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
             chi(idir) = chi(idir) - contrib2
          ENDIF
          !
          ! correction: +dk_ii * 2[C0|(r-dk)_iiiB | d_iii(C1(piiB))] +\
          !             -dk_iii * 2[C0|(r-dk)_iiiB | d_ii(C1(piiB))]
          ! the factor 2 should be already included in the matrix elements
          contrib = 0.0_dp
          CALL cp_sm_fm_multiply(op_mom_der_ao(iiiB,iii)%matrix,psi_p1_istate,&
               &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
          chi(idir) = chi(idir) + 2.0_dp*dk(ii)*contrib
          !
          contrib2 = 0.0_dp
          CALL cp_sm_fm_multiply(op_mom_der_ao(iiiB,ii)%matrix,psi_p1_istate,&
               &                 momxpsi_istate,ncol=1,alpha=1.e0_dp,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
          chi(idir) = chi(idir) - 2.0_dp*dk(iii)*contrib2

       ENDDO  ! idir

    ENDDO  ! istate

    IF (output_unit>0) WRITE(output_unit,*) "CHI TENSOR : iB =", iB
    DO idir = 1,3
       epr_env%bind0(idir,iB) = epr_env%bind0(idir,iB)+& !R
            &                        chi(idir) !R *epr_env%chi_factor 
       IF (output_unit>0) WRITE(output_unit,*) idir, epr_env%bind0(idir,iB) !R
    ENDDO
    !
    ! deallocate the sparse matrices
    CALL deallocate_matrix_set(op_mom_ao,error=error)
    CALL deallocate_matrix_set(op_mom_der_ao,error=error)
    CALL cp_fm_release (psi0_istate, error=error )
    CALL cp_fm_release (psi_rxp_istate, error=error )
    CALL cp_fm_release (psi_p1_istate, error=error )
    CALL cp_fm_release (psi_p2_istate, error=error )
    CALL cp_fm_release (momxpsi_istate, error=error )
    DEALLOCATE(vecbuf_0,vecbuf_1_rxp,vecbuf_1_d,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO idir = 1,3
       DEALLOCATE(vecbuf_1_p(idir)%array,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDDO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE chicorr_analytic


!R!  *****************************************************************************************
!R
!R!!****f* qs_linres_epr_shift%interpolate_shift_pwgrid
!R!!
!R!!  NAME
!R!!      epr_interpolate_shift_pwgrid
!R!!
!R!!  FUNCTION
!R!!      interpolate the shift calculated on the PW grid in order to ger
!R!!      the value on arbitrary points in real space
!R!!
!R!!  ARGUMENTS
!R!!      nme_env : to get the shift tensor and the list of additional points
!R!!      particle_set : for the atomic position
!R!!      cell : to take into account the pbs, and to have the volume
!R!!      shift_pw_rspace : specific component of the shift tensor on the pw grid
!R!!      i_B : component of the magnetic field for which the shift is calculated (row)
!R!!      idir : component of the vector \int_{r}[ ((r-r') x j(r))/|r-r'|^3 ] = Bind(r')
!R!!
!R!!  NOTES
!R!!
!R!!  AUTHOR
!R!!    MI
!R!!
!R!!*** **********************************************************************
!R
!R  SUBROUTINE interpolate_shift_pwgrid(epr_env,pw_env,particle_set,cell,shift_pw_rspace,&
!R       &                              i_B,idir,epr_section,error)
!R
!R    TYPE(epr_env_type)                       :: epr_env
!R    TYPE(pw_env_type), POINTER               :: pw_env
!R    TYPE(particle_type), DIMENSION(:), &
!R         POINTER                             :: particle_set
!R    TYPE(cell_type), POINTER                 :: cell
!R    TYPE(pw_p_type)                          :: shift_pw_rspace
!R    INTEGER, INTENT(IN)                      :: i_B, idir
!R    TYPE(section_vals_type), POINTER         :: epr_section
!R    TYPE(cp_error_type), INTENT(INOUT)       :: error
!R
!R    CHARACTER(LEN=*), PARAMETER :: routineN = 'interpolate_shift_pwgrid', &
!R         routineP = moduleN//':'//routineN
!R    INTEGER                                  :: aint_precond, iat, iatom, &
!R                                                max_iter, n_nics, natom, &
!R                                                precond_kind
!R    INTEGER, DIMENSION(:), POINTER           :: cs_atom_list
!R    LOGICAL                                  :: failure, success
!R    REAL(dp)                                 :: eps_r, eps_x, R_iatom(3), &
!R                                                ra(3), shift_val
!R    REAL(dp), DIMENSION(:, :), POINTER       :: r_nics
!R    REAL(dp), DIMENSION(:, :, :), POINTER    :: chemical_shift, &
!R                                                chemical_shift_nics
!R    TYPE(pw_p_type)                          :: shiftspl
!R    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
!R    TYPE(pw_spline_precond_type), POINTER    :: precond
!R    TYPE(section_vals_type), POINTER         :: interp_section
!R
!R#ifdef __parallel
!R    CALL stop_program(routineP,"Doesnt support parallel yet.")
!R#endif
!R
!R    failure = .FALSE.
!R    NULLIFY (interp_section)
!R    NULLIFY (auxbas_pw_pool,precond)
!R    NULLIFY (cs_atom_list,chemical_shift,chemical_shift_nics,r_nics)
!R
!R    CPPrecondition(ASSOCIATED(shift_pw_rspace%pw),cp_failure_level,routineP,error,failure)
!R    IF(.NOT. failure) THEN
!R
!R      interp_section =>  section_vals_get_subs_vals(epr_section,&
!R                        "INTERPOLATOR",error=error)
!R      CALL section_vals_val_get(interp_section,"aint_precond", &
!R                                i_val=aint_precond, error=error)
!R      CALL section_vals_val_get(interp_section,"precond",i_val=precond_kind, error=error)
!R      CALL section_vals_val_get(interp_section,"max_iter",i_val=max_iter, error=error)
!R      CALL section_vals_val_get(interp_section,"eps_r",r_val=eps_r, error=error)
!R      CALL section_vals_val_get(interp_section,"eps_x",r_val=eps_x, error=error)
!R
!R      ! calculate spline coefficients
!R      CALL pw_env_get(pw_env,  auxbas_pw_pool=auxbas_pw_pool,error=error)
!R      CALL pw_pool_create_pw(auxbas_pw_pool,shiftspl%pw, &
!R                             use_data=REALDATA3D,in_space=REALSPACE,&
!R                             error=error)
!R
!R      CALL pw_spline_precond_create(precond,precond_kind=aint_precond,&
!R          pool=auxbas_pw_pool,pbc=.TRUE.,transpose=.FALSE.,error=error)
!R      CALL pw_spline_do_precond(precond,shift_pw_rspace%pw,shiftspl%pw,error=error)
!R      CALL pw_spline_precond_set_kind(precond,precond_kind,error=error)
!R      success=find_coeffs(values=shift_pw_rspace%pw,coeffs=shiftspl%pw,&
!R             linOp=spl3_pbc,preconditioner=precond,pool=auxbas_pw_pool, &
!R             eps_r=eps_r,eps_x=eps_x,max_iter=max_iter, &
!R             print_section=interp_section,error=error)
!R      CPPostconditionNoFail(success,cp_warning_level,routineP,error)
!R      CALL pw_spline_precond_release(precond,error=error)
!R
!R
!R      CALL get_epr_env(epr_env=epr_env, cs_atom_list=cs_atom_list,&
!R           chemical_shift=chemical_shift,&
!R           chemical_shift_nics=chemical_shift_nics,&
!R           n_nics=n_nics,r_nics=r_nics,error=error)
!R
!R      IF(ASSOCIATED(cs_atom_list)) THEN
!R        natom = SIZE(cs_atom_list,1)
!R      ELSE
!R        natom = -1
!R      ENDIF
!R
!R      DO iat = 1,natom
!R         iatom = cs_atom_list(iat)
!R         R_iatom = pbc(particle_set(iatom)%r,cell)
!R         shift_val = Eval_Interp_Spl3_pbc(R_iatom,shiftspl%pw,error=error)
!R         chemical_shift(idir,i_B,iatom)= chemical_shift(idir,i_B,iatom)+&
!R                                         epr_env%shift_factor*shift_val
!R      END DO
!R
!R      IF(epr_env%do_nics) THEN
!R        DO iatom = 1,n_nics
!R          ra(1:3) = r_nics(1:3,iatom)
!R          R_iatom = pbc(ra,cell)
!R          shift_val = Eval_Interp_Spl3_pbc(R_iatom,shiftspl%pw,error)
!R          chemical_shift_nics(idir,i_B,iatom)= chemical_shift_nics(idir,i_B,iatom)+&
!R                                               epr_env%shift_factor*shift_val
!R        END DO
!R      END IF
!R
!R      CALL pw_pool_give_back_pw(auxbas_pw_pool,shiftspl%pw,error=error)
!R
!R    END IF
!R
!R  END SUBROUTINE interpolate_shift_pwgrid
!R
!R  SUBROUTINE gsum_shift_pwgrid(epr_env,pw_env,particle_set,cell,shift_pw_gspace,&
!R       &                       i_B,idir,epr_section,error)
!R    TYPE(epr_env_type)                       :: epr_env
!R    TYPE(pw_env_type), POINTER               :: pw_env
!R    TYPE(particle_type), DIMENSION(:), &
!R      POINTER                                :: particle_set
!R    TYPE(cell_type), POINTER                 :: cell
!R    TYPE(pw_p_type)                          :: shift_pw_gspace
!R    INTEGER, INTENT(IN)                      :: i_B, idir
!R    TYPE(section_vals_type), POINTER         :: epr_section
!R    TYPE(cp_error_type), INTENT(INOUT)       :: error
!R
!R    CHARACTER(LEN=*), PARAMETER :: routineN = 'gsum_shift_pwgrid', &
!R         routineP = moduleN//':'//routineN
!R
!R    INTEGER                                  :: iat, iatom, n_nics, natom
!R    INTEGER, DIMENSION(:), POINTER           :: cs_atom_list
!R    LOGICAL                                  :: failure
!R    REAL(dp)                                 :: R_iatom(3), ra(3)
!R    REAL(dp), DIMENSION(:, :), POINTER       :: r_nics
!R    REAL(dp), DIMENSION(:, :, :), POINTER    :: chemical_shift, &
!R                                                chemical_shift_nics
!R
!R    !
!R
!R    failure = .FALSE.
!R    NULLIFY(cs_atom_list,chemical_shift,chemical_shift_nics,r_nics)
!R    CPPrecondition(ASSOCIATED(shift_pw_gspace%pw),cp_failure_level,routineP,error,failure)
!R    IF(.NOT. failure) THEN
!R       !
!R       CALL get_epr_env(epr_env=epr_env, cs_atom_list=cs_atom_list,&
!R            &           chemical_shift=chemical_shift,&
!R            &           chemical_shift_nics=chemical_shift_nics,&
!R            &           n_nics=n_nics,r_nics=r_nics,error=error)
!R       !
!R       IF(ASSOCIATED(cs_atom_list)) THEN
!R          natom = SIZE(cs_atom_list,1)
!R       ELSE
!R          natom = -1
!R       ENDIF
!R       !
!R       ! compute the chemical shift
!R       DO iat = 1,natom
!R          iatom = cs_atom_list(iat)
!R          R_iatom = pbc(particle_set(iatom)%r,cell)
!R          chemical_shift(idir,i_B,iatom)= chemical_shift(idir,i_B,iatom)+&
!R               & epr_env%shift_factor*REAL(gsumr(R_iatom,shift_pw_gspace%pw),dp)
!R       ENDDO
!R       !
!R       ! compute nics
!R       IF(epr_env%do_nics) THEN
!R          DO iat = 1,n_nics
!R             ra = pbc(r_nics(:,iat),cell)
!R             chemical_shift_nics(idir,i_B,iat)= chemical_shift_nics(idir,i_B,iat)+&
!R                  & epr_env%shift_factor*REAL(gsumr(ra,shift_pw_gspace%pw),dp)
!R          ENDDO
!R       ENDIF
!R    ENDIF
!R  END SUBROUTINE gsum_shift_pwgrid
!R
!R  FUNCTION gsumr(r,pw)
!R    REAL(dp), INTENT(IN)        :: r(3)
!R    TYPE(pw_type), POINTER      :: pw
!R    TYPE(pw_grid_type), POINTER :: grid
!R    COMPLEX(dp)                 :: gsumr
!R
!R    COMPLEX(dp)                 :: rg
!R    INTEGER                     :: ig
!R    TYPE(pw_grid_type), POINTER              :: grid 
!R   
!R    grid => pw%pw_grid
!R    gsumr = CMPLX(0.0_dp,0.0_dp,KIND=dp)
!R    DO ig = grid%first_gne0,grid%ngpts_cut_local
!R       rg = ( grid%g(1,ig)*r(1) + grid%g(2,ig)*r(2) + grid%g(3,ig)*r(3) ) * gaussi
!R       gsumr = gsumr + pw%cc(ig)*EXP(rg)
!R    ENDDO
!R    IF(grid%have_g0) gsumr = gsumr + pw%cc(1)
!R    CALL mp_sum(gsumr,grid%para%group)
!R  END FUNCTION gsumr

!!****f* qs_linres_epr_shift%mult_G_ov_G2_grid
!!
!!  NAME
!!      mult_G_ov_G2_grid
!!
!!  FUNCTION
!!      Given the current density on the PW grid in reciprcal space
!!      (obtained by FFT), calculate the integral
!!        \int_{r}[ ((r-r') x j(r))/|r-r'|^3 ] = Bind(r')
!!      which in reciprcal space reads  (for G/=0)
!!         i G/|G|^2 x J(G)
!!
!!  ARGUMENTS
!!
!!  NOTES
!!     The G=0 component is not comnputed here, but can be evaluated
!!     through the susceptibility and added to the shift in a second time
!!
!!     This method would not work for a non periodic system
!!     It should be generalized like the calculation of Hartree
!!
!!  AUTHOR
!!    MI
!!
!!*** **********************************************************************


  SUBROUTINE  mult_G_ov_G2_grid(cell,pw_pool,rho_gspace,funcG_times_rho,idir,my_chi,error)


    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(pw_p_type), POINTER                 :: rho_gspace
    TYPE(pw_p_type)                          :: funcG_times_rho
    INTEGER, INTENT(IN)                      :: idir
    REAL(dp), INTENT(IN)                     :: my_chi
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mult_G_ov_G2_grid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ig, ng
    LOGICAL                                  :: failure
    REAL(dp)                                 :: g2
    TYPE(pw_grid_type), POINTER              :: grid
    TYPE(pw_type), POINTER                   :: frho, influence_fn

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)

    CALL pw_pool_create_pw ( pw_pool, influence_fn,&
               use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE ,error=error)

    grid => influence_fn % pw_grid
    DO ig = grid % first_gne0, grid % ngpts_cut_local
       g2 = grid % gsq ( ig )
       influence_fn%cc(ig) = gaussi * grid % g(idir,ig)/g2
    END DO  ! ig
    IF ( grid % have_g0 )  influence_fn%cc ( 1 ) = 0.0_dp

    frho => funcG_times_rho%pw
    CALL pw_transfer (rho_gspace%pw,frho)

    ng = SIZE(grid % gsq)
    frho%cc(1:ng) = frho%cc(1:ng)*influence_fn % cc  ( 1 : ng )
    IF ( grid % have_g0 ) frho%cc(1) = my_chi

    CALL pw_pool_give_back_pw(pw_pool,influence_fn,&
         accept_non_compatible=.TRUE.,error=error)

  END SUBROUTINE  mult_G_ov_G2_grid

!R! *****************************************************************************
!R!!****f* qs_linres_epr_shift%epr_print_shift
!R!!
!R!!  NAME
!R!!      epr_print_shift
!R!!
!R!!  FUNCTION
!R!!      Shielding tensor and Chi are printed into a file
!R!!      if required from input
!R!!      It is possible to print only for a subset of atoms or
!R!!      or points in non-ionic positions
!R!!
!R!!  ARGUMENTS
!R!!
!R!!  NOTES
!R!!
!R!!
!R!!  AUTHOR
!R!!    MI
!R!!
!R!!*** **********************************************************************
!R
!R
!R  SUBROUTINE epr_print_shift(epr_env,qs_env,epr_section,error)
!R    TYPE(epr_env_type)                       :: epr_env
!R    TYPE(qs_environment_type), POINTER       :: qs_env
!R    TYPE(section_vals_type), POINTER         :: epr_section
!R    TYPE(cp_error_type), INTENT(INOUT)       :: error
!R
!R    CHARACTER(len=*), PARAMETER :: routineN = 'epr_print_shift', &
!R         routineP = moduleN//':'//routineN
!R
!R    CHARACTER(LEN=2)                         :: element_symbol
!R    CHARACTER(LEN=default_string_length)     :: name, title
!R    INTEGER                                  :: i, iatom, ir, istat, j, &
!R                                                n_nics, nat_print, natom, &
!R                                                output_unit, unit_atoms, &
!R                                                unit_nics
!R    INTEGER, DIMENSION(:), POINTER           :: cs_atom_list
!R    LOGICAL                                  :: failure, gapw, ionode
!R    REAL(dp) :: chi_aniso, chi_iso, chi_sym(3,3), chi_sym_tot(3,3), &
!R      chi_tensor(3,3), chi_tensor_loc(3,3), chi_tmp(3,3), eigenv(3), &
!R      eigenv_tmp(3), rpos(3), shift_nosym(3,3), shift_sym(3,3), shift_tmp(3,3)
!R    REAL(dp), DIMENSION(:, :), POINTER       :: eigenv_chi, r_nics
!R    REAL(dp), DIMENSION(:, :, :), POINTER    :: chemical_shift, &
!R      chemical_shift_loc, chemical_shift_loc_nics, chemical_shift_nics, &
!R      chemical_shift_tot, chemical_shift_nics_tot
!R    REAL(dp), EXTERNAL                       :: DDOT
!R    TYPE(atomic_kind_type), POINTER          :: atom_kind
!R    TYPE(cp_logger_type), POINTER            :: logger
!R    TYPE(dft_control_type), POINTER          :: dft_control
!R    TYPE(particle_type), DIMENSION(:), &
!R      POINTER                                :: particle_set
!R
!R    failure = .FALSE.
!R
!R    NULLIFY (chemical_shift, chemical_shift_nics, r_nics)
!R    NULLIFY (chemical_shift_loc, chemical_shift_loc_nics)
!R    NULLIFY (logger, particle_set,atom_kind,dft_control)
!R    NULLIFY (eigenv_chi)
!R    logger => cp_error_get_logger(error)
!R    ionode = logger%para_env%mepos==logger%para_env%source 
!R    output_unit = -1 
!R    IF (ionode) THEN 
!R    output_unit = cp_logger_get_default_unit_nr(logger)
!R    END IF
!R
!R    CALL get_epr_env(epr_env, chemical_shift=chemical_shift, &
!R         chemical_shift_nics=chemical_shift_nics, &
!R         chemical_shift_loc=chemical_shift_loc, &
!R         chemical_shift_loc_nics=chemical_shift_loc_nics, &
!R         chi_tensor=chi_tensor, chi_tensor_loc=chi_tensor_loc,&
!R         cs_atom_list=cs_atom_list, n_nics=n_nics, r_nics=r_nics,error=error)
!R
!R    CALL get_qs_env(qs_env, dft_control=dft_control, particle_set=particle_set,error=error)
!R    natom = SIZE(particle_set,1)
!R    gapw = dft_control%qs_control%gapw
!R    nat_print=SIZE(cs_atom_list,1)
!R
!R    ALLOCATE(eigenv_chi(3,natom), STAT=istat)
!R    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R    ALLOCATE(chemical_shift_tot(3,3,nat_print),STAT=istat)
!R    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R    IF(epr_env%do_nics) THEN
!R       ALLOCATE(chemical_shift_nics_tot(3,3,n_nics),STAT=istat)
!R       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R    ENDIF
!R    ! Finalize Chi calculation
!R    IF(gapw) THEN
!R       ! Symmetrize
!R       DO i = 1,3
!R          DO j = 1,3
!R             chi_sym(i,j) = (chi_tensor(i,j)+chi_tensor(j,i))/2.0_dp
!R             chi_sym_tot(i,j) = (chi_tensor(i,j)+chi_tensor(j,i)+&
!R                                 chi_tensor_loc(i,j)+chi_tensor_loc(j,i))/2.0_dp
!R          END DO
!R       END DO
!R
!R    ELSE
!R       ! Symmetrize
!R       DO i = 1,3
!R          DO j = 1,3
!R             chi_sym(i,j) = (chi_tensor(i,j)+chi_tensor(j,i))/2.0_dp
!R          END DO
!R       END DO
!R       CALL  dcopy(9,chi_sym,1,chi_sym_tot,1)
!R    END IF
!R    CALL  dcopy(9,chi_sym_tot,1,chi_tmp,1)
!R    CALL  diamat_all(chi_tmp,eigenv)
!R    chi_iso   = (eigenv(1)+ eigenv(2)+eigenv(3))/3.e0_dp
!R    chi_aniso =  eigenv(3)-(eigenv(2)+eigenv(1))/2.e0_dp
!R    !
!R    IF(output_unit > 0) THEN
!R       WRITE(output_unit,'(T2,A,E24.16)') 'CheckSum Chi =',&
!R            & SQRT(DDOT(9,chi_tensor(1,1),1,chi_tensor(1,1),1))
!R    ENDIF
!R    !
!R    IF (BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
!R              "PRINT%CHI_TENSOR",error=error),cp_p_file)) THEN
!R
!R       unit_atoms=cp_print_key_unit_nr(logger,epr_section,"PRINT%CHI_TENSOR",&
!R                  extension=".data",middle_name="CHI",log_filename=.FALSE.,&
!R                  error=error)
!R
!R       WRITE(title,'(A)') "Magnetic Susceptibility Tensor "
!R       IF(unit_atoms > 0) THEN
!R          WRITE(unit_atoms,'(T2,A)') title
!R          IF(gapw) THEN
!R             WRITE(unit_atoms,'(T30,A)') " CHI from SOFT J in 10^-30 J/T^2 units"
!R             WRITE(unit_atoms,'(3(A,f15.6))') ' XX=',chi_tensor(1,1),&
!R                  ' XY=',chi_tensor(1,2),' XZ=',chi_tensor(1,3)
!R             WRITE(unit_atoms,'(3(A,f15.6))') ' YX=',chi_tensor(2,1),&
!R                  ' YY=',chi_tensor(2,2),' YZ=',chi_tensor(2,3)
!R             WRITE(unit_atoms,'(3(A,f15.6))') ' ZX=',chi_tensor(3,1),&
!R                  ' ZY=',chi_tensor(3,2),' ZZ=',chi_tensor(3,3)
!R             WRITE(unit_atoms,'(T30,A)') " CHI from LOCAL J in 10^-30 J/T^2 units"
!R             WRITE(unit_atoms,'(3(A,f15.6))') ' XX=',chi_tensor_loc(1,1),&
!R                  ' XY=',chi_tensor_loc(1,2),' XZ=',chi_tensor_loc(1,3)
!R             WRITE(unit_atoms,'(3(A,f15.6))') ' YX=',chi_tensor_loc(2,1),&
!R                  ' YY=',chi_tensor_loc(2,2),' YZ=',chi_tensor_loc(2,3)
!R             WRITE(unit_atoms,'(3(A,f15.6))') ' ZX=',chi_tensor_loc(3,1),&
!R                 ' ZY=',chi_tensor_loc(3,2),' ZZ=',chi_tensor_loc(3,3)
!R          ELSE
!R             WRITE(unit_atoms,'(T30,A)') " CHI from J in atomic units"
!R             WRITE(unit_atoms,'(3(A,f15.6))')' XX=',chi_tensor(1,1),&
!R                  ' XY=',chi_tensor(1,2),' XZ=',chi_tensor(1,3)
!R             WRITE(unit_atoms,'(3(A,f15.6))')' YX=',chi_tensor(2,1),&
!R                  ' YY=',chi_tensor(2,2),' YZ=',chi_tensor(2,3)
!R             WRITE(unit_atoms,'(3(A,f15.6))')' ZX=',chi_tensor(3,1),&
!R                  ' ZY=',chi_tensor(3,2),' ZZ=',chi_tensor(3,3)
!R          ENDIF
!R          WRITE (unit_atoms,'(/T2,3(A,F15.5))')' PV1=',eigenv(1),&
!R               ' PV2=',eigenv(2),' PV3=',eigenv(3)
!R          WRITE (unit_atoms,'(T2,A,F15.5,10X,A,F15.5)')&
!R               '  iso=',chi_iso,  'aniso=',chi_aniso
!R          CALL dscal(3*3,epr_env%chi_SI2ppmcgs,chi_sym_tot,1)
!R          WRITE(unit_atoms,'(T30,A)') " Total CHI in ppm-cgs units"
!R          WRITE(unit_atoms,'(3(A,f15.6))')' XX=',chi_sym_tot(1,1),&
!R               ' XY=',chi_sym_tot(1,2),' XZ=',chi_sym_tot(1,3)
!R          WRITE(unit_atoms,'(3(A,f15.6))')' YX=',chi_sym_tot(2,1),&
!R               ' YY=',chi_sym_tot(2,2),' YZ=',chi_sym_tot(2,3)
!R          WRITE(unit_atoms,'(3(A,f15.6))')' ZX=',chi_sym_tot(3,1),&
!R               ' ZY=',chi_sym_tot(3,2),' ZZ=',chi_sym_tot(3,3)
!R          WRITE (unit_atoms, '(/T2,3(A,F15.5))') &
!R               ' PV1=',epr_env%chi_SI2ppmcgs*eigenv(1),&
!R               ' PV2=',epr_env%chi_SI2ppmcgs*eigenv(2),&
!R               ' PV3=',epr_env%chi_SI2ppmcgs*eigenv(3)
!R          WRITE (unit_atoms,'(T2,A,F15.5,10X,A,F15.5)')&
!R               '  iso=',epr_env%chi_SI2ppmcgs*chi_iso, &
!R               'aniso=',epr_env%chi_SI2ppmcgs*chi_aniso
!R       ENDIF
!R       
!R       CALL cp_print_key_finished_output(unit_atoms, logger,epr_section,&
!R            &                            "PRINT%CHI_TENSOR", error=error)
!R    ENDIF ! print chi
!R    !
!R    ! Add the chi part to the shifts
!R    chemical_shift_tot = 0.0_dp
!R    DO ir = 1,nat_print
!R       iatom = cs_atom_list(ir)
!R       rpos(1:3) = particle_set(iatom)%r(1:3)
!R       atom_kind => particle_set(iatom)%atomic_kind
!R       CALL get_atomic_kind(atom_kind,name=name,element_symbol=element_symbol)
!R       DO i = 1,3
!R          DO j = 1,3
!R             ! here we need to multiply chi and the shift by some factors
!R             chemical_shift_tot(i,j,ir) = chi_tensor(i,j) * epr_env%chi_SI2shiftppm &
!R                  & + chemical_shift(i,j,iatom) * twopi**2
!R             IF(gapw) THEN
!R                !vw there is a sign problem here wv!
!R                chemical_shift_tot(i,j,ir) = chemical_shift_tot(i,j,ir)-chemical_shift_loc(i,j,iatom)
!R             ENDIF
!R!write(*,*) 'soft',chi_tensor(i,j) * epr_env%chi_SI2shiftppm &
!R!                  & + chemical_shift(i,j,iatom) * twopi**2,&
!R!           'locl',chemical_shift_loc(i,j,iatom)
!R          ENDDO
!R       ENDDO
!R    END DO  ! ir
!R    IF(output_unit > 0) THEN
!R       WRITE(output_unit,'(T2,A,E24.16)') 'CheckSum Shifts =',&
!R            & SQRT(DDOT(9*SIZE(chemical_shift_tot,3),chemical_shift_tot(1,1,1),1,&
!R            &           chemical_shift_tot(1,1,1),1))
!R    ENDIF
!R    !
!R    ! print shifts
!R    IF (BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
!R              "PRINT%SHIELDING_TENSOR",error=error),cp_p_file)) THEN
!R
!R       unit_atoms=cp_print_key_unit_nr(logger,epr_section,"PRINT%SHIELDING_TENSOR",&
!R            &                          extension=".data",middle_name="SHIFT",&
!R            &                          log_filename=.FALSE.,error=error)
!R
!R       nat_print = SIZE(cs_atom_list,1)
!R       IF(unit_atoms > 0) THEN
!R          WRITE(title,'(A,1X,I5)') "Shielding atom at atomic positions. # tensors printed ", nat_print
!R          WRITE(unit_atoms,'(T2,A)') "WARNING! SHIFTS ARE MULTIPLIED BY A MAGIC CONSTANT [4*Pi^2]!"
!R          WRITE(unit_atoms,'(T2,A)') "WARNING! THE SHIFTS CONTAIN THE CHI CORRECTION!"
!R          WRITE(unit_atoms,'(T2,A)') title
!R          DO ir = 1,nat_print
!R             iatom = cs_atom_list(ir)
!R             rpos(1:3) = particle_set(iatom)%r(1:3)
!R             atom_kind => particle_set(iatom)%atomic_kind
!R             CALL get_atomic_kind(atom_kind,name=name,element_symbol=element_symbol)
!R             DO i = 1,3
!R                DO j = 1,3
!R                   ! here we need to multiply chi and the shift by some factors
!R                   shift_nosym(i,j) = chi_tensor(i,j) * epr_env%chi_SI2shiftppm &
!R                        &           + chemical_shift(i,j,iatom) * twopi**2
!R                   shift_sym(i,j) = chi_sym(i,j)+(chemical_shift(i,j,iatom)+&
!R                        &                         chemical_shift(j,i,iatom))/2.0_dp
!R                   IF(gapw) THEN
!R                     !vw there is a sign problem here wv!
!R                     shift_nosym(i,j) = shift_nosym(i,j)-chemical_shift_loc(i,j,iatom)
!R                     shift_sym(i,j) = shift_sym(i,j)-(chemical_shift_loc(i,j,iatom)+&
!R                                      chemical_shift_loc(j,i,iatom))/2.0_dp
!R                   ENDIF
!R                ENDDO
!R             ENDDO
!R             CALL  dcopy(9,shift_sym,1,shift_tmp,1)
!R             CALL  diamat_all(shift_tmp,eigenv_tmp)
!R             eigenv_chi(1:3, iatom) = eigenv_tmp(1:3)
!R             WRITE(unit_atoms,'(T2,I5,A,2X,A2,2X,3f15.6)') iatom,name,element_symbol,rpos(1:3)
!R             ! here we set eigenv_tmp to huge... need to be consistent with the units used above.
!R             WRITE(unit_atoms,'(A,3f15.6,5x,f15.6)') 'S_x = ',shift_nosym(1, 1:3), HUGE(0.0_dp)!eigenv_tmp(1)
!R             WRITE(unit_atoms,'(A,3f15.6,5x,f15.6)') 'S_y = ',shift_nosym(2, 1:3), HUGE(0.0_dp)!eigenv_tmp(2)
!R             WRITE(unit_atoms,'(A,3f15.6,5x,f15.6)') 'S_z = ',shift_nosym(3, 1:3), HUGE(0.0_dp)!eigenv_tmp(3)
!R          ENDDO  ! ir
!R       ENDIF
!R       CALL cp_print_key_finished_output(unit_atoms,logger,epr_section,&
!R            &                            "PRINT%SHIELDING_TENSOR",error=error)
!R
!R       IF(epr_env%do_nics) THEN
!R          !
!R          ! Add the chi part to the nics
!R          chemical_shift_nics_tot = 0.0_dp
!R          DO ir = 1,n_nics
!R             DO i = 1,3
!R                DO j = 1,3
!R                   ! here we need to multiply chi and the shift by some factors
!R                   chemical_shift_nics_tot(i,j,ir) = chi_tensor(i,j) * epr_env%chi_SI2shiftppm &
!R                        &                          + chemical_shift_nics(i,j,ir) * twopi**2
!R                   IF(gapw) THEN
!R                      CALL stop_program(routineP,"GAPW and NICS doesnt work yet.")
!R                   ENDIF
!R                ENDDO
!R             ENDDO
!R          END DO  ! ir
!R          IF(output_unit > 0) THEN
!R             WRITE(output_unit,'(T2,A,E24.16)') 'CheckSum NICS =',&
!R                  & SQRT(DDOT(9*SIZE(chemical_shift_nics_tot,3),chemical_shift_nics_tot(1,1,1),&
!R                  &           1,chemical_shift_nics_tot(1,1,1),1))
!R          ENDIF
!R          !
!R          unit_nics=cp_print_key_unit_nr(logger,epr_section,"PRINT%SHIELDING_TENSOR",&
!R               &                         extension=".data",middle_name="NICS",&
!R               &                         log_filename=.FALSE.,error=error)
!R          IF(unit_nics > 0) THEN
!R             WRITE(title,'(A,1X,I5)') "Shielding atom at atomic positions. # tensors printed ", n_nics
!R             WRITE(unit_nics,'(T2,A)') title
!R             DO ir = 1,n_nics
!R                WRITE(unit_nics,'(T2,I5,2X,3f15.6)') ir,r_nics(1:3,ir)
!R                IF(gapw) THEN
!R                   CALL stop_program(routineP,"GAPW and NICS doesnt work yet.")
!R                   !WRITE(unit_nics,'(3f15.6,5x,3f15.6)') chemical_shift_nics(1, 1:3, iatom),&
!R                   !     chemical_shift_loc_nics(1,1:3, iatom)
!R                   !WRITE(unit_nics,'(3f15.6,5x,3f15.6)') chemical_shift_nics(2, 1:3, iatom),&
!R                   !     chemical_shift_loc_nics(2,1:3, iatom)
!R                   !WRITE(unit_nics,'(3f15.6,5x,3f15.6)') chemical_shift_nics(3, 1:3, iatom),&
!R                   !     chemical_shift_loc_nics(3,1:3, iatom)
!R                ELSE
!R                   WRITE(unit_nics,'(3f15.6)') chemical_shift_nics_tot(1, 1:3, ir)
!R                   WRITE(unit_nics,'(3f15.6)') chemical_shift_nics_tot(2, 1:3, ir)
!R                   WRITE(unit_nics,'(3f15.6)') chemical_shift_nics_tot(3, 1:3, ir)
!R                ENDIF
!R             ENDDO
!R          ENDIF
!R          CALL cp_print_key_finished_output(unit_nics,logger,epr_section,&
!R               &                            "PRINT%SHIELDING_TENSOR",error=error)
!R       ENDIF
!R    ENDIF  ! print shift
!R    !
!R    ! clean up
!R    DEALLOCATE(eigenv_chi,STAT=istat)
!R    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R    DEALLOCATE(chemical_shift_tot, STAT=istat)
!R    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R    IF(epr_env%do_nics) THEN
!R       DEALLOCATE(chemical_shift_nics_tot,STAT=istat)
!R       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
!R    ENDIF
!R    !
!R  END SUBROUTINE epr_print_shift

END MODULE qs_linres_epr_bind0



