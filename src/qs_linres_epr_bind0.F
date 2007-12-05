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
  PUBLIC :: chicorr_analytic, mult_G_ov_G2_grid

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
    !gapw=.FALSE.

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
       epr_env%bind0(idir,iB) = epr_env%bind0(idir,iB)+&
            &                        chi(idir)
       IF (output_unit>0) WRITE(output_unit,*) idir, epr_env%bind0(idir,iB)
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

END MODULE qs_linres_epr_bind0



