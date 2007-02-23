
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****m* cp2k/qs_linres_epr_bind0 *
!!
!!   NAME 
!!     qs_linres_epr_current
!!
!!   FUNCTION
!!     from the response current density calculates the shift tensor
!!     and the susceptibility
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

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE coefficient_types,               ONLY: coeff_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_trace
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply,&
                                             copy_sm_to_fm 
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_get_rval,&
                                             section_get_lval,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: gaussi
  USE mathlib,                         ONLY: diamat_all
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_coeff,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_spline_utils,                 ONLY: Eval_Interp_Spl3_pbc,&
                                             find_coeffs,&
                                             pw_spline_do_precond,&
                                             pw_spline_precond_create,&
                                             pw_spline_precond_release,&
                                             pw_spline_precond_set_kind,&
                                             pw_spline_precond_type,&
                                             spl3_pbc
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_transfer,&
                                             pw_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_epr_op,                ONLY: ind_m2,&
                                             set_vecp
  USE qs_linres_types,                 ONLY: get_epr_env,&
                                             epr_env_type
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
  PUBLIC :: bind0_soft_analytic, mult_G_ov_G2_grid !R

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_epr_shift'

!!***
! **************

CONTAINS

  SUBROUTINE bind0_soft_analytic(epr_env,qs_env,ispin,iB,error)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: ispin, iB
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'bind0_soft_analytic', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, idir, idir2, ii, iiB, iii, iiiB, istat, istate, nao, nstates(2)
    LOGICAL :: failure, gapw
    REAL(dp)      :: chi(3), contrib, contrib2, dk(3)
    REAL(dp), DIMENSION(:), POINTER          :: occupation
    REAL(dp), DIMENSION(:,:), POINTER        :: vecbuf_0
    REAL(dp), DIMENSION(:,:), POINTER        :: vecbuf_1_rxp
    REAL(dp), DIMENSION(:,:), POINTER        :: vecbuf_1_d
    TYPE(cp_2d_r_p_type)                     :: vecbuf_1_p(3)
    TYPE(cp_2d_r_p_type), DIMENSION(:),&
      POINTER                                :: centers_set
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER :: fm_work2, mo_coeff, momxpsi_istate, &
      psi0_istate, psi_p1_istate, psi_p2_istate, psi_rxp_istate
    TYPE(cp_fm_p_type), DIMENSION(:,:), &
      POINTER                                :: psi1_p, psi1_rxp, psi1_D
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, op_mom_ao
    TYPE(real_matrix_p_type), DIMENSION(:,:), &
      POINTER                                :: op_mom_der_ao
!dbg
!     INTEGER :: fnum, i, j
!     TYPE(cp_fm_type), POINTER :: fm_work2
!dbg

    CALL timeset(routineN,"I"," ",handle)

    failure = .FALSE.

    NULLIFY(dft_control,matrix_s,mos,para_env)
    NULLIFY(mo_coeff,occupation)
    NULLIFY(op_mom_ao,op_mom_der_ao,centers_set)
    NULLIFY(psi0_istate,psi_rxp_istate,psi_p1_istate,psi_p2_istate,momxpsi_istate)
    NULLIFY(psi1_p, psi1_rxp, psi1_D)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,matrix_s=matrix_s,&
         mos=mos,&
         para_env=para_env,&
         error=error)
    gapw = dft_control%qs_control%gapw

    CALL get_epr_env(epr_env=epr_env, nao=nao, nstates=nstates, &
                     centers_set=centers_set,&
                     psi1_p=psi1_p, psi1_rxp=psi1_rxp, psi1_D=psi1_D  )

!   Allocate sparse matric for dipole, quadrupole and their derivatives => 9x3
!   Remember the derivatives are antisymmetric
    CALL allocate_matrix_set(op_mom_ao,9,error=error)
    CALL allocate_matrix_set(op_mom_der_ao,9,3,error=error)
    CALL set_up_op_sm(op_mom_ao(1)%matrix,qs_env,symmetry="none",&
          name="op_mom",error=error)    
    DO idir2=1,3
       CALL replicate_matrix_structure(op_mom_ao(1)%matrix, &
            op_mom_der_ao(1,idir2)%matrix,&
            "op_mom_der_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir2))),&
            target_symmetry="none",error=error)
       CALL set_matrix(op_mom_der_ao(1,idir2)%matrix,0.0_dp)
    END DO

    DO idir = 2,SIZE(op_mom_ao,1)
      CALL replicate_matrix_structure(op_mom_ao(1)%matrix, &
           op_mom_ao(idir)%matrix,&
           "op_mom_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir))),&
           target_symmetry="none",error=error)
      CALL set_matrix(op_mom_ao(idir)%matrix,0.0_dp)
      DO idir2=1,3
        CALL replicate_matrix_structure(op_mom_ao(1)%matrix, &
             op_mom_der_ao(idir,idir2)%matrix,&
             "op_mom_der_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir*idir2))),&
             target_symmetry="none",error=error)
        CALL set_matrix(op_mom_der_ao(idir,idir2)%matrix,0.0_dp)
      END DO
    END DO

!   get iiB and iiiB
    CALL set_vecp(iB,iiB,iiiB)

!   get ground state MOS
    CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,&
                    occupation_numbers=occupation)

!   Allocate full matrices for only one vector
    NULLIFY(tmp_fm_struct)
    CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
         ncol_global=1,para_env=para_env,context=mo_coeff%matrix_struct%context,error=error)
    CALL cp_fm_create (psi0_istate, tmp_fm_struct, error=error )
    CALL cp_fm_create (psi_rxp_istate, tmp_fm_struct, error=error )
    CALL cp_fm_create (psi_p1_istate, tmp_fm_struct, error=error )
    CALL cp_fm_create (psi_p2_istate, tmp_fm_struct, error=error )
    CALL cp_fm_create (momxpsi_istate, tmp_fm_struct, error=error )
    CALL cp_fm_struct_release ( tmp_fm_struct, error=error )
!dbg
!    CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
!          ncol_global=nao,para_env=para_env,context=mo_coeff%matrix_struct%context)
!    CALL cp_fm_create (fm_work2, tmp_fm_struct )
!    CALL cp_fm_struct_release ( tmp_fm_struct )
!dbg
  

!   Allocate buffer vectors
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
    END DO 
    
!   Initialize the temporary vector chi
    chi = 0.0_dp

!   Start loop over the occupied  states
    DO istate = 1,nstates(ispin)

!     Get the Wannier center of the istate-th ground state orbital
      dk(1:3) = centers_set(ispin)%array(1:3,istate)

!     Compute the multipole integrals for the state istate,
!     using as reference center the corresponding Wannier center
      DO idir = 1,9
         CALL set_matrix(op_mom_ao(idir)%matrix, 0.0_dp)
         DO idir2 = 1,3
           CALL set_matrix(op_mom_der_ao(idir,idir2)%matrix, 0.0_dp)
         END DO
      END DO
!dbg
!  write(*,*) ' state ', istate, dk(1:3)
      CALL rRc_xyz_der_ao(op_mom_ao,op_mom_der_ao,qs_env,dk,order=2,&
           minimum_image=.TRUE.,soft=gapw,error=error)
!  stop 'check derivatives' 
!dbg
!dbg
!       DO idir= 1,9      
!         fnum = 100+idir
!         CALL copy_sm_to_fm(op_mom_ao(idir)%matrix,fm_work2,error=error)
!         write(fnum,*) ' op_mom_ao', idir
!         DO  i=1,nao
!           write(fnum,'(i5,50f10.6)') i, (fm_work2%local_data(i,j),j=1,nao)
!         END DO
!         DO idir2 = 1,3
!           fnum = 200+idir
!           CALL copy_sm_to_fm(op_mom_der_ao(idir,idir2)%matrix,fm_work2,error=error)
!           write(fnum,*) ' op_mom_der_ao', idir, idir2
!           DO  i=1,nao
!             write(fnum,'(i5,50f10.6)') i, (fm_work2%local_data(i,j),j=1,nao)
!           END DO
!         END DO
!       END DO
!       stop 'check operators'
!dbg
!     EXTRACT from the ground state full MOS matrix the coefficients of istate 
      CALL cp_fm_get_submatrix(mo_coeff,&
           vecbuf_0,1,istate,nao,1,transpose=.TRUE.,error=error)
      vecbuf_0(1,1:nao) = vecbuf_0(1,1:nao) * occupation(istate)
      CALL cp_fm_set_submatrix(psi0_istate,vecbuf_0,1,1,&
           nao,1,transpose=.TRUE.,error=error)


!     Extract from the full MOS matrices of the response functions the coefficients
!     referring to istate
!     psi1_rxp_iB_istate
      CALL cp_fm_get_submatrix(psi1_rxp(ispin,iB)%matrix,&
           vecbuf_1_rxp,1,istate,nao,1,transpose=.TRUE.,&
           error=error)

      IF(epr_env%full_epr) THEN
!     psi1_d_iB_istate
        CALL cp_fm_get_submatrix(psi1_D(ispin,iB)%matrix,&
             vecbuf_1_d,1,istate,nao,1,transpose=.TRUE.,&
             error=error)
        vecbuf_1_rxp(1,1:nao) = vecbuf_1_rxp(1,1:nao) - vecbuf_1_d(1,1:nao)
      END IF
      CALL cp_fm_set_submatrix(psi_rxp_istate,vecbuf_1_rxp,1,1,&
           nao,1,transpose=.TRUE.,error=error)


!     psi1_p_iiB_istate and psi1_p_iiiB_istate
      CALL cp_fm_get_submatrix(psi1_p(ispin,iiB)%matrix,&
           vecbuf_1_p(iiB)%array,1,istate,nao,1,transpose=.TRUE.,&
           error=error)
      CALL cp_fm_set_submatrix(psi_p1_istate,vecbuf_1_p(iiB)%array,1,1,&
           nao,1,transpose=.TRUE.,error=error)
      CALL cp_fm_get_submatrix(psi1_p(ispin,iiiB)%matrix,&
           vecbuf_1_p(iiiB)%array,1,istate,nao,1,transpose=.TRUE.,&
           error=error)
      CALL cp_fm_set_submatrix(psi_p2_istate,vecbuf_1_p(iiiB)%array,1,1,&
           nao,1,transpose=.TRUE.,error=error)


!     Multuply left and right by the appropriate coefficients and sum into the
!     correct component of the chi tensor using the appropriate multiplicative factor
!     (don't forget the occupation number)
!     Loop over the cartesian components of the tensor 
!     The loop over the components of the external field is external, thereby
!     only one column of the chi tensor is computed here 
      DO idir = 1,3
!dbg
!   write(*,*) 'istate = ', istate, ' idir =', idir 
  !     get ii and iii
        CALL set_vecp(idir,ii,iii)

!       term: 2[C0| r_ii |d_iii(C1(rxp-D))]-2[C0| r_iii |d_ii(C1(rxp-D))]
!       the factor 2 should be already included in the matrix elements
        contrib = 0.0_dp
        CALL cp_sm_fm_multiply(op_mom_der_ao(ii,iii)%matrix,psi_rxp_istate,&
             momxpsi_istate,ncol=1,alpha=1.d0,error=error)
        CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
        chi(idir) = chi(idir) + contrib
!dbg
!  write(*,'(2(A,I5),2f12.6)') ' first term a  r',ii,' der',iii, contrib, chi(idir)
!dbg
        contrib = 0.0_dp
        CALL cp_sm_fm_multiply(op_mom_der_ao(iii,ii)%matrix,psi_rxp_istate,&
             momxpsi_istate,ncol=1,alpha=1.d0,error=error)
        CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
        chi(idir) = chi(idir) - contrib
!dbg
!  write(*,'(2(A,I5),2f12.6)') ' first term b  r',iii,' der',ii,-contrib, chi(idir)
!dbg

!       term: -2[C0| r_ii r_iiB | d_iii(C1(piiiB))] + 2[C0| r_iii r_iiB | d_ii(C1(piiiB))] 
!       the factor 2 should be already included in the matrix elements
        contrib = 0.0_dp
        idir2 = ind_m2(ii,iiB)
        CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,iii)%matrix,psi_p2_istate,&
             momxpsi_istate,ncol=1,alpha=1.d0,error=error)
        CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
        chi(idir) = chi(idir) - contrib
        contrib2 = 0.0_dp
        IF(iiB==iii) THEN
          CALL cp_sm_fm_multiply(op_mom_ao(ii)%matrix,psi_p2_istate,&
               momxpsi_istate,ncol=1,alpha=1.d0,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
          chi(idir) = chi(idir) - contrib2
        END IF
!dbg
!  write(*,'(4(A,I5),3f12.6)') ' second term a  idir2 ', idir2, ' r',ii,' r',iiB,' der',iii,-contrib,-contrib2, chi(idir)
!dbg

        contrib = 0.0_dp
        idir2 = ind_m2(iii,iiB)
        CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,ii)%matrix,psi_p2_istate,&
             momxpsi_istate,ncol=1,alpha=1.d0,error=error)
        CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
        chi(idir) = chi(idir) + contrib
        contrib2 = 0.0_dp
        IF(iiB==ii) THEN
          CALL cp_sm_fm_multiply(op_mom_ao(iii)%matrix,psi_p2_istate,&
               momxpsi_istate,ncol=1,alpha=1.d0,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
          chi(idir) = chi(idir) + contrib2
        END IF
!dbg
!  write(*,'(4(A,I5),3f12.6)') ' second term b  idir2 ', idir2, ' r',iii,' r',iiB,' der',ii,contrib, contrib2, chi(idir)
!dbg

!       term: 2[C0| r_ii   r_iiiB | d_iii(C1(piiB))] - 2[C0| r_iii r_iiiB | d_ii(C1(piiB))] 
!       the factor 2 should be already included in the matrix elements
        contrib = 0.0_dp
        idir2 = ind_m2(ii,iiiB)
        CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,iii)%matrix,psi_p1_istate,&
             momxpsi_istate,ncol=1,alpha=1.d0,error=error)
        CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
        chi(idir) = chi(idir) + contrib
        contrib2 = 0.0_dp
        IF(iiB==iii) THEN
          CALL cp_sm_fm_multiply(op_mom_ao(ii)%matrix,psi_p1_istate,&
               momxpsi_istate,ncol=1,alpha=1.d0,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
          chi(idir) = chi(idir) + contrib2
        END IF
!dbg
!  write(*,'(4(A,I5),2f12.6)') ' third term a  idir2 ', idir2, ' r',ii,' r',iiiB,' der',iii,contrib, contrib2, chi(idir)
!dbg
        idir2 = ind_m2(iii,iiiB)
        CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,ii)%matrix,psi_p1_istate,&
             momxpsi_istate,ncol=1,alpha=1.d0,error=error)
        CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
        chi(idir) = chi(idir) - contrib
        contrib2 = 0.0_dp
        IF(iiB==ii) THEN
          CALL cp_sm_fm_multiply(op_mom_ao(iii)%matrix,psi_p1_istate,&
               momxpsi_istate,ncol=1,alpha=1.d0,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib2,error=error)
          chi(idir) = chi(idir) - contrib2
        END IF

!dbg
!  write(*,'(4(A,I5),2f12.6)') ' third term b  idir2 ', idir2, ' r',iii,' r',iiiB,' der',ii,-contrib, chi(idir)
!dbg
!dbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbg
   IF(.FALSE.) THEN
        IF(iii==iiB) THEN
!         term: -[C0 | r_ii | C1(piiiB)]   
          CALL cp_sm_fm_multiply(op_mom_ao(ii)%matrix,psi_p2_istate,&
               momxpsi_istate,ncol=1,alpha=1.d0,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
          chi(idir) = chi(idir) - contrib
!dbg
!  write(*,'(2(A,I5),2f12.6)') ' fourth term iii=iiB,  r',ii,' C',iiiB,-contrib, chi(idir)
!dbg
        ELSEIF(iii==iiiB) THEN
!         term: [C0| r_ii | C1(piiB)]
          CALL cp_sm_fm_multiply(op_mom_ao(ii)%matrix,psi_p1_istate,&
               momxpsi_istate,ncol=1,alpha=1.d0,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
          chi(idir) = chi(idir) + contrib
!dbg
!  write(*,'(2(A,I5),2f12.6)') ' fourth term iii=iiiB,  r',ii,' C',iiB,contrib,chi(idir) 
!dbg
        END IF
        IF(ii==iiB) THEN
!         term: [C0| r_iii | C1(piiiB)]
          CALL cp_sm_fm_multiply(op_mom_ao(iii)%matrix,psi_p2_istate,&
               momxpsi_istate,ncol=1,alpha=1.d0,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
          chi(idir) = chi(idir) + contrib
!dbg
!  write(*,'(2(A,I5),2f12.6)') ' fourth term ii=iiB,  r',iii,' C',iiiB,contrib, chi(idir)
!dbg
        ELSEIF(ii==iiiB) THEN
!         term: -[C0| r_iii |C1(piiB)]
          CALL cp_sm_fm_multiply(op_mom_ao(iii)%matrix,psi_p1_istate,&
               momxpsi_istate,ncol=1,alpha=1.d0,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib,error=error)
          chi(idir) = chi(idir) - contrib
!dbg
!  write(*,'(2(A,I5),2f12.6)') ' fourth term ii=iiiB, r',iii,' C',iiB,-contrib, chi(idir)
!dbg
        END IF
   END IF
!dbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbgdbg
      END DO  ! idir


    END DO  ! istate
!R  WRITE(*,*) "CHI TENSOR : iB =", iB 
    DO idir = 1,3
      epr_env%bind0(idir,iB) = chi(idir)
!R    epr_env%chi_tensor(idir,iB) =  epr_env%chi_tensor(idir,iB)+&
!R                                   chi(idir)*epr_env%chi_factor/2.0_dp
!R    WRITE(*,*)  idir, epr_env%chi_tensor(idir,iB)
    END DO 


!   deallocate the sparse matrices
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
    END DO 

    CALL timestop(0.0_dp,handle)
    
  END SUBROUTINE bind0_soft_analytic


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
!R             i_B,idir,epr_section,error)
!R
!R    TYPE(epr_env_type)                       :: epr_env
!R    TYPE(pw_env_type), POINTER               :: pw_env
!R    TYPE(particle_type), DIMENSION(:), &
!R      POINTER                                :: particle_set
!R    TYPE(cell_type), POINTER                 :: cell
!R    TYPE(coeff_type)                         :: shift_pw_rspace
!R    INTEGER, INTENT(IN)                      :: i_B, idir
!R    TYPE(section_vals_type), POINTER         :: epr_section
!R    TYPE(cp_error_type), INTENT(INOUT), &
!R      OPTIONAL                               :: error
!R
!R    CHARACTER(LEN=*), PARAMETER :: routineN = 'interpolate_shift_pwgrid', &
!R      routineP = moduleN//':'//routineN
!R
!R    INTEGER                                  :: aint_precond, iat, iatom, &
!R                                                max_iter, n_nics, natom, &
!R                                                precond_kind
!R    INTEGER, DIMENSION(:), POINTER           :: cs_atom_list
!R    LOGICAL                                  :: failure, success
!R    REAL(dp)                                 :: eps_r, eps_x, R_iatom(3), &
!R                                               ra(3), shift_val
!R    REAL(dp), DIMENSION(:, :), POINTER       :: r_nics
!R    REAL(dp), DIMENSION(:, :, :), POINTER    :: chemical_shift, &
!R                                                chemical_shift_nics
!R    TYPE(coeff_type)                         :: shiftspl
!R    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
!R    TYPE(pw_spline_precond_type), POINTER    :: precond
!R    TYPE(section_vals_type), POINTER         :: interp_section
!R
!R!dbg
!R    INTEGER :: bo(2,3), i, k, j, ig
!R    REAL(dp) :: dr(3), rgrid(3), val_grid, val_inter
!R!dbg
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
!R!dbg
!R!  write(*,*) precond_kind, max_iter, eps_r, eps_x
!R!dbg
!R
!R      ! calculate spline coefficients
!R      CALL pw_env_get(pw_env,  auxbas_pw_pool=auxbas_pw_pool)
!R      CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=shiftspl, &
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
!R           n_nics=n_nics,r_nics=r_nics)
!R
!R      IF(ASSOCIATED(cs_atom_list)) THEN
!R        natom = SIZE(cs_atom_list,1)
!R      ELSE
!R        natom = -1
!R      ENDIF
!R  !dbg
!R    WRITE(*,*) "natom ", natom, " idir ", idir , " i_B " , i_B
!R  !dbg
!R
!R      DO iat = 1,natom
!R        iatom = cs_atom_list(iat)
!R        R_iatom = pbc(particle_set(iatom)%r,cell)
!R        shift_val = Eval_Interp_Spl3_pbc(R_iatom,shiftspl%pw,error=error)
!R        chemical_shift(idir,i_B,iatom)= chemical_shift(idir,i_B,iatom)+&
!R                                        epr_env%shift_factor*shift_val
!R  !dbg
!R    WRITE(*,*) " CS iatom ", iatom, chemical_shift(idir,i_B,iatom)
!R  !dbg
!R      END DO 
!R!dbg check interpolation
!R!      bo = shift_pw_rspace%pw%pw_grid%bounds_local
!R!      dr = shift_pw_rspace%pw%pw_grid%dr
!R!      ig = 0
!R!      DO k =  bo(1,3),bo(2,3)
!R!        rgrid(3) = REAL(k-bo(1,3))*dr(3)
!R!        DO j = bo(1,2), bo(2,2)
!R!          rgrid(2) = REAL(j-bo(1,2))*dr(2)
!R!          DO i =  bo(1,1), bo(2,1) 
!R!            rgrid(1) = REAL(i-bo(1,1))*dr(1)
!R!            val_grid = shift_pw_rspace%pw%cr3d(i,j,k)
!R!            val_inter = Eval_Interp_Spl3_pbc(rgrid,shiftspl%pw,error=error)
!R!            ig = ig + 1
!R!            write(200,'(I10,3f12.6,2f16.8,ES12.4)') ig , rgrid(1:3), &
!R!                  val_grid, val_inter, (val_grid-val_inter) 
!R!          END DO
!R!        END DO 
!R!      END DO 
!R!dbg
!R
!R      IF(epr_env%do_nics) THEN
!R
!R        DO iatom = 1,n_nics
!R          ra(1:3) = r_nics(1:3,iatom)
!R          R_iatom = pbc(ra,cell)
!R          shift_val = Eval_Interp_Spl3_pbc(R_iatom,shiftspl%pw,error)
!R          chemical_shift_nics(idir,i_B,iatom)= chemical_shift_nics(idir,i_B,iatom)+&
!R                                               epr_env%shift_factor*shift_val
!R        END DO
!R      END IF
!R
!R      CALL pw_pool_give_back_coeff(auxbas_pw_pool,shiftspl,error=error)
!R
!R    END IF
!R
!R  END SUBROUTINE interpolate_shift_pwgrid

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
    TYPE(coeff_type), POINTER                :: rho_gspace
    TYPE(coeff_type)                         :: funcG_times_rho
    INTEGER, INTENT(IN)                      :: idir
    REAL(dp), INTENT(IN)                     :: my_chi
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
               use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, error=error )

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
!R
!R    TYPE(epr_env_type)                       :: epr_env
!R    TYPE(qs_environment_type), POINTER       :: qs_env
!R    TYPE(section_vals_type), POINTER         :: epr_section
!R    TYPE(cp_error_type), INTENT(INOUT), &
!R      OPTIONAL                               :: error
!R
!R    CHARACTER(len=*), PARAMETER :: routineN = 'epr_print_shift', &
!R      routineP = moduleN//':'//routineN
!R
!R    CHARACTER(LEN=2)                         :: element_symbol
!R    CHARACTER(LEN=default_string_length)     :: name, title
!R    INTEGER                                  :: i, iatom, ir, istat, j, &
!R                                                n_nics, nat_print, natom, &
!R                                                unit_atoms, unit_nics
!R    INTEGER, DIMENSION(:), POINTER           :: cs_atom_list
!R    LOGICAL                                  :: failure, gapw
!R    REAL(dp) :: chi_aniso, chi_iso, chi_sym(3,3), chi_sym_tot(3,3), &
!R      chi_tensor(3,3), chi_tensor_loc(3,3), chi_tmp(3,3), eigenv(3), &
!R      eigenv_tmp(3), rpos(3), shift_nosym(3,3), shift_sym(3,3), shift_tmp(3,3)
!R    REAL(dp), DIMENSION(:, :), POINTER       :: eigenv_chi, r_nics
!R    REAL(dp), DIMENSION(:, :, :), POINTER    :: chemical_shift, &
!R                                                chemical_shift_loc, &
!R                                                chemical_shift_loc_nics, &
!R                                                chemical_shift_nics
!R    TYPE(atomic_kind_type), POINTER          :: atom_kind
!R    TYPE(cp_logger_type), POINTER            :: logger
!R    TYPE(dft_control_type), POINTER          :: dft_control
!R    TYPE(particle_type), DIMENSION(:), &
!R      POINTER                                :: particle_set
!R
!R    NULLIFY (chemical_shift, chemical_shift_nics, r_nics)
!R    NULLIFY (chemical_shift_loc, chemical_shift_loc_nics)
!R    NULLIFY (logger, particle_set,atom_kind,dft_control)
!R    NULLIFY (eigenv_chi)
!R!    ionode = logger%para_env%mepos==logger%para_env%source
!R    logger => cp_error_get_logger(error)
!R
!R    CALL get_epr_env(epr_env, chemical_shift=chemical_shift, &
!R         chemical_shift_nics=chemical_shift_nics, &
!R         chemical_shift_loc=chemical_shift_loc, &
!R         chemical_shift_loc_nics=chemical_shift_loc_nics, &
!R         chi_tensor=chi_tensor, chi_tensor_loc=chi_tensor_loc,&
!R         cs_atom_list=cs_atom_list, n_nics=n_nics, r_nics=r_nics)
!R
!R    CALL get_qs_env(qs_env, dft_control=dft_control, particle_set=particle_set)
!R    natom = SIZE(particle_set,1)
!R    gapw = dft_control%qs_control%gapw
!R
!R    ALLOCATE (eigenv_chi(3,natom), STAT=istat)
!R
!R
!R    ! Finalize Chi calculation
!R    IF(gapw) THEN
!R      ! Symmetrize
!R      DO i = 1,3
!R        DO j = 1,3
!R          chi_sym(i,j) = (chi_tensor(i,j)+chi_tensor(j,i))/2.0_dp
!R          chi_sym_tot(i,j) = (chi_tensor(i,j)+chi_tensor(j,i)+&
!R                          chi_tensor_loc(i,j)+chi_tensor_loc(j,i))/2.0_dp
!R        END DO 
!R      END DO
!R
!R    ELSE
!R
!R      ! Symmetrize
!R      DO i = 1,3
!R        DO j = 1,3
!R          chi_sym(i,j) = (chi_tensor(i,j)+chi_tensor(j,i))/2.0_dp
!R        END DO 
!R      END DO
!R
!R      CALL  dcopy(9,chi_sym,1,chi_sym_tot,1)
!R
!R    END IF
!R
!R    CALL  dcopy(9,chi_sym_tot,1,chi_tmp,1)
!R    CALL  diamat_all(chi_tmp,eigenv)
!R
!R    chi_iso   = (eigenv(1)+ eigenv(2)+eigenv(3))/3.d0
!R    chi_aniso =  eigenv(3)-(eigenv(2)+eigenv(1))/2.d0
!R
!R    IF (BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
!R              "PRINT%CHI_TENSOR",error=error),cp_p_file)) THEN
!R
!R      unit_atoms=cp_print_key_unit_nr(logger,epr_section,"PRINT%CHI_TENSOR",&
!R               extension=".data",middle_name="",log_filename=.FALSE.,&
!R               error=error)
!R
!R      WRITE(title,'(A)') "Magnetic Susceptibility Tensor "
!R      IF(unit_atoms > 0) THEN
!R
!R        WRITE(unit_atoms,'(T2,A)') title
!R        IF(gapw) THEN
!R
!R           WRITE(unit_atoms,'(T30,A)') " CHI from SOFT J in 10^-30 J/T^2 units"
!R           WRITE(unit_atoms,'(3(A,f15.6))') ' XX=',chi_tensor(1,1),&
!R                 ' XY=',chi_tensor(1,2),' XZ=',chi_tensor(1,3)
!R           WRITE(unit_atoms,'(3(A,f15.6))') ' YX=',chi_tensor(2,1),&
!R                 ' YY=',chi_tensor(2,2),' YZ=',chi_tensor(2,3)
!R           WRITE(unit_atoms,'(3(A,f15.6))') ' ZX=',chi_tensor(3,1),&
!R                 ' ZY=',chi_tensor(3,2),' ZZ=',chi_tensor(3,3) 
!R           WRITE(unit_atoms,'(T30,A)') " CHI from LOCAL J in 10^-30 J/T^2 units"
!R           WRITE(unit_atoms,'(3(A,f15.6))') ' XX=',chi_tensor_loc(1,1),&
!R                 ' XY=',chi_tensor_loc(1,2),' XZ=',chi_tensor_loc(1,3)
!R           WRITE(unit_atoms,'(3(A,f15.6))') ' YX=',chi_tensor_loc(2,1),&
!R                 ' YY=',chi_tensor_loc(2,2),' YZ=',chi_tensor_loc(2,3)
!R           WRITE(unit_atoms,'(3(A,f15.6))') ' ZX=',chi_tensor_loc(3,1),&
!R                 ' ZY=',chi_tensor_loc(3,2),' ZZ=',chi_tensor_loc(3,3) 
!R
!R        ELSE
!R
!R           WRITE(unit_atoms,'(T30,A)') " CHI from J in 10^-30 J/T^2 units"
!R           WRITE(unit_atoms,'(3(A,f15.6))')' XX=',chi_tensor(1,1),&
!R                 ' XY=',chi_tensor(1,2),' XZ=',chi_tensor(1,3)
!R           WRITE(unit_atoms,'(3(A,f15.6))')' YX=',chi_tensor(2,1),&
!R                 ' YY=',chi_tensor(2,2),' YZ=',chi_tensor(2,3)
!R           WRITE(unit_atoms,'(3(A,f15.6))')' ZX=',chi_tensor(3,1),&
!R                 ' ZY=',chi_tensor(3,2),' ZZ=',chi_tensor(3,3)
!R
!R        END IF
!R
!R        write (unit_atoms,'(/T2,3(A,F15.5))')' PV1=',eigenv(1),&
!R                 ' PV2=',eigenv(2),' PV3=',eigenv(3)
!R        write (unit_atoms,'(T2,A,F15.5,10X,A,F15.5)')&
!R                 '  iso=',chi_iso,  'aniso=',chi_aniso
!R
!R        CALL dscal(3*3,epr_env%chi_SI2ppmcgs,chi_sym_tot,1)
!R        WRITE(unit_atoms,'(T30,A)') " Total CHI in ppm-cgs units"
!R        WRITE(unit_atoms,'(3(A,f15.6))')' XX=',chi_sym_tot(1,1),&
!R                 ' XY=',chi_sym_tot(1,2),' XZ=',chi_sym_tot(1,3)
!R        WRITE(unit_atoms,'(3(A,f15.6))')' YX=',chi_sym_tot(2,1),&
!R                 ' YY=',chi_sym_tot(2,2),' YZ=',chi_sym_tot(2,3)
!R        WRITE(unit_atoms,'(3(A,f15.6))')' ZX=',chi_sym_tot(3,1),&
!R                 ' ZY=',chi_sym_tot(3,2),' ZZ=',chi_sym_tot(3,3)
!R        WRITE (unit_atoms, '(/T2,3(A,F15.5))') &
!R                 ' PV1=',epr_env%chi_SI2ppmcgs*eigenv(1),&
!R                 ' PV2=',epr_env%chi_SI2ppmcgs*eigenv(2),&
!R                 ' PV3=',epr_env%chi_SI2ppmcgs*eigenv(3)
!R        WRITE (unit_atoms,'(T2,A,F15.5,10X,A,F15.5)')&
!R                 '  iso=',epr_env%chi_SI2ppmcgs*chi_iso, &
!R                 'aniso=',epr_env%chi_SI2ppmcgs*chi_aniso
!R      END IF
!R
!R      CALL cp_print_key_finished_output(unit_atoms, logger,epr_section,&
!R               "PRINT%CHI_TENSOR", error=error)
!R
!R    END IF  ! print chi
!R
!R
!R    !Finalize the Shift calculation
!R
!R    IF (BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
!R              "PRINT%SHIELDING_TENSOR",error=error),cp_p_file)) THEN
!R
!R      unit_atoms=cp_print_key_unit_nr(logger,epr_section,"PRINT%SHIELDING_TENSOR",&
!R               extension=".data",middle_name="CS",log_filename=.FALSE.,&
!R               error=error)
!R
!R      nat_print = SIZE(cs_atom_list,1)
!R      WRITE(title,'(A,1X,I5)') "Shielding atom at atomic positions. # tensors printed ", nat_print
!R      IF(unit_atoms > 0) THEN
!R
!R        WRITE(unit_atoms,'(T2,A)') title
!R        DO ir = 1,nat_print
!R          iatom = cs_atom_list(ir)
!R          rpos(1:3) = particle_set(iatom)%r(1:3)
!R          atom_kind => particle_set(iatom)%atomic_kind
!R          CALL get_atomic_kind(atom_kind,name=name,element_symbol=element_symbol)
!R
!R          DO i = 1,3
!R            DO j = 1,3
!R              shift_nosym(i,j) = chi_tensor(i,j)+chemical_shift(i,j,iatom)
!R              shift_sym(i,j) = chi_sym(i,j)+(chemical_shift(i,j,iatom)+&
!R                               chemical_shift(j,i,iatom))/2.0_dp
!R
!R              IF(gapw) THEN
!R                shift_nosym(i,j) = shift_nosym(i,j)+chemical_shift_loc(i,j,iatom)
!R                shift_sym(i,j) = shift_sym(i,j)+(chemical_shift_loc(i,j,iatom)+&
!R                                 chemical_shift_loc(j,i,iatom))/2.0_dp
!R              END IF
!R            END DO 
!R          END DO 
!R
!R          CALL  dcopy(9,shift_sym,1,shift_tmp,1)
!R          CALL  diamat_all(shift_tmp,eigenv_tmp)
!R          eigenv_chi(1:3, iatom) = eigenv_tmp(1:3)
!R          
!R          WRITE(unit_atoms,'(T2,I5,A,2X,A2,2X,3f15.6)') iatom,name,element_symbol,rpos(1:3)
!R          WRITE(unit_atoms,'(A,3f15.6,5x,f15.6)') 'S_x = ',shift_nosym(1, 1:3), eigenv_tmp(1)
!R          WRITE(unit_atoms,'(A,3f15.6,5x,f15.6)') 'S_y = ',shift_nosym(2, 1:3), eigenv_tmp(2)
!R          WRITE(unit_atoms,'(A,3f15.6,5x,f15.6)') 'S_z = ',shift_nosym(3, 1:3), eigenv_tmp(3)
!R        END DO  ! ir
!R
!R      END IF
!R
!R      CALL cp_print_key_finished_output(unit_atoms, logger,epr_section,&
!R               "PRINT%SHIELDING_TENSOR", error=error)
!R
!R      IF(epr_env%do_nics) THEN
!R        unit_nics=cp_print_key_unit_nr(logger,epr_section,"PRINT%SHIELDING_TENSOR",&
!R                 extension=".data",middle_name="NICS",log_filename=.FALSE.,&
!R                 error=error)
!R
!R        WRITE(title,'(A,1X,I5)') "Shielding atom at atomic positions. # tensors printed ", n_nics
!R        IF(unit_nics > 0) THEN
!R          WRITE(unit_nics,'(T2,A)') title
!R          DO ir = 1,nat_print
!R            rpos(1:3) = r_nics(1:3,iatom)
!R          
!R            WRITE(unit_nics,'(T2,I5,2X,3f15.6)') iatom,rpos(1:3)
!R            IF(gapw) THEN
!R              WRITE(unit_nics,'(3f15.6,5x,3f15.6)') chemical_shift_nics(1, 1:3, iatom),& 
!R                    chemical_shift_loc_nics(1,1:3, iatom)
!R              WRITE(unit_nics,'(3f15.6,5x,3f15.6)') chemical_shift_nics(2, 1:3, iatom),& 
!R                    chemical_shift_loc_nics(2,1:3, iatom)
!R              WRITE(unit_nics,'(3f15.6,5x,3f15.6)') chemical_shift_nics(3, 1:3, iatom),& 
!R                    chemical_shift_loc_nics(3,1:3, iatom)
!R            ELSE
!R              WRITE(unit_nics,'(3f15.6)') chemical_shift_nics(1, 1:3, iatom)
!R              WRITE(unit_nics,'(3f15.6)') chemical_shift_nics(2, 1:3, iatom)
!R              WRITE(unit_nics,'(3f15.6)') chemical_shift_nics(3, 1:3, iatom)
!R            END  IF
!R          END DO
!R        END IF  
!R
!R        CALL cp_print_key_finished_output(unit_nics, logger,epr_section,&
!R                 "PRINT%SHIELDING_TENSOR", error=error)
!R      END IF
!R
!R    END IF  ! print shift
!R
!R
!R
!R  END SUBROUTINE epr_print_shift

END MODULE qs_linres_epr_bind0



