
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****m* cp2k/qs_linres_nmr_shift *
!!
!!   NAME 
!!     qs_linres_nmr_current
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

MODULE qs_linres_nmr_shift

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
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply,copy_sm_to_fm 
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_get_rval,&
                                             section_get_lval,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_coeff,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_spline_utils,                 ONLY: Eval_Interp_Spl3_pbc,&
                                             find_coeffs,&
                                             pw_spline_precond_type,&
                                             pw_spline_precond_create,&
                                             pw_spline_precond_release,&
                                             pw_spline_do_precond,&
                                             pw_spline_precond_set_kind,&
                                             spl3_pbc
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             RECIPROCALSPACE,&
                                             REALSPACE,&
                                             pw_transfer,&
                                             pw_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_nmr_op,                ONLY: fac_vecp,&
                                             ind_m2,&
                                             set_vecp,&
                                             set_vecp_rev
  USE qs_linres_types,                 ONLY: linres_control_type,&
                                             nmr_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_operators_ao,                 ONLY: rRc_xyz_der_ao,set_up_op_sm
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
  PUBLIC :: chi_soft_analytic, interpolate_shift_pwgrid, mult_G_ov_G2_grid

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_nmr_shift'

!!***
! **************

CONTAINS

  SUBROUTINE chi_soft_analytic(nmr_env,qs_env,ispin,iB,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    INTEGER, INTENT(IN)                      :: ispin
    INTEGER, INTENT(IN)                      :: iB
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'chi_soft_analytic', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, idir, idir2, ii, iiB, iii, iiiB, istat, istate, nao
              
    LOGICAL :: failure, gapw
    REAL(dp)      :: chi(3), contrib, dk(3)
    REAL(dp), DIMENSION(:), POINTER          :: occupation
    REAL(dp), DIMENSION(:,:), POINTER        :: vecbuf_0
    REAL(dp), DIMENSION(:,:), POINTER        :: vecbuf_1_rxp
    REAL(dp), DIMENSION(:,:), POINTER        :: vecbuf_1_d
    TYPE(cp_2d_r_p_type)                     :: vecbuf_1_p(3)
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_fm_type), POINTER                :: psi0_istate
    TYPE(cp_fm_type), POINTER                :: psi_rxp_istate
    TYPE(cp_fm_type), POINTER                :: psi_p1_istate
    TYPE(cp_fm_type), POINTER                :: psi_p2_istate
    TYPE(cp_fm_type), POINTER                :: momxpsi_istate
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: op_mom_ao
    TYPE(real_matrix_p_type), DIMENSION(:,:), &
      POINTER                                :: op_mom_der_ao
!dbg
     INTEGER :: fnum, i, j
     TYPE(cp_fm_type), POINTER :: fm_work2
!dbg
 

    CALL timeset(routineN,"I"," ",handle)

    failure = .FALSE.

    NULLIFY(dft_control,matrix_s,mos,para_env)
    NULLIFY(mo_coeff,occupation)
    NULLIFY(op_mom_ao,op_mom_der_ao)
    NULLIFY(psi0_istate,psi_rxp_istate,psi_p1_istate,psi_p2_istate,momxpsi_istate)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,matrix_s=matrix_s,&
         mos=mos,&
         para_env=para_env)
    gapw = dft_control%qs_control%gapw

!   Allocate sparse matric for dipole, quadrupole and their derivatives => 9x3
!   Remember the derivatives are antisymmetric
    CALL allocate_matrix_set(op_mom_ao,9)
    CALL allocate_matrix_set(op_mom_der_ao,9,3)
    CALL set_up_op_sm(op_mom_ao(1)%matrix,qs_env,symmetry="none",&
          name="op_mom",error=error)    
    DO idir2=1,3
       CALL replicate_matrix_structure(op_mom_ao(1)%matrix, &
            op_mom_der_ao(1,idir2)%matrix,&
            "op_mom_der_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir2))),&
            target_symmetry="none")
       CALL set_matrix(op_mom_der_ao(1,idir2)%matrix,0.0_dp)
    END DO

    DO idir = 2,SIZE(op_mom_ao,1)
      CALL replicate_matrix_structure(op_mom_ao(1)%matrix, &
           op_mom_ao(idir)%matrix,&
           "op_mom_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir))),&
           target_symmetry="none")
      CALL set_matrix(op_mom_ao(idir)%matrix,0.0_dp)
      DO idir2=1,3
        CALL replicate_matrix_structure(op_mom_ao(1)%matrix, &
             op_mom_der_ao(idir,idir2)%matrix,&
             "op_mom_der_ao"//"-"//TRIM(ADJUSTL(cp_to_string(idir*idir2))),&
             target_symmetry="none")
        CALL set_matrix(op_mom_der_ao(idir,idir2)%matrix,0.0_dp)
      END DO
    END DO

!   get iiB and iiiB
    CALL set_vecp(iB,iiB,iiiB)

!   get ground state MOS
    CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,&
                    occupation_numbers=occupation, nao=nao)

!   Allocate full matrices for only one vector
    NULLIFY(tmp_fm_struct)
    CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
         ncol_global=1,para_env=para_env,context=mo_coeff%matrix_struct%context)
    CALL cp_fm_create (psi0_istate, tmp_fm_struct )
    CALL cp_fm_create (psi_rxp_istate, tmp_fm_struct )
    CALL cp_fm_create (psi_p1_istate, tmp_fm_struct )
    CALL cp_fm_create (psi_p2_istate, tmp_fm_struct )
    CALL cp_fm_create (momxpsi_istate, tmp_fm_struct )
    CALL cp_fm_struct_release ( tmp_fm_struct )
!dbg
    CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
          ncol_global=nao,para_env=para_env,context=mo_coeff%matrix_struct%context)
    CALL cp_fm_create (fm_work2, tmp_fm_struct )
    CALL cp_fm_struct_release ( tmp_fm_struct )
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
    DO istate = 1,nmr_env%nstates(ispin)

!     Get the Wannier center of the istate-th ground state orbital
      dk(1:3) = nmr_env%centers_set(ispin)%array(1:3,istate)

!     Compute the multipole integrals for the state istate,
!     using as reference center the corresponding Wannier center
      DO idir = 1,9
         CALL set_matrix(op_mom_ao(idir)%matrix, 0.0_dp)
         DO idir2 = 1,3
           CALL set_matrix(op_mom_der_ao(idir,idir2)%matrix, 0.0_dp)
         END DO
      END DO
      CALL rRc_xyz_der_ao(op_mom_ao,op_mom_der_ao,qs_env,dk,order=2,&
           minimum_image=.TRUE.,soft=gapw,error=error)

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
      CALL cp_fm_get_submatrix(nmr_env%psi1_rxp(ispin,iB)%matrix,&
           vecbuf_1_rxp,1,istate,nao,1,transpose=.TRUE.,&
           error=error)

      IF(nmr_env%full_nmr) THEN
!     psi1_d_iB_istate
        CALL cp_fm_get_submatrix(nmr_env%psi1_D(ispin,iB)%matrix,&
             vecbuf_1_d,1,istate,nao,1,transpose=.TRUE.,&
             error=error)
        vecbuf_1_rxp(1,1:nao) = vecbuf_1_rxp(1,1:nao) - vecbuf_1_d(1,1:nao)
      END IF
      CALL cp_fm_set_submatrix(psi_rxp_istate,vecbuf_1_rxp,1,1,&
           nao,1,transpose=.TRUE.,error=error)


!     psi1_p_iiB_istate and psi1_p_iiiB_istate
      CALL cp_fm_get_submatrix(nmr_env%psi1_p(ispin,iiB)%matrix,&
           vecbuf_1_p(iiB)%array,1,istate,nao,1,transpose=.TRUE.,&
           error=error)
      CALL cp_fm_set_submatrix(psi_p1_istate,vecbuf_1_p(iiB)%array,1,1,&
           nao,1,transpose=.TRUE.,error=error)
      CALL cp_fm_get_submatrix(nmr_env%psi1_p(ispin,iiiB)%matrix,&
           vecbuf_1_p(iiiB)%array,1,istate,nao,1,transpose=.TRUE.,&
           error=error)
      CALL cp_fm_set_submatrix(psi_p2_istate,vecbuf_1_p(iiiB)%array,1,1,&
           nao,1,transpose=.TRUE.,error=error)


!     Multuply left and right for the appropriate coefficients and sum into the
!     correct component of the chi tensor using the appropriate multiplicative factor
!     (don't forget the occupation number)
!     Loop over the cartesian components of the tensor 
!     The loop over the components of the external field is external, thereby
!     only one column of the chi tensor is computed here 
      DO idir = 1,3
!dbg
   write(*,*) 'istate = ', istate, ' idir =', idir 
  !     get ii and iii
        CALL set_vecp(idir,ii,iii)

!       term: 2[C0| r_ii |d_iii(C1(rxp-D))]-2[C0| r_iii |d_ii(C1(rxp-D))]
!       the factor 2 should be already included in the matrix elements
        CALL cp_sm_fm_multiply(op_mom_der_ao(ii,iii)%matrix,psi_rxp_istate,&
             momxpsi_istate,ncol=1,alpha=1.d0,error=error)
        CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib)
        chi(idir) = chi(idir) + contrib
!dbg
  write(*,'(2(A,I5),2f12.6)') ' first term a  r',ii,' der',iii, contrib, chi(idir)
!dbg
        CALL cp_sm_fm_multiply(op_mom_der_ao(iii,ii)%matrix,psi_rxp_istate,&
             momxpsi_istate,ncol=1,alpha=1.d0,error=error)
        CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib)
        chi(idir) = chi(idir) - contrib
!dbg
  write(*,'(2(A,I5),2f12.6)') ' first term b  r',iii,' der',ii,-contrib, chi(idir)
!dbg

!       term: -2[C0| r_ii r_iiB | d_iii(C1(piiiB))] + 2[C0| r_iii r_iiB | d_ii(C1(piiiB))] 
!       the factor 2 should be already included in the matrix elements
        idir2 = ind_m2(ii,iiB)
        CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,iii)%matrix,psi_p2_istate,&
             momxpsi_istate,ncol=1,alpha=1.d0,error=error)
        CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib)
        chi(idir) = chi(idir) - contrib
!dbg
  write(*,'(4(A,I5),2f12.6)') ' second term a  idir2 ', idir2, ' r',ii,' r',iiB,' der',iii,-contrib, chi(idir)
!dbg
        idir2 = ind_m2(iii,iiB)
        CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,ii)%matrix,psi_p2_istate,&
             momxpsi_istate,ncol=1,alpha=1.d0,error=error)
        CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib)
        chi(idir) = chi(idir) + contrib
!dbg
  write(*,'(4(A,I5),2f12.6)') ' second term b  idir2 ', idir2, ' r',iii,' r',iiB,' der',ii,contrib, chi(idir)
!dbg

!       term: 2[C0| r_ii   r_iiiB | d_iii(C1(piiB))] - 2[C0| r_iii r_iiiB | d_ii(C1(piiB))] 
!       the factor 2 should be already included in the matrix elements
        idir2 = ind_m2(ii,iiiB)
        CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,iii)%matrix,psi_p1_istate,&
             momxpsi_istate,ncol=1,alpha=1.d0,error=error)
        CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib)
        chi(idir) = chi(idir) + contrib
!dbg
  write(*,'(4(A,I5),2f12.6)') ' third term a  idir2 ', idir2, ' r',ii,' r',iiiB,' der',iii,contrib, chi(idir)
!dbg
        idir2 = ind_m2(iii,iiiB)
        CALL cp_sm_fm_multiply(op_mom_der_ao(idir2,ii)%matrix,psi_p1_istate,&
             momxpsi_istate,ncol=1,alpha=1.d0,error=error)
        CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib)
        chi(idir) = chi(idir) - contrib!
!dbg
  write(*,'(4(A,I5),2f12.6)') ' third term b  idir2 ', idir2, ' r',iii,' r',iiiB,' der',ii,-contrib, chi(idir)
!dbg

        IF(iii==iiB) THEN
!         term: -[C0 | r_ii | C1(piiiB)]   
          CALL cp_sm_fm_multiply(op_mom_ao(ii)%matrix,psi_p2_istate,&
               momxpsi_istate,ncol=1,alpha=1.d0,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib)
          chi(idir) = chi(idir) - contrib
!dbg
  write(*,'(2(A,I5),2f12.6)') ' fourth term iii=iiB,  r',ii,' C',iiiB,-contrib, chi(idir)
!dbg
        ELSEIF(iii==iiiB) THEN
!         term: [C0| r_ii | C1(piiB)]
          CALL cp_sm_fm_multiply(op_mom_ao(ii)%matrix,psi_p1_istate,&
               momxpsi_istate,ncol=1,alpha=1.d0,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib)
          chi(idir) = chi(idir) + contrib
!dbg
  write(*,'(2(A,I5),2f12.6)') ' fourth term iii=iiiB,  r',ii,' C',iiB,contrib,chi(idir) 
!dbg
        END IF
        IF(ii==iiB) THEN
!         term: [C0| r_iii | C1(piiiB)]
          CALL cp_sm_fm_multiply(op_mom_ao(iii)%matrix,psi_p2_istate,&
               momxpsi_istate,ncol=1,alpha=1.d0,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib)
          chi(idir) = chi(idir) + contrib
!dbg
  write(*,'(2(A,I5),2f12.6)') ' fourth term ii=iiB,  r',iii,' C',iiiB,contrib, chi(idir)
!dbg
        ELSEIF(ii==iiiB) THEN
!         term: -[C0| r_iii |C1(piiB)]
          CALL cp_sm_fm_multiply(op_mom_ao(iii)%matrix,psi_p1_istate,&
               momxpsi_istate,ncol=1,alpha=1.d0,error=error)
          CALL cp_fm_trace(psi0_istate,momxpsi_istate,contrib)
          chi(idir) = chi(idir) - contrib
!dbg
  write(*,'(2(A,I5),2f12.6)') ' fourth term ii=iiiB, r',iii,' C',iiB,-contrib, chi(idir)
!dbg
        END IF

      END DO  ! idir


    END DO  ! istate
    write(*,*) "CHI TENSOR : iB =", iB 
    DO idir = 1,3
      nmr_env%chi_tensor(idir,iB) = chi(idir)*nmr_env%chi_factor
      write(*,*)  idir, nmr_env%chi_tensor(idir,iB)
    END DO 

!   deallocate the sparse matrices
    CALL deallocate_matrix_set(op_mom_ao)
    CALL deallocate_matrix_set(op_mom_der_ao)
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
    
  END SUBROUTINE chi_soft_analytic


!  *****************************************************************************************

!!****f* qs_linres_nmr_shift%interpolate_shift_pwgrid
!!
!!  NAME 
!!      nmr_interpolate_shift_pwgrid
!!
!!  FUNCTION
!!      interpolate the shift calculated on the PW grid in order to ger 
!!      the value on arbitrary points in real space
!!
!!  ARGUMENTS
!!      nme_env : to get the shift tensor and the list of additional points 
!!      particle_set : for the atomic position
!!      cell : to take into account the pbs, and to have the volume
!!      shift_pw_rspace : specific component of the shift tensor on the pw grid
!!      i_B : component of the magnetic field for which the shift is calculated (row)
!!      idir : component of the vector \int_{r}[ ((r-r') x j(r))/|r-r'|^3 ] = Bind(r')
!!
!!  NOTES
!!
!!  AUTHOR
!!    MI
!!
!!*** **********************************************************************

  SUBROUTINE interpolate_shift_pwgrid(nmr_env,pw_env,particle_set,cell,shift_pw_rspace,&
             i_B,idir,nmr_section,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type)                         :: shift_pw_rspace
    INTEGER, INTENT(IN)                      :: i_B, idir
    TYPE(section_vals_type), POINTER         :: nmr_section
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'interpolate_shift_pwgrid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: aint_precond, iatom, max_iter,&
                                                natom, precond_kind
    LOGICAL                                  :: failure, success
    REAL(dp)                                 :: eps_r, eps_x, R_iatom(3),&
                                                ra(3), shift_val
    TYPE(coeff_type),  POINTER               :: shiftspl
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_spline_precond_type), POINTER    :: precond
    TYPE(section_vals_type), POINTER         :: interp_section

    failure = .FALSE.
    NULLIFY (interp_section)
    NULLIFY (auxbas_pw_pool,precond,shiftspl)

    CPPrecondition(ASSOCIATED(shift_pw_rspace%pw),cp_failure_level,routineP,error,failure)
    IF(.NOT. failure) THEN

      interp_section =>  section_vals_get_subs_vals(nmr_section,&
                        "INTERPOLATOR",error=error)
      CALL section_vals_val_get(interp_section,"aint_precond", &
                                i_val=aint_precond, error=error)
      CALL section_vals_val_get(interp_section,"precond",i_val=precond_kind, error=error)
      CALL section_vals_val_get(interp_section,"max_iter",i_val=max_iter, error=error)
      CALL section_vals_val_get(interp_section,"eps_r",r_val=eps_r, error=error)
      CALL section_vals_val_get(interp_section,"eps_x",r_val=eps_x, error=error)


      ! calculate spline coefficients
      CALL pw_env_get(pw_env,  auxbas_pw_pool=auxbas_pw_pool)
      CALL pw_pool_init_coeff(pool=auxbas_pw_pool,coeff=shiftspl, &
                             use_data=REALDATA3D,in_space=REALSPACE,&
                             error=error)
        
      CALL pw_spline_precond_create(precond,precond_kind=aint_precond,&
          pool=auxbas_pw_pool,pbc=.TRUE.,transpose=.FALSE.,error=error)
      CALL pw_spline_do_precond(precond,shift_pw_rspace%pw,shiftspl%pw,error=error)
      CALL pw_spline_precond_set_kind(precond,precond_kind,error=error)
      success=find_coeffs(values=shift_pw_rspace%pw,coeffs=shiftspl%pw,&
             linOp=spl3_pbc,preconditioner=precond,pool=auxbas_pw_pool, &
             eps_r=eps_r,eps_x=eps_x,max_iter=max_iter, &
             print_section=interp_section,error=error)
      CPPostconditionNoFail(success,cp_warning_level,routineP,error)
      CALL pw_spline_precond_release(precond,error=error)


      natom = SIZE(particle_set,1)

      DO iatom = 1,natom
        R_iatom = pbc(particle_set(iatom)%r,cell)
        shift_val = Eval_Interp_Spl3_pbc(R_iatom,shiftspl%pw,error)
        nmr_env%chemical_shift(i_B,idir,iatom)= nmr_env%chemical_shift(i_B,idir,iatom)+&
            nmr_env%shift_factor*shift_val
      END DO 

      IF(nmr_env%do_nics) THEN

        DO iatom = 1,nmr_env%n_nics
          ra(1:3) = nmr_env%r_nics(1:3,iatom)
          R_iatom = pbc(ra,cell)
          shift_val = Eval_Interp_Spl3_pbc(R_iatom,shiftspl%pw,error)
          nmr_env%chemical_shift_nics(i_B,idir,iatom) = &
             nmr_env%chemical_shift_nics(i_B,idir,iatom) + &
             nmr_env%shift_factor*shift_val
        END DO
      END IF

      CALL pw_pool_give_back_coeff(auxbas_pw_pool,shiftspl,error=error)

    END IF

  END SUBROUTINE interpolate_shift_pwgrid

!!****f* qs_linres_nmr_shift%mult_G_ov_G2_grid
!!
!!  NAME 
!!      mult_G_ov_G2_grid
!!
!!  FUNCTION
!!      Given the current density on the PW grid in reciprcal space
!!      (obtained by FFT), calculate the integral  
!!        \int_{r}[ ((r-r') x j(r))/|r-r'|^3 ] = Bind(r')
!!      which in reciprcal space reads  (for G/=0)
!!         G/|G|^2 x J(G)    
!!
!!  ARGUMENTS
!!
!!  NOTES
!!     The G=0 component is not comnputed here, but can be evaluated
!!     through the susceptibility and added to the shift in a second time
!!
!!  AUTHOR
!!    MI
!!
!!*** **********************************************************************


  SUBROUTINE  mult_G_ov_G2_grid(cell,pw_pool,rho_gspace,funcG_times_rho,idir,error)


    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_pool_type), POINTER              :: pw_pool
    TYPE(coeff_type), POINTER                :: rho_gspace
    TYPE(coeff_type)                         :: funcG_times_rho
    INTEGER, INTENT(IN)                      :: idir
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
               use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE )

    grid => influence_fn % pw_grid

    DO ig = grid % first_gne0, grid % ngpts_cut_local
       g2 = grid % gsq ( ig )
       influence_fn%cc(ig) = grid % g(idir,ig)/g2
    END DO  ! ig
    IF ( grid % have_g0 )  influence_fn%cc ( 1 ) = 0.0_dp

    frho => funcG_times_rho%pw
    CALL pw_transfer (rho_gspace%pw,frho)
    
    ng = SIZE(grid % gsq)
    frho%cc(1:ng) = frho%cc(1:ng)*influence_fn % cc  ( 1 : ng )

    CALL pw_pool_give_back_pw(pw_pool,influence_fn,&
         accept_non_compatible=.TRUE.,error=error)
 
  END SUBROUTINE  mult_G_ov_G2_grid

END MODULE qs_linres_nmr_shift



