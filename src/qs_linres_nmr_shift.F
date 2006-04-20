
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
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: coset,&
                                             ncoset
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_spline_utils,                 ONLY: Eval_Interp_Spl3_pbc
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             RECIPROCALSPACE,&
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
  USE qs_operators_ao,                 ONLY: rRc_xyz_ao
  USE sparse_matrix_types,             ONLY: set_matrix

#include "cp_common_uses.h"

  IMPLICIT NONE


  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: chi_soft_analytic, interpolate_shift_pwgrid, mult_G_ov_G2_grid

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_nmr_shift'

!!***
! **************

CONTAINS


  SUBROUTINE chi_soft_analytic(nmr_env,qs_env,psi1,p_psi1,ispin,i_B,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: psi1, p_psi1

    INTEGER, INTENT(IN)                      :: ispin
    INTEGER, INTENT(IN)                      :: i_B
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'chi_soft_analytic', &
      routineP = moduleN//':'//routineN


    INTEGER :: homo, i,j,iao, idir, idir2, idir3, ii_B, iii_B, imom, &
               istat, istate, nao, nmo 
    LOGICAL :: failure, gapw, uni_occ
    REAL(dp) :: alpha, contrib,d1dk(3),d1(3), dk(3),maxocc,scale_fac
    REAL(dp), DIMENSION(:), POINTER          :: occupation
    REAL(KIND=dp), DIMENSION(:), POINTER     :: d1mdk_i,d1mdk_ii
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: d1mdk
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_fm_type), POINTER                :: fm_work1,fm_work2
    TYPE(cp_fm_type), POINTER                :: fm_work3
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos


    failure = .FALSE.

    NULLIFY(cell,dft_control,linres_control,mo_coeff,mos,occupation,para_env)
    NULLIFY(fm_work1,fm_work2,fm_work3)

    CALL set_vecp(i_B,ii_B,iii_B)

    CALL get_qs_env(qs_env=qs_env,cell=cell,dft_control=dft_control,&
         linres_control=linres_control,&
         mos=mos,para_env=para_env)
    gapw = (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc)
    CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,&
         occupation_numbers=occupation, homo=homo, nmo=nmo,&
         nao=nao, uniform_occupation=uni_occ,maxocc=maxocc)

    CALL cp_fm_to_fm(nmr_env%psi1_rxp(ispin,i_B)%matrix,psi1(ispin)%matrix)
    IF(nmr_env%full_nmr) THEN
        CALL  cp_fm_scale_and_add(1.0_dp,psi1(ispin)%matrix,&
              -1.0_dp,nmr_env%psi1_D(ispin,i_B)%matrix,error=error)
    END IF

    IF ( .NOT. uni_occ ) THEN
       alpha = 1.0_dp
       CALL cp_fm_column_scale(psi1(ispin)%matrix,occupation(1:homo))
    ELSE
       alpha = maxocc
    END IF

!dbg
! write(*,*) 'size ', nao, nmo, homo
!dbg
!dbg
! write(*,*) 'PSI1 '
! do i = 1,nao
!     write(*,'(I4,4f12.6)') i,(psi1(ispin)%matrix%local_data(i,j),j=1,4)
! enddo
! write(*,*)
!dbg
    ALLOCATE(d1mdk(3,homo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(d1mdk_i(homo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(d1mdk_ii(homo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
  
    d1 = linres_control%localized_wfn_control%centers_set(1)%array(1:3,1)
    DO istate = 1,homo
      dk = linres_control%localized_wfn_control%centers_set(ispin)%array(1:3,istate)
      d1dk = pbc(dk,d1,cell)
      d1mdk(1:3,istate) = d1dk(1:3)
    END DO 

    ! create a new matrix 
    NULLIFY(tmp_fm_struct)
    CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
           ncol_global=homo,para_env=para_env,context=mo_coeff%matrix_struct%context)
    CALL cp_fm_create (fm_work1, tmp_fm_struct )
    CALL cp_fm_set_all(fm_work1,0.0_dp,error=error)
    CALL cp_fm_create (fm_work2, tmp_fm_struct )
    CALL cp_fm_set_all(fm_work2,0.0_dp,error=error)
    CALL cp_fm_create (fm_work3, tmp_fm_struct )
    CALL cp_fm_set_all(fm_work3,0.0_dp,error=error)
    CALL cp_fm_struct_release ( tmp_fm_struct )

!     Lets say that idir gives the j_rho component and compute the relative contributions
    DO idir = 1,3

   ! Apply the component idir of the p operator
       CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,psi1(ispin)%matrix,&
              p_psi1(ispin)%matrix,ncol=homo,alpha=+1.0_dp,error=error)
       ! These additional terms should give zero contribution due to
       ! the orthogonality of the response functions with respect to the unperturbed set
       IF(idir==ii_B) THEN
         CALL  cp_fm_scale_and_add(1.0_dp,p_psi1(ispin)%matrix,-1.0_dp,&
               nmr_env%psi1_p(ispin,iii_B)%matrix,error=error)
       ELSEIF(idir==iii_B)THEN
         CALL  cp_fm_scale_and_add(1.0_dp,p_psi1(ispin)%matrix,+1.0_dp,&
               nmr_env%psi1_p(ispin,ii_B)%matrix,error=error)
       END IF

       IF ( .NOT. uni_occ ) THEN
          CALL cp_fm_column_scale(p_psi1(ispin)%matrix,occupation(1:homo))
       END IF

!dbg
 write(*,*) 'P_PSI1 dir=', idir
! do i = 1,nao
!     write(*,'(I4,4f12.6)') i,(p_psi1(ispin)%matrix%local_data(i,j),j=1,4)
! enddo
! write(*,*)
!dbg

!      idir2 is the component of the position operator
       DO idir2 = 1,3
         IF(idir /= idir2) THEN

!          idir3 is the component of CHI to which this term contribute
             CALL set_vecp_rev(idir,idir2,idir3)
             scale_fac=fac_vecp(idir3,idir2,idir)
             d1mdk_i(1:homo) = d1mdk(idir2,1:homo)

! first term ppsi0 times p1
             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(idir2)%matrix,psi1(ispin)%matrix,&
                  fm_work1,ncol=homo,alpha=1.0_dp,error=error)


             CALL cp_fm_to_fm(psi1(ispin)%matrix,fm_work2,ncol=homo)
             CALL cp_fm_column_scale(fm_work2,d1mdk_i(1:homo))
             CALL  cp_fm_scale_and_add(1.0_dp,fm_work1,1.0_dp,fm_work2,error=error)

             CALL cp_fm_to_fm(nmr_env%p_psi0(ispin,idir)%matrix,&
                 fm_work3,ncol=homo)
             CALL cp_fm_trace(fm_work3,fm_work1,contrib,error=error)
             nmr_env%chi_tensor(i_B,idir3) = nmr_env%chi_tensor(i_B,idir3) +&
                      scale_fac*contrib*alpha

      write(*,*) idir3, idir2, idir, contrib 

! first term psi0 times pp1
             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(idir2)%matrix,p_psi1(ispin)%matrix,&
                  fm_work1,ncol=homo,alpha=1.0_dp,error=error)

             CALL cp_fm_to_fm(p_psi1(ispin)%matrix,fm_work2,ncol=homo)
             CALL cp_fm_column_scale(fm_work2,d1mdk_i(1:homo))
             CALL  cp_fm_scale_and_add(1.0_dp,fm_work1,1.0_dp,fm_work2,error=error)

             CALL cp_fm_to_fm(mo_coeff, fm_work3,ncol=homo)
             CALL cp_fm_trace(fm_work3,fm_work1,contrib,error=error)
             nmr_env%chi_tensor(i_B,idir3) = nmr_env%chi_tensor(i_B,idir3) +&
                      scale_fac*contrib*alpha

      write(*,*) idir3, idir2, idir, contrib 

             imom = ind_m2(idir2,iii_B)
             d1mdk_ii(1:homo) = d1mdk(iii_B,1:homo)
!dbg
!    write(*,*) 1 , 1 , ind_m2(1,1), coset(2,0,0)
!    write(*,*) 2 , 2 , ind_m2(2,2),coset(0,2,0)
!    write(*,*) 3 , 3 , ind_m2(3,3),coset(0,0,2)
!    write(*,*) 1 , 2 , ind_m2(1,2),coset(1,1,0)
!    write(*,*) 1 , 3 , ind_m2(1,3),coset(1,0,1)
!    write(*,*) 2 , 3 , ind_m2(2,3),coset(0,1,1)
!dbg


! second term ppsi0 times p2
             CALL cp_fm_set_all(fm_work1,0.0_dp,error=error)
             CALL cp_fm_set_all(fm_work2,0.0_dp,error=error)
             CALL cp_fm_set_all(fm_work3,0.0_dp,error=error)
             CALL cp_fm_to_fm(nmr_env%psi1_p(ispin,ii_B)%matrix,fm_work2,ncol=homo)
             IF ( .NOT. uni_occ ) THEN
               CALL cp_fm_column_scale(fm_work2,occupation(1:homo))
             END IF
             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(imom)%matrix,&
                  fm_work2,fm_work3,ncol=homo,alpha=1.0_dp,error=error)

             CALL cp_fm_to_fm(fm_work2,fm_work1,ncol=homo)
             CALL cp_fm_column_scale(fm_work1,d1mdk_i(1:homo))
             CALL cp_fm_column_scale(fm_work1,d1mdk_ii(1:homo))
             CALL  cp_fm_scale_and_add(1.0_dp,fm_work3,1.0_dp,fm_work1,error=error)

             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(idir2)%matrix,&
                  fm_work2,fm_work1,ncol=homo,alpha=1.0_dp,error=error)
             CALL cp_fm_column_scale(fm_work1,d1mdk_ii(1:homo))
             CALL  cp_fm_scale_and_add(1.0_dp,fm_work3,1.0_dp,fm_work1,error=error)

             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(iii_B)%matrix,&
                  fm_work2,fm_work1,ncol=homo,alpha=1.0_dp,error=error)
             CALL cp_fm_column_scale(fm_work1,d1mdk_i(1:homo))
             CALL  cp_fm_scale_and_add(1.0_dp,fm_work3,1.0_dp,fm_work1,error=error)

             CALL cp_fm_to_fm(nmr_env%p_psi0(ispin,idir)%matrix,&
                  fm_work1,ncol=homo)
             CALL cp_fm_trace(fm_work1,fm_work3,contrib,error=error)
             nmr_env%chi_tensor(i_B,idir3) = nmr_env%chi_tensor(i_B,idir3) +&
                      scale_fac*contrib*alpha


      write(*,*) idir3, idir2, idir, contrib 

! second term psi0 times pp2
             CALL cp_fm_set_all(fm_work3,0.0_dp,error=error)
             CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,fm_work2,&
                  fm_work1,ncol=homo,alpha=1.0_dp,error=error)
             
             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(imom)%matrix,&
                  fm_work1,fm_work3,ncol=homo,alpha=1.0_dp,error=error)
             
             CALL cp_fm_to_fm(fm_work1,fm_work2,ncol=homo)
             CALL cp_fm_column_scale(fm_work2,d1mdk_i(1:homo))
             CALL cp_fm_column_scale(fm_work2,d1mdk_ii(1:homo))
             CALL  cp_fm_scale_and_add(1.0_dp,fm_work3,1.0_dp,fm_work2,error=error)

             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(idir2)%matrix,&
                  fm_work1,fm_work2,ncol=homo,alpha=1.0_dp,error=error)
             CALL cp_fm_column_scale(fm_work2,d1mdk_ii(1:homo))
             CALL  cp_fm_scale_and_add(1.0_dp,fm_work3,+1.0_dp,fm_work2,error=error)

             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(iii_B)%matrix,&
                  fm_work1,fm_work2,ncol=homo,alpha=1.0_dp,error=error)
             CALL cp_fm_column_scale(fm_work2,d1mdk_i(1:homo))
             CALL cp_fm_scale_and_add(1.0_dp,fm_work3,+1.0_dp,fm_work2,error=error)

             CALL cp_fm_to_fm(mo_coeff,fm_work1,ncol=homo)
             CALL cp_fm_trace(fm_work1,fm_work3,contrib,error=error)
             nmr_env%chi_tensor(i_B,idir3) = nmr_env%chi_tensor(i_B,idir3) +&
                      scale_fac*contrib*alpha

      write(*,*) idir3, idir2, idir, contrib 

             imom = ind_m2(idir2,ii_B)
             d1mdk_ii(1:homo) = d1mdk(ii_B,1:homo)
! third term ppsi0 times p3
             CALL cp_fm_set_all(fm_work1,0.0_dp,error=error)
             CALL cp_fm_set_all(fm_work2,0.0_dp,error=error)
             CALL cp_fm_set_all(fm_work3,0.0_dp,error=error)
             CALL cp_fm_to_fm(nmr_env%psi1_p(ispin,iii_B)%matrix,fm_work2,ncol=homo)
             IF ( .NOT. uni_occ ) THEN
               CALL cp_fm_column_scale(fm_work2,occupation(1:homo))
             END IF
             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(imom)%matrix,&
                  fm_work2,fm_work3,ncol=homo,alpha=-1.0_dp,error=error)

             CALL cp_fm_to_fm(fm_work2,fm_work1,ncol=homo)
             CALL cp_fm_column_scale(fm_work1,d1mdk_i(1:homo))
             CALL cp_fm_column_scale(fm_work1,d1mdk_ii(1:homo))
             CALL  cp_fm_scale_and_add(1.0_dp,fm_work3,-1.0_dp,fm_work1,error=error)

             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(idir2)%matrix,&
                  fm_work2,fm_work1,ncol=homo,alpha=1.0_dp,error=error)
             CALL cp_fm_column_scale(fm_work1,d1mdk_ii(1:homo))
             CALL  cp_fm_scale_and_add(1.0_dp,fm_work3,-1.0_dp,fm_work1,error=error)

             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(ii_B)%matrix,&
                  fm_work2,fm_work1,ncol=homo,alpha=1.0_dp,error=error)
             CALL cp_fm_column_scale(fm_work1,d1mdk_i(1:homo))
             CALL  cp_fm_scale_and_add(1.0_dp,fm_work3,-1.0_dp,fm_work1,error=error)

             CALL cp_fm_to_fm(nmr_env%p_psi0(ispin,idir)%matrix,&
                  fm_work1,ncol=homo)
             CALL cp_fm_trace(fm_work1,fm_work3,contrib,error=error)
             nmr_env%chi_tensor(i_B,idir3) = nmr_env%chi_tensor(i_B,idir3) +&
                      scale_fac*contrib*alpha

      write(*,*)  idir3, idir2, idir, contrib 
! second term psi0 times pp3
             CALL cp_fm_set_all(fm_work3,0.0_dp,error=error)
             CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,fm_work2,&
                  fm_work1,ncol=homo,alpha=1.0_dp,error=error)
             
             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(imom)%matrix,&
                  fm_work1,fm_work3,ncol=homo,alpha=-1.0_dp,error=error)
             
             CALL cp_fm_to_fm(fm_work1,fm_work2,ncol=homo)
             CALL cp_fm_column_scale(fm_work2,d1mdk_i(1:homo))
             CALL cp_fm_column_scale(fm_work2,d1mdk_ii(1:homo))
             CALL  cp_fm_scale_and_add(1.0_dp,fm_work3,-1.0_dp,fm_work2,error=error)

             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(idir2)%matrix,&
                  fm_work1,fm_work2,ncol=homo,alpha=1.0_dp,error=error)
             CALL cp_fm_column_scale(fm_work2,d1mdk_ii(1:homo))
             CALL  cp_fm_scale_and_add(1.0_dp,fm_work3,-1.0_dp,fm_work2,error=error)

             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(ii_B)%matrix,&
                  fm_work1,fm_work2,ncol=homo,alpha=1.0_dp,error=error)
             CALL cp_fm_column_scale(fm_work2,d1mdk_i(1:homo))
             CALL cp_fm_scale_and_add(1.0_dp,fm_work3,-1.0_dp,fm_work2,error=error)

             CALL cp_fm_to_fm(mo_coeff,fm_work1,ncol=homo)
             CALL cp_fm_trace(fm_work1,fm_work3,contrib,error=error)
             nmr_env%chi_tensor(i_B,idir3) = nmr_env%chi_tensor(i_B,idir3) +&
                      scale_fac*contrib*alpha

      write(*,*) idir3, idir2, idir, contrib 

         END IF
       END DO  ! idir2

    END DO  ! idir

    write(*,*) 'chi ', nmr_env%chi_tensor(i_B,1:3)
 
! stop 'vediamo'


    CALL cp_fm_release(fm_work1,error=error)
    CALL cp_fm_release(fm_work2,error=error)
    CALL cp_fm_release(fm_work3,error=error)

    DEALLOCATE(d1mdk,d1mdk_i,d1mdk_ii,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
  
  END SUBROUTINE chi_soft_analytic

!  *****************************************************************************************

  SUBROUTINE chi_soft_analytic_b(nmr_env,qs_env,psi1,p_psi1,ispin,i_B,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: psi1, p_psi1

    INTEGER, INTENT(IN)                      :: ispin
    INTEGER, INTENT(IN)                      :: i_B
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'chi_soft_analytic', &
      routineP = moduleN//':'//routineN


    INTEGER :: homo, i,j,iao, idir, idir2, idir3, ii_B, iii_B, imom, &
               istat, istate, nao, nmo 
    LOGICAL :: failure, gapw, uni_occ
    REAL(dp) :: alpha, contrib,d1dk(3),d1(3), dk(3),maxocc,scale_fac
    REAL(dp), DIMENSION(:), POINTER          :: occupation
    REAL(KIND=dp), DIMENSION(:), POINTER     :: d1mdk_i,d1mdk_ii
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: d1mdk
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_fm_type), POINTER                :: fm_work1,fm_work2
    TYPE(cp_fm_type), POINTER                :: fm_work3
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos


    failure = .FALSE.

    NULLIFY(cell,dft_control,linres_control,mo_coeff,mos,occupation,para_env)
    NULLIFY(fm_work1,fm_work2,fm_work3)

    CALL set_vecp(i_B,ii_B,iii_B)

    CALL get_qs_env(qs_env=qs_env,cell=cell,dft_control=dft_control,&
         linres_control=linres_control,&
         mos=mos,para_env=para_env)
    gapw = (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc)
    CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,&
         occupation_numbers=occupation, homo=homo, nmo=nmo,&
         nao=nao, uniform_occupation=uni_occ,maxocc=maxocc)

    ALLOCATE(d1mdk(3,homo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(d1mdk_i(homo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(d1mdk_ii(homo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
  
    d1 = linres_control%localized_wfn_control%centers_set(1)%array(1:3,1)
    DO istate = 1,homo
      dk = linres_control%localized_wfn_control%centers_set(ispin)%array(1:3,istate)
      d1dk = pbc(dk,d1,cell)
      d1mdk(1:3,istate) = d1dk(1:3)
    END DO 

    ! create a new matrix 
    NULLIFY(tmp_fm_struct)
    CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
           ncol_global=homo,para_env=para_env,context=mo_coeff%matrix_struct%context)
    CALL cp_fm_create (fm_work1, tmp_fm_struct )
    CALL cp_fm_set_all(fm_work1,0.0_dp,error=error)
!    CALL cp_fm_create (fm_work2, tmp_fm_struct )
!    CALL cp_fm_set_all(fm_work2,0.0_dp,error=error)
!    CALL cp_fm_create (fm_work3, tmp_fm_struct )
!    CALL cp_fm_set_all(fm_work3,0.0_dp,error=error)
    CALL cp_fm_struct_release ( tmp_fm_struct )

! Generate right part no gradient (the left part is r x p_psi0)

!   rxp term
    CALL cp_fm_to_fm(nmr_env%psi1_rxp(ispin,i_B)%matrix,psi1(ispin)%matrix)
!   (dl-dk)xp term
    IF(nmr_env%full_nmr) THEN
        CALL  cp_fm_scale_and_add(1.0_dp,psi1(ispin)%matrix,&
              -1.0_dp,nmr_env%psi1_D(ispin,i_B)%matrix,error=error)
    END IF
!   (r-dk)xp term
    CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(ii_B)%matrix,&
         nmr_env%psi1_p(ispin,iii_B)%matrix,psi1(ispin)%matrix,ncol=homo,&
         alpha=-1.0_dp,beta=1.0_dp,error=error)
    CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(iii_B)%matrix,&
         nmr_env%psi1_p(ispin,ii_B)%matrix,psi1(ispin)%matrix,ncol=homo,&
         alpha=1.0_dp,beta=1.0_dp,error=error)
    ! ii_B x iii_B component
    d1mdk_ii(1:homo) = -d1mdk(ii_B,1:homo)
    CALL cp_fm_to_fm(nmr_env%psi1_p(ispin,iii_B)%matrix,fm_work1,ncol=homo)
    CALL cp_fm_column_scale(fm_work1,d1mdk_ii(1:homo))
    CALL  cp_fm_scale_and_add(1.0_dp,psi1(ispin)%matrix,1.0_dp,fm_work1,error=error)
    ! iii_B x ii_B component
    d1mdk_ii(1:homo) = d1mdk(iii_B,1:homo)
    CALL cp_fm_to_fm(nmr_env%psi1_p(ispin,ii_B)%matrix,fm_work1,ncol=homo)
    CALL cp_fm_column_scale(fm_work1,d1mdk_ii(1:homo))
    CALL  cp_fm_scale_and_add(1.0_dp,psi1(ispin)%matrix,1.0_dp,fm_work1,error=error)

!  Right part no gradient done

    IF ( .NOT. uni_occ ) THEN
       alpha = 1.0_dp
       CALL cp_fm_column_scale(psi1(ispin)%matrix,occupation(1:homo))
    ELSE
       alpha = maxocc
    END IF



!     Lets say that idir gives the j_rho component and compute the relative contributions
    DO idir = 1,3


! Generate right part with gradient
       CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,psi1(ispin)%matrix,&
              p_psi1(ispin)%matrix,ncol=homo,alpha=+1.0_dp,error=error)
       ! These additional terms should give zero contribution due to
       ! the orthogonality of the response functions with respect to the unperturbed set
       IF(idir==ii_B) THEN
         CALL  cp_fm_scale_and_add(1.0_dp,p_psi1(ispin)%matrix,-1.0_dp,&
               nmr_env%psi1_p(ispin,iii_B)%matrix,error=error)
       ELSEIF(idir==iii_B)THEN
         CALL  cp_fm_scale_and_add(1.0_dp,p_psi1(ispin)%matrix,+1.0_dp,&
               nmr_env%psi1_p(ispin,ii_B)%matrix,error=error)
       END IF

!  Right part with gradient done
       IF ( .NOT. uni_occ ) THEN
          CALL cp_fm_column_scale(p_psi1(ispin)%matrix,occupation(1:homo))
       END IF

!dbg
! write(*,*) 'P_PSI1 dir=', idir
! do i = 1,nao
!     write(*,'(I4,4f12.6)') i,(p_psi1(ispin)%matrix%local_data(i,j),j=1,4)
! enddo
! write(*,*)
!dbg

!      idir2 is the component of the position operator
       DO idir2 = 1,3
         IF(idir /= idir2) THEN

!          idir3 is the component of CHI to which this term contributes
             CALL set_vecp_rev(idir,idir2,idir3)
             scale_fac=fac_vecp(idir3,idir2,idir)
             d1mdk_i(1:homo) = d1mdk(idir2,1:homo)

! first term p_psi0 times psi1
             CALL cp_fm_to_fm(nmr_env%p_psi0(ispin,idir)%matrix,fm_work1,ncol=homo)
             CALL cp_fm_column_scale(fm_work1,d1mdk_i(1:homo))
             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(idir2)%matrix,&
                  nmr_env%p_psi0(ispin,idir)%matrix,&
                  fm_work1,ncol=homo,alpha=+1.0_dp,beta=1.0_dp,error=error)
!dbg
! write(*,*) 'r PSI1 dir=', idir2
! do i = 1,nao
!     write(*,'(I4,4f12.6)') i,(fm_work1%local_data(i,j),j=1,4)
! enddo
! write(*,*)
!dbg
             CALL cp_fm_trace(fm_work1,psi1(ispin)%matrix,contrib,error=error)
             nmr_env%chi_tensor(i_B,idir3) = nmr_env%chi_tensor(i_B,idir3) +&
                      scale_fac*contrib*alpha

      write(*,'(A,3I4,3f12.6)') 'r p_psi0 psi1 ', idir3, idir2 , idir, scale_fac, contrib , scale_fac*contrib*alpha
! stop 'chi cont'



! second term psi0 times p_psi1
             CALL cp_fm_to_fm(mo_coeff,fm_work1,ncol=homo)
             CALL cp_fm_column_scale(fm_work1,d1mdk_i(1:homo))
             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(idir2)%matrix,&
                  mo_coeff,&
                  fm_work1,ncol=homo,alpha=1.0_dp,beta=1.0_dp,error=error)

             CALL cp_fm_trace(fm_work1,p_psi1(ispin)%matrix,contrib,error=error)
             nmr_env%chi_tensor(i_B,idir3) = nmr_env%chi_tensor(i_B,idir3) +&
                      scale_fac*contrib*alpha

      write(*,'(A,3I4,3f12.6)') 'r psi0 p_psi1 ', idir3, idir2 , idir, scale_fac, contrib , scale_fac*contrib*alpha

         END IF
       END DO  ! idir2

    END DO  ! idir

    write(*,*) 'chi ', nmr_env%chi_tensor(i_B,1:3)
 
! stop 'vediamo'

    CALL cp_fm_release(fm_work1,error=error)
!    CALL cp_fm_release(fm_work2,error=error)
!    CALL cp_fm_release(fm_work3,error=error)

    DEALLOCATE(d1mdk,d1mdk_i,d1mdk_ii,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
  
  END SUBROUTINE chi_soft_analytic_b

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

  SUBROUTINE interpolate_shift_pwgrid(nmr_env,particle_set,cell,shift_pw_rspace,&
             i_B,idir,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type)                         :: shift_pw_rspace
    INTEGER, INTENT(IN)                      :: i_B, idir
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'interpolate_shift_pwgrid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, natom
    LOGICAL                                  :: failure
    REAL(dp)                                 :: R_iatom(3), ra(3), shift_val

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(shift_pw_rspace%pw),cp_failure_level,routineP,error,failure)
    IF(.NOT. failure) THEN
      natom = SIZE(particle_set,1)

      DO iatom = 1,natom
        R_iatom = pbc(particle_set(iatom)%r,cell)
        shift_val = Eval_Interp_Spl3_pbc(R_iatom,shift_pw_rspace%pw,error)
        nmr_env%chemical_shift(i_B,idir,iatom)= nmr_env%chemical_shift(i_B,idir,iatom)+&
            nmr_env%shift_factor*shift_val
      END DO 

      IF(nmr_env%do_nics) THEN

        DO iatom = 1,nmr_env%n_nics
          ra(1:3) = nmr_env%r_nics(1:3,iatom)
          R_iatom = pbc(ra,cell)
          shift_val = Eval_Interp_Spl3_pbc(R_iatom,shift_pw_rspace%pw,error)
          nmr_env%chemical_shift_nics(i_B,idir,iatom) = &
             nmr_env%chemical_shift_nics(i_B,idir,iatom) + &
             nmr_env%shift_factor*shift_val
        END DO
      END IF

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



