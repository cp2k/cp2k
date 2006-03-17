
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
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE kinds,                           ONLY: dp
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


    INTEGER :: homo, iao, idir, idir2, idir3, ii_B, iii_B, imom, istat, istate, nao, nmo 
    LOGICAL :: failure, gapw
    REAL(dp) :: contrib,dk(3),scale_fac
    REAL(dp), DIMENSION(:), POINTER          :: occupation
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vec_psi0
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vec_p_psi0
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_fm_type), POINTER                :: vec1,p1,p2,p3
    TYPE(cp_fm_type), POINTER                :: vec2,pp1,pp2,pp3
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos


    failure = .FALSE.

    NULLIFY(dft_control,linres_control,mo_coeff,mos,occupation,para_env)
    NULLIFY(vec1,p1,p2,p3,vec2,pp1,pp2,pp3,vecbuffer)

    CALL set_vecp(i_B,ii_B,iii_B)
    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,linres_control=linres_control,&
         mos=mos,para_env=para_env)
    gapw = (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc)
    CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,&
         occupation_numbers=occupation, homo=homo, nmo=nmo)

    CALL cp_fm_to_fm(nmr_env%psi1_rxp(ispin,i_B)%matrix,psi1(ispin)%matrix)
    IF(nmr_env%full_nmr) THEN
        CALL  cp_fm_scale_and_add(1.0_dp,psi1(ispin)%matrix,&
              -1.0_dp,nmr_env%psi1_D(ispin,i_B)%matrix,error=error)
    END IF

    ALLOCATE(vecbuffer(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(vec_psi0(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(vec_p_psi0(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
  
    ! create a new matrix structure nao x 1
    NULLIFY(tmp_fm_struct)
    CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
           ncol_global=1,para_env=para_env,context=mo_coeff%matrix_struct%context)
    CALL cp_fm_create (vec1, tmp_fm_struct )
    CALL cp_fm_create (p1, tmp_fm_struct )
    CALL cp_fm_create (p2, tmp_fm_struct )
    CALL cp_fm_create (p3, tmp_fm_struct )
    CALL cp_fm_create (vec2, tmp_fm_struct )
    CALL cp_fm_create (pp1, tmp_fm_struct )
    CALL cp_fm_create (pp2, tmp_fm_struct )
    CALL cp_fm_create (pp3, tmp_fm_struct )
    CALL cp_fm_struct_release ( tmp_fm_struct )

    DO istate = 1,homo

!     Extract the vector istate from the different set of mos
      CALL cp_fm_get_submatrix(mo_coeff,vec_psi0,1,istate,nao,1,transpose=.TRUE.,error=error)

      CALL cp_fm_get_submatrix(psi1(ispin)%matrix,vecbuffer,1,&
           istate,nao,1,transpose=.TRUE.,error=error)
      CALL cp_fm_set_submatrix(p1,vecbuffer,1,1,nao,1,transpose=.TRUE.,error=error)

      CALL cp_fm_get_submatrix(nmr_env%psi1_p(ispin,ii_B)%matrix,vecbuffer,1,&
           istate,nao,1,transpose=.TRUE.,error=error)
      CALL cp_fm_set_submatrix(p2,vecbuffer,1,1,nao,1,transpose=.TRUE.,error=error)
      CALL cp_fm_get_submatrix(nmr_env%psi1_p(ispin,iii_B)%matrix,vecbuffer,1,&
           istate,nao,1,transpose=.TRUE.,error=error)
      CALL cp_fm_set_submatrix(p3,vecbuffer,1,1,nao,1,transpose=.TRUE.,error=error)

    ! Calculate the x y z x2 xy xz y2 yz z2 operators (in this order)
      dk =  linres_control%localized_wfn_control%centers_set(ispin)%array(1:3,istate)
      DO idir = 1,9
        CALL set_matrix(nmr_env%op_mom_ao(idir)%matrix,0.0_dp)
      END DO
      CALL rRc_xyz_ao(nmr_env%op_mom_ao,qs_env,dk,order=2,&
           minimum_image=.TRUE.,soft=gapw,error=error)

!     Lets say that idir gives the j_rho component and compute the relative contributions
      DO idir = 1,3

!       Extract the vector istate from the component idir of p_psi0 
        CALL cp_fm_get_submatrix(nmr_env%p_psi0(ispin,idir)%matrix,vec_p_psi0,1,&
             istate,nao,1,transpose=.TRUE.,error=error)
  
      ! Apply the component idir of the p operator
         CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,p1,&
              pp1,ncol=1,alpha=-1.0_dp,error=error)
         IF(idir==ii_B) THEN
           CALL  cp_fm_scale_and_add(1.0_dp,pp1,-1.0_dp,p3,error=error)
         ELSEIF(idir==iii_B)THEN
           CALL  cp_fm_scale_and_add(1.0_dp,pp1,+1.0_dp,p2,error=error)
         END IF

         CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,p2,&
              pp2,ncol=1,alpha=-1.0_dp,error=error)

         CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,p3,&
              pp3,ncol=1,alpha=-1.0_dp,error=error)

!        idir2 is the component of the position operator
         DO idir2 = 1,3
           IF(idir /= idir2) THEN

  !          idir3 is the component of CHI to which this term contribute
             CALL set_vecp_rev(idir,idir2,idir3)
             scale_fac=fac_vecp(idir3,idir2,idir)

             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(idir2)%matrix,p1,&
                  vec1,ncol=1,alpha=1.0_dp,error=error)
             imom = ind_m2(idir2,iii_B)
             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(imom)%matrix,p2,&
                  vec2,ncol=1,alpha=1.0_dp,error=error)
             CALL  cp_fm_scale_and_add(1.0_dp,vec1,+1.0_dp,vec2,error=error)
             imom = ind_m2(idir2,ii_B)
             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(imom)%matrix,p3,&
                  vec2,ncol=1,alpha=1.0_dp,error=error)
             CALL  cp_fm_scale_and_add(1.0_dp,vec1,-1.0_dp,vec2,error=error)
             CALL cp_fm_get_submatrix(vec1,vecbuffer,1,1,nao,1,transpose=.TRUE.,error=error)
             contrib = 0.0_dp
             DO iao=1,nao
               contrib = contrib+vecbuffer(1,iao)*vec_p_psi0(1,iao)
             end do 
             nmr_env%chi_tensor(i_B,idir3) = nmr_env%chi_tensor(i_B,idir3) +&
                      scale_fac*contrib*occupation(istate)

             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(idir2)%matrix,pp1,&
                  vec1,ncol=1,alpha=-1.0_dp,error=error)
             imom = ind_m2(idir2,iii_B)
             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(imom)%matrix,pp2,&
                  vec2,ncol=1,alpha=1.0_dp,error=error)
             CALL  cp_fm_scale_and_add(1.0_dp,vec1,-1.0_dp,vec2,error=error)
             imom = ind_m2(idir2,ii_B)
             CALL cp_sm_fm_multiply(nmr_env%op_mom_ao(imom)%matrix,pp3,&
                  vec2,ncol=1,alpha=1.0_dp,error=error)
             CALL  cp_fm_scale_and_add(1.0_dp,vec1,1.0_dp,vec2,error=error)
             CALL cp_fm_get_submatrix(vec1,vecbuffer,1,1,nao,1,transpose=.TRUE.,error=error)
             contrib = 0.0_dp
             DO iao=1,nao
               contrib = contrib+vecbuffer(1,iao)*vec_psi0(1,iao)
             end do 
             nmr_env%chi_tensor(i_B,idir3) = nmr_env%chi_tensor(i_B,idir3) +&
                      scale_fac*contrib*occupation(istate)

           END IF
         END DO  ! idir2

      END DO  ! idir
    END DO  ! istate

    CALL cp_fm_release(vec1,error=error)
    CALL cp_fm_release(p1,error=error)
    CALL cp_fm_release(p2,error=error)
    CALL cp_fm_release(p3,error=error)
    CALL cp_fm_release(vec2,error=error)
    CALL cp_fm_release(pp1,error=error)
    CALL cp_fm_release(pp2,error=error)
    CALL cp_fm_release(pp3,error=error)

    DEALLOCATE(vecbuffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(vec_psi0,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(vec_p_psi0,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
  
  END SUBROUTINE chi_soft_analytic

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



