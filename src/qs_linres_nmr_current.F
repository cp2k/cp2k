!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****m* cp2k/qs_linres_methods *
!!
!!   NAME 
!!     qs_linres_nmr_current
!!
!!   FUNCTION
!!     given the response wavefunctions obtained by the application
!!     of the (rxp), p, and ((dk-dl)xp) operators,
!!     here the current density vector (jx, jy, jz)
!!     is computed for the 3 directions of the magnetic field (Bx, By, Bz)   
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

MODULE qs_linres_nmr_current



  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE coefficient_types,               ONLY: coeff_zero,&
                                             coeff_type
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_column_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_transpose
  USE cp_fm_types,                     ONLY: cp_fm_get_submatrix,&
                                             cp_fm_p_type,&
                                             cp_fm_type,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply,&
                                             cp_sm_plus_fm_fm_t
  USE kinds,                           ONLY: dp
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_linres_types,                 ONLY: nmr_env_type
  USE qs_linres_nmr_op,                ONLY: set_vecp
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE sparse_matrix_types,             ONLY: real_matrix_type,&
                                             set_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop

#include "cp_common_uses.h"

  IMPLICIT NONE


  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: nmr_response_current

  CHARACTER(LEN=*), PARAMETER :: moduleN='qs_linres_nmr_current'

!!***
! *****************************************************************************

CONTAINS


  SUBROUTINE nmr_response_current(nmr_env,qs_env,psi1,p_psi1,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: psi1, p_psi1
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_response_current', &
      routineP = moduleN//'/'//routineN

    INTEGER                                  :: homo, i_B, ii_B, iii_B, iao,&
                                                idir, ispin, istate, nmo, nspins
    LOGICAL                                  :: failure, gapw, uni_occ
    REAL(dp)                                 :: alpha, maxocc
    REAL(dp)                                 :: dk(3),rmu(3),rmu_dk(3)
    REAL(dp), DIMENSION(:), POINTER          :: occupation
    REAL(dp), DIMENSION(:, :), POINTER       :: vecbuf_psi1

    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type), POINTER                :: rho_rspace, rho_gspace
    TYPE(cp_2d_r_p_type), DIMENSION(3)       :: vecbuf
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(real_matrix_type), POINTER          :: density_matrix
    TYPE(real_matrix_type), POINTER          :: density_matrix_2_ii
    TYPE(real_matrix_type), POINTER          :: density_matrix_2_iii


    failure = .FALSE.
    gapw=.FALSE.
    NULLIFY(rho_rspace,rho_gspace)
    NULLIFY(density_matrix,density_matrix_2_ii,density_matrix_2_iii)
    NULLIFY(cell,dft_control,mos,rho_struct)
    nspins = SIZE(psi1,1)
   
   
    CALL get_qs_env(qs_env=qs_env,rho=rho_struct,&
         cell=cell, dft_control=dft_control,mos=mos)

! Loop on the field direction 
    DO i_B = 1,3
       DO ispin = 1,dft_control%nspins
         CALL cp_fm_set_all(psi1(ispin)%matrix,0.0_dp,error=error)
         CALL cp_fm_set_all(p_psi1(ispin)%matrix,0.0_dp,error=error)
       END DO

       CALL set_vecp(i_B,ii_B,iii_B)
       DO ispin = 1,nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,&
              occupation_numbers=occupation, homo=homo, nmo=nmo,&
              uniform_occupation=uni_occ,maxocc=maxocc)

         !First contribution from the response psi1_p
         DO istate = 1,nmr_env%nstates(ispin)
           CALL cp_fm_get_submatrix(nmr_env%psi1_p(ispin,ii_B)%matrix,&
                vecbuf(ii_B)%array,1,istate,nmr_env%nao,1,transpose=.TRUE.,&
                error=error)
           CALL cp_fm_get_submatrix(nmr_env%psi1_p(ispin,iii_B)%matrix,&
                vecbuf(iii_B)%array,1,istate,nmr_env%nao,1,transpose=.TRUE.,&
                error=error)
           dk(1:3) = nmr_env%centers_set(ispin)%array(1:3,istate)
           DO iao = 1,nmr_env%nao
             Rmu(1:3) =  nmr_env%basisfun_center(1:3,iao)
             Rmu_dk = pbc(dk,Rmu,cell)
             vecbuf_psi1(1,iao) = Rmu_dk(ii_B)*vecbuf(iii_B)%array(1,iao)-&
                                  Rmu_dk(iii_B)*vecbuf(ii_B)%array(1,iao)
           END DO  ! iao
           ! Copy the vector in the full matrix psi1
           CALL cp_fm_set_submatrix(psi1(ispin)%matrix,vecbuf_psi1,&
                1,istate,nmr_env%nao,1,transpose=.TRUE.,error=error)
         END DO  ! istate

         !Add to psi the psi1_rxp contribution
         CALL  cp_fm_scale_and_add(-1.0_dp,psi1(ispin)%matrix,1.0_dp,nmr_env%psi1_rxp(ispin,i_B)%matrix,error=error)
         IF(nmr_env%full_nmr) THEN
            CALL  cp_fm_scale_and_add(1.0_dp,psi1(ispin)%matrix,-1.0_dp,nmr_env%psi1_D(ispin,i_B)%matrix,error=error)
         END IF

         DO idir = 1,3
         ! Apply the p operator
           CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,psi1(ispin)%matrix,&
                p_psi1(ispin)%matrix,ncol=nmo,alpha=-1.0_dp,error=error)
           density_matrix => nmr_env%jp1_ao(idir,ispin,i_B)%matrix
           CALL set_matrix(density_matrix,0.0_dp)

         ! Contribution to the current density matrix that does not depend 
         !   on the position r where J(r) is computed
           IF ( .NOT. uni_occ ) THEN
             alpha = 1.0_dp
             CALL cp_fm_column_scale(psi1(ispin)%matrix,occupation(1:homo))
             CALL cp_fm_column_scale(p_psi1(ispin)%matrix,occupation(1:homo))
           ELSE
             alpha = maxocc
           END IF
           CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix,&
                            matrix_v=nmr_env%p_psi0(ispin,idir)%matrix,&
                            matrix_g=psi1(ispin)%matrix,ncol=homo,&
                            alpha=alpha)

           CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix,&
                            matrix_v=mo_coeff,&
                            matrix_g=p_psi1(ispin)%matrix,ncol=homo,&
                            alpha=-alpha)
          
           ! For the position dependent term two more density matrix are needed
           density_matrix_2_ii => nmr_env%jp2_ao(1,ispin)%matrix
           CALL set_matrix(density_matrix_2_ii,0.0_dp)
           density_matrix_2_iii => nmr_env%jp2_ao(2,ispin)%matrix
           CALL set_matrix(density_matrix_2_iii,0.0_dp)

           CALL cp_fm_to_fm(nmr_env%psi1_p(ispin,iii_B)%matrix,psi1(ispin)%matrix)
           CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,psi1(ispin)%matrix,&
                p_psi1(ispin)%matrix,ncol=nmo,alpha=-1.0_dp,error=error)
           IF ( .NOT. uni_occ ) THEN
             alpha = 1.0_dp
             CALL cp_fm_column_scale(psi1(ispin)%matrix,occupation(1:homo))
             CALL cp_fm_column_scale(p_psi1(ispin)%matrix,occupation(1:homo))
           ELSE
             alpha = maxocc
           END IF
           CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix_2_ii,&
                            matrix_v=nmr_env%p_psi0(ispin,idir)%matrix,&
                            matrix_g=psi1(ispin)%matrix,ncol=homo,&
                            alpha=alpha)

           CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix_2_ii,&
                            matrix_v=mo_coeff,&
                            matrix_g=p_psi1(ispin)%matrix,ncol=homo,&
                            alpha=-alpha)

           CALL cp_fm_to_fm(nmr_env%psi1_p(ispin,ii_B)%matrix,psi1(ispin)%matrix)
           CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,psi1(ispin)%matrix,&
                p_psi1(ispin)%matrix,ncol=nmo,alpha=-1.0_dp,error=error)
           IF ( .NOT. uni_occ ) THEN
             alpha = 1.0_dp
             CALL cp_fm_column_scale(psi1(ispin)%matrix,occupation(1:homo))
             CALL cp_fm_column_scale(p_psi1(ispin)%matrix,occupation(1:homo))
           ELSE
             alpha = maxocc
           END IF
           CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix_2_iii,&
                            matrix_v=nmr_env%p_psi0(ispin,idir)%matrix,&
                            matrix_g=psi1(ispin)%matrix,ncol=homo,&
                            alpha=alpha)

           CALL cp_sm_plus_fm_fm_t(sparse_matrix=density_matrix_2_iii,&
                            matrix_v=mo_coeff,&
                            matrix_g=p_psi1(ispin)%matrix,ncol=homo,&
                            alpha=-alpha)

        ! Calculate the current density on the pw grid (only soft if GAPW)
        ! This is vector component idir of the response current density 
        ! generated by the magnetic field in cartesian direction i_B
        ! Use the qs_rho_type already  used for rho during the scf

            rho_rspace => rho_struct%rho_r(ispin)
            rho_gspace => rho_struct%rho_g(ispin)

 !           CALL calculate_jrho_resp(density_matrix, density_matrix_2_ii,&
 !                density_matrix_2_iii, rho_rspace, rho_gspace,qs-env,&
 !                gapw, error=error)

            IF(gapw) THEN
            ! compute the atomic response current densities on the spherical grids
            END IF

         END DO  ! idir
         
       END DO  !  ispin
    END DO

  END SUBROUTINE nmr_response_current
 
END MODULE qs_linres_nmr_current
