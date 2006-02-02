!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2005  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****m* cp2k/qs_linres_nmr_op
!!
!!   NAME
!!     qs_linres_nmr_op
!!
!!   FUNCTION
!!     Calculate the operators p rxp and D needed in the optimization
!!     of the different contribution of the firs order response orbitals
!!     in a nmr calculation
!!
!!   NOTE
!!     The interactions are considered only within the minimum image convention
!!  
!!   AUTHOR
!!      MI 
!!
!!   MODIFICATION HISTORY
!!      created 07-2005 [MI]
!!
!!   SOURCE
!!****

MODULE qs_linres_nmr_op
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_types,                     ONLY: cp_fm_get_submatrix,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_types,                 ONLY: linres_control_type,&
                                             nmr_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_operators_ao,                 ONLY: p_xyz_ao,&
                                             rmc_x_p_xyz_ao
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: nmr_operators, set_vecp

  CHARACTER(LEN=*), PARAMETER :: moduleN='qs_linres_nmr_utils'

!!***
! *****************************************************************************
CONTAINS
                                      
! *****************************************************************************
!!****f* qs_linres_nmr_op/nmr_operators
!!
!!   NAME
!!     nmr_operators
!!
!!   FUNCTION
!!     Calculate the first order hamiltonian applied to the ao
!!     and then apply them to the ground state orbitals, 
!!     the h1_psi1 full matrices are then ready to solve the
!!     non-homogeneous linear equations that give the psi1
!!     linear response orbitals. 
!!
!!   NOTES
!!     For the operators rxp and D the h1 depends on the psi0 to which
!!     is applied, or better the center of charge of the psi0 is 
!!     used to define the position operator
!!     The centers of the orbitals result form the orbital localization procedure
!!     that typicaly uses the berry phase operator to define the Wannier centers.
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     07.2005 created [MI]
!!
!****************************************************************************

  SUBROUTINE nmr_operators(nmr_env,qs_env,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_operators', &
      routineP = moduleN//'/'//routineN

    INTEGER                                  :: homo, iao, iatom, idir, &
                                                idir2, ii, iii, ispin, istat, &
                                                istate, nao, natom, nmo, &
                                                nsgf, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    LOGICAL                                  :: failure, ionode
    REAL(dp)                                 :: ck(3), ckdk(3), dk(3)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_r_p_type), DIMENSION(3)       :: vecbuf_cmd, vecbuf_p0, &
                                                vecbuf_rxp0
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: lr_section

    failure = .FALSE.

    NULLIFY(atomic_kind_set,cell,dft_control,linres_control,logger,mos,particle_set,&
         lr_section)
 
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)
    IF (output_unit>0) THEN

      WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
         "Calculation of the p and (r-d)xp operators applied to psi0"

    END IF
    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, cell=cell, &
                    dft_control=dft_control, linres_control=linres_control,&
                    mos=mos, particle_set=particle_set)

    ! Calculate the p operator
    CALL p_xyz_ao(nmr_env%op_p_ao,qs_env,minimum_image=.TRUE.,error=error)
    ! Multiply the p operator to the psi0
    DO idir = 1,3
      DO ispin = 1, dft_control%nspins
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nao=nao, nmo=nmo)
         CALL cp_fm_set_all(nmr_env%p_psi0(ispin,idir)%matrix,0.0_dp,error=error)
         CALL cp_sm_fm_multiply(nmr_env%op_p_ao(idir)%matrix,mo_coeff,&
              nmr_env%p_psi0(ispin,idir)%matrix,ncol=nmo,alpha=-1.0_dp,error=error)
      END DO 
    END DO 

    DO idir = 1,3 
      NULLIFY(vecbuf_p0(idir)%array,vecbuf_rxp0(idir)%array,vecbuf_cmd(idir)%array)
      ALLOCATE(vecbuf_p0(idir)%array(1,nao),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(vecbuf_rxp0(idir)%array(1,nao),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(vecbuf_cmd(idir)%array(1,nao),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END DO 

    ! Calculate the (r - dk)xp operator applied to psi0k
    ! We use the distributive property of the vector product and calculatr
    ! (r-c)xp + (c-d)xp 
    ! where c depends on the contracted functions and not on the states
    ! d is the center of a specific state and a loop over states is needed
    ! the second term can be added in a second moment as a correction
    ! notice: (r-c) and p are operators, whereas (c-d) is a multiplicative factor

    !First term: operator matrix elements
    CALL rmc_x_p_xyz_ao(nmr_env%op_rmd_ao,qs_env,minimum_image=.TRUE.,error=error)

    ! Multiply the (r-c)xp operator to the psi0
    DO idir = 1,3
      DO ispin = 1, dft_control%nspins
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nao=nao, nmo=nmo)
         CALL cp_fm_set_all(nmr_env%rxp_psi0(ispin,idir)%matrix,0.0_dp,error=error)
         CALL cp_sm_fm_multiply(nmr_env%op_rmd_ao(idir)%matrix,mo_coeff,&
              nmr_env%rxp_psi0(ispin,idir)%matrix,ncol=nmo,alpha=-1.0_dp,error=error)
      END DO 
    END DO 

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf)

    natom = SIZE(particle_set,1)
    ALLOCATE (first_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf(natom),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf,&
                          last_sgf=last_sgf)

    !Build 3 arrays where for each contracted basis function 
    !the x y and z coordinates of the center are given
    DO iatom = 1,natom
      DO iao = first_sgf(iatom),last_sgf(iatom)
        DO idir = 1,3
          vecbuf_cmd(idir)%array(1,iao) = particle_set(iatom)%r(idir)
          nmr_env%basisfun_center(idir,iao) = particle_set(iatom)%r(idir)
        END DO
      END DO
    END DO

    
    !Calculate the second term of the operator proceeding state by state
    DO ispin = 1, dft_control%nspins
       CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nao=nao, nmo=nmo,homo=homo)
      
! write(*,*) 'homo ', homo, ' nao ', nao
       DO istate = 1,homo
         !center of the localized psi0 state istate 
         dk(1:3) = linres_control%localized_wfn_control%centers_set(ispin)%array(1:3,istate)
! write(*,*) 'dk ', istate, dk
         !get the corresponding vector from the full matrix (p psi0)
         DO idir = 1,3
           CALL cp_fm_get_submatrix(nmr_env%p_psi0(ispin,idir)%matrix,vecbuf_p0(idir)%array,&
                1,istate,nao,1,transpose=.TRUE.,error=error)
! write(*,*) 'idir p0 ', idir
! write(*,*) vecbuf_p0(idir)%array(1,:)
         END DO
         DO idir = 1,3
           CALL set_vecp(idir,ii,iii)
           !get the full matrix vector containing only the first term, 
           !to which now the second term will be added
           CALL cp_fm_get_submatrix(nmr_env%rxp_psi0(ispin,idir)%matrix,vecbuf_rxp0(idir)%array,&
                1,istate,nao,1,transpose=.TRUE.,error=error)

!  This loop should be distributed over the processors
           DO iao = 1,nao
             DO idir2 = 1,3
               ck(idir2) = vecbuf_cmd(idir2)%array(1,iao)
!dbg
               ck(idir2) = nmr_env%basisfun_center(idir2,iao)
!dbg
             END DO
! write(*,*) 'iao', iao, 'ck ', ck
             ckdk = pbc(dk,ck,cell)
! write(*,*) 'c-dk ', ckdk
             !Add the second term element by element (not very efficient indeed)
! write(*,*) 'rxp0 terms',vecbuf_rxp0(idir)%array(1,iao),ckdk(ii)*vecbuf_p0(iii)%array(1,iao) - &
!                                ckdk(iii)*vecbuf_p0(ii)%array(1,iao)
             vecbuf_rxp0(idir)%array(1,iao) =  vecbuf_rxp0(idir)%array(1,iao) + & 
                                ckdk(ii)*vecbuf_p0(iii)%array(1,iao) - &
                                ckdk(iii)*vecbuf_p0(ii)%array(1,iao)
! write(*,*) 'rxp0 fin',vecbuf_rxp0(idir)%array(1,iao)
           END DO
           !Give back the modified vector to the full matrix
           CALL cp_fm_set_submatrix(nmr_env%rxp_psi0(ispin,idir)%matrix,vecbuf_rxp0(idir)%array,&
                1,istate,nao,1,transpose=.TRUE.,error=error)
         END DO 

       END DO  ! istate
    END DO  ! ispin

    DEALLOCATE (first_sgf,last_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)


    DO idir = 1,3
      DEALLOCATE(vecbuf_p0(idir)%array,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(vecbuf_rxp0(idir)%array,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(vecbuf_cmd(idir)%array,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END DO

  END SUBROUTINE nmr_operators

!****************************************************************************
  SUBROUTINE set_vecp(i1,i2,i3)

    INTEGER, INTENT(IN)                      :: i1
    INTEGER, INTENT(OUT)                     :: i2, i3

    IF(i1==1) THEN
      i2=2
      i3=3
    ELSEIF(i1==2) THEN
      i2=3
      i3=1
    ELSEIF(i1==3) THEN
      i2=1
      i3=2
    ELSE
    END IF

  END SUBROUTINE set_vecp 
!****************************************************************************
END MODULE qs_linres_nmr_op
