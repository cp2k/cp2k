!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate the operators p rxp and D needed in the optimization
!>      of the different contribution of the firs order response orbitals
!>      in a nmr calculation
!> \note
!>      The interactions are considered only within the minimum image convention
!> \par History
!>       created 07-2005 [MI]
!> \author MI
! *****************************************************************************
MODULE qs_linres_nmr_op
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
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
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE orbital_pointers,                ONLY: coset
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_types,                 ONLY: get_nmr_env,&
                                             linres_control_type,&
                                             nmr_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_operators_ao,                 ONLY: p_xyz_ao,&
                                             rmc_x_p_xyz_ao
  USE sparse_matrix_types,             ONLY: checksum_matrix,&
                                             real_matrix_p_type,&
                                             set_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: nmr_operators, set_vecp, set_vecp_rev, fac_vecp, ind_m2

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_nmr_op'

CONTAINS

! *****************************************************************************
!> \brief Calculate the first order hamiltonian applied to the ao
!>      and then apply them to the ground state orbitals,
!>      the h1_psi1 full matrices are then ready to solve the
!>      non-homogeneous linear equations that give the psi1
!>      linear response orbitals.
!> \note
!>      For the operators rxp and D the h1 depends on the psi0 to which
!>      is applied, or better the center of charge of the psi0 is
!>      used to define the position operator
!>      The centers of the orbitals result form the orbital localization procedure
!>      that typicaly uses the berry phase operator to define the Wannier centers.
!> \par History
!>      07.2005 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE nmr_operators(nmr_env,qs_env,error)

    TYPE(nmr_env_type)                       :: nmr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'nmr_operators', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homo, i, iao, iatom, &
                                                idir, ii, iii, ispin, istat, &
                                                istate, nao, natom, nmo, &
                                                nsgf, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: first_sgf, last_sgf
    LOGICAL                                  :: failure, ionode
    REAL(dp)                                 :: chk(3), ck(3), ckdk(3), dk(3)
    REAL(dp), DIMENSION(:, :), POINTER       :: basisfun_center, vecbuf_c0
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_2d_r_p_type), DIMENSION(3)       :: vecbuf_RmdC0
    TYPE(cp_2d_r_p_type), DIMENSION(:), &
      POINTER                                :: centers_set
    TYPE(cp_fm_p_type), DIMENSION(3)         :: fm_Rmd_mos
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: p_psi0, rxp_psi0
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: fm_work1, mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(linres_control_type), POINTER       :: linres_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: op_p_ao, op_rmd_ao
    TYPE(section_vals_type), POINTER         :: lr_section

    CALL timeset(routineN,"I"," ",handle)

    failure = .FALSE.

    NULLIFY(atomic_kind_set,cell,dft_control,linres_control)
    NULLIFY(logger,mos,particle_set,lr_section,fm_work1)
    NULLIFY(basisfun_center, centers_set, p_psi0, rxp_psi0, op_p_ao, op_rmd_ao)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    lr_section => section_vals_get_subs_vals(qs_env%input,&
                  "PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
                  extension=".linresLog",error=error)
    IF (output_unit>0) THEN
      WRITE (UNIT=output_unit,FMT="(T10,A,/)")&
         "Calculation of the p and (r-d)xp operators applied to psi0"
    END IF

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, cell=cell, &
                    dft_control=dft_control, linres_control=linres_control,&
                    mos=mos, para_env=para_env, &
                    particle_set=particle_set,error=error)

    CALL get_nmr_env(nmr_env=nmr_env, nao=nao, centers_set=centers_set,&
                     basisfun_center=basisfun_center, &
                     p_psi0=p_psi0, rxp_psi0=rxp_psi0,&
                     op_p_ao=op_p_ao, op_rmd_ao=op_rmd_ao,error=error)

    NULLIFY(vecbuf_c0)
    ALLOCATE(vecbuf_c0(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO idir = 1,3
      NULLIFY(vecbuf_Rmdc0(idir)%array)
      ALLOCATE(vecbuf_Rmdc0(idir)%array(1,nao),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END DO

    ! Calculate the (r - dk)xp operator applied to psi0k
    ! We use the distributive property of the vector product and calculatr
    ! (r-c)xp + (c-d)xp
    ! where c depends on the contracted functions and not on the states
    ! d is the center of a specific state and a loop over states is needed
    ! the second term can be added in a second moment as a correction
    ! notice: (r-c) and p are operators, whereas (c-d) is a multiplicative factor

!!$    !First term: operator matrix elements
!!$    CALL rmc_x_p_xyz_ao(op_rmd_ao,qs_env,minimum_image=.FALSE.,error=error)
!************************************************************
    IF(output_unit>0) THEN
       WRITE(output_unit,*) 'WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING'
       WRITE(output_unit,*) 'rmc_x_p_xyz_ao: L is computed for each states'
       WRITE(output_unit,*) 'WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING'
    ENDIF

    DO ispin=1,dft_control%nspins
       p_psi0(ispin,1)%matrix%local_data=0e0_dp
       p_psi0(ispin,2)%matrix%local_data=0e0_dp
       p_psi0(ispin,3)%matrix%local_data=0e0_dp
       CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nmo=nmo,homo=homo)
       DO istate = 1,homo    
          CALL set_matrix(op_rmd_ao(1)%matrix,0.0_dp)
          CALL set_matrix(op_rmd_ao(2)%matrix,0.0_dp)
          CALL set_matrix(op_rmd_ao(3)%matrix,0.0_dp)

          CALL rmc_x_p_xyz_ao(op_rmd_ao,qs_env,minimum_image=.FALSE.,&
               & wancen=centers_set(ispin)%array(1:3,istate),error=error)

          DO idir = 1,3
             CALL cp_fm_set_all(rxp_psi0(ispin,idir)%matrix,0.0_dp,error=error)
             CALL cp_sm_fm_multiply(op_rmd_ao(idir)%matrix,mo_coeff,&
                  rxp_psi0(ispin,idir)%matrix,ncol=nmo,alpha=-1.0_dp,error=error)
             CALL cp_fm_to_fm(rxp_psi0(ispin,idir)%matrix,&
                  &           p_psi0(ispin,idir)%matrix,1,istate,istate)
          END DO
       ENDDO
       rxp_psi0(ispin,1)%matrix%local_data = p_psi0(ispin,1)%matrix%local_data
       rxp_psi0(ispin,2)%matrix%local_data = p_psi0(ispin,2)%matrix%local_data
       rxp_psi0(ispin,3)%matrix%local_data = p_psi0(ispin,3)%matrix%local_data
    ENDDO
    !
    ! print checksums
    chk(1) = checksum_matrix(op_rmd_ao(1)%matrix,para_env)
    chk(2) = checksum_matrix(op_rmd_ao(2)%matrix,para_env)
    chk(3) = checksum_matrix(op_rmd_ao(3)%matrix,para_env)
    IF(output_unit>0) THEN
       WRITE(output_unit,'(A,E23.16)') 'nmr_operators: CheckSum L_x =',chk(1)
       WRITE(output_unit,'(A,E23.16)') 'nmr_operators: CheckSum L_y =',chk(2)
       WRITE(output_unit,'(A,E23.16)') 'nmr_operators: CheckSum L_z =',chk(3)
    ENDIF
    !
    ! Calculate the px py pz operators
    CALL p_xyz_ao(op_p_ao,qs_env,minimum_image=.FALSE.,error=error)
    !
    ! print checksums
    chk(1) = checksum_matrix(op_p_ao(1)%matrix,para_env)
    chk(2) = checksum_matrix(op_p_ao(2)%matrix,para_env)
    chk(3) = checksum_matrix(op_p_ao(3)%matrix,para_env)
    IF(output_unit>0) THEN
       WRITE(output_unit,'(A,E23.16)') 'nmr_operators: CheckSum P_x =',chk(1)
       WRITE(output_unit,'(A,E23.16)') 'nmr_operators: CheckSum P_y =',chk(2)
       WRITE(output_unit,'(A,E23.16)') 'nmr_operators: CheckSum P_z =',chk(3)
    ENDIF
    ! Apply the p operator to the psi0
    DO idir = 1,3
      DO ispin = 1, dft_control%nspins
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nmo=nmo)
         CALL cp_fm_set_all(p_psi0(ispin,idir)%matrix,0.0_dp,error=error)
         CALL cp_sm_fm_multiply(op_p_ao(idir)%matrix,mo_coeff,&
                p_psi0(ispin,idir)%matrix,ncol=nmo,alpha=-1.0_dp,error=error)
      END DO
    END DO
    !
    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)
    !
    ! Apply the (r-c)xp operator to the psi0
    !DO ispin = 1, dft_control%nspins
    !  CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nmo=nmo, homo=homo)
    !  DO idir = 1,3
    !     CALL cp_fm_set_all(rxp_psi0(ispin,idir)%matrix,0.0_dp,error=error)
    !     CALL cp_sm_fm_multiply(op_rmd_ao(idir)%matrix,mo_coeff,&
    !            rxp_psi0(ispin,idir)%matrix,ncol=nmo,alpha=-1.0_dp,error=error)
    !  END DO
    !END DO

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
          nmr_env%basisfun_center(idir,iao) = particle_set(iatom)%r(idir)
        END DO
      END DO
    END DO

    !Calculate the second term of the operator state by state
    DO ispin = 1, -1!dft_control%nspins
       CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, nmo=nmo,homo=homo)

      !   Allocate full matrices as working storage in the calculation
      !   of the rxp operator matrix. 3 matrices for the 3 Cartesian direction
      !   plus one to apply the momentum oprator to the modified mos fm
       NULLIFY(tmp_fm_struct)
       CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nao,&
            ncol_global=homo,para_env=para_env,context=mo_coeff%matrix_struct%context,error=error)
       DO idir = 1,3
          NULLIFY(fm_Rmd_mos(idir)%matrix)
          CALL cp_fm_create (fm_Rmd_mos(idir)%matrix,tmp_fm_struct,error=error)
       END DO
       CALL cp_fm_create (fm_work1,tmp_fm_struct,error=error)
       CALL cp_fm_struct_release ( tmp_fm_struct ,error=error)

   ! This part should be done better, using the full matrix distribution
       DO istate = 1,homo
         CALL cp_fm_get_submatrix(mo_coeff,vecbuf_c0,1,istate,nao,1,transpose=.TRUE.,error=error)
!!         !center of the localized psi0 state istate
         dk(1:3) = centers_set(ispin)%array(1:3,istate)
         DO idir = 1,3
!!!  This loop should be distributed over the processors
           DO iao = 1,nao
             ck(1:3) = basisfun_center(1:3,iao)
             ckdk = pbc(dk,ck,cell)
             vecbuf_Rmdc0(idir)%array(1,iao) = vecbuf_c0(1,iao)*ckdk(idir)
           END DO ! iao
           CALL cp_fm_set_submatrix(fm_Rmd_mos(idir)%matrix,vecbuf_Rmdc0(idir)%array,&
                1,istate,nao,1,transpose=.TRUE.,error=error)
         END DO  ! idir
       END DO  ! istate

       DO idir = 1,3
           CALL set_vecp(idir,ii,iii)

           !Add the second term to the idir component
           CALL cp_fm_set_all(fm_work1,0.0_dp,error=error)
           CALL cp_sm_fm_multiply(op_p_ao(iii)%matrix,fm_Rmd_mos(ii)%matrix,&
                fm_work1,ncol=homo,alpha=-1.0_dp,error=error)
           CALL cp_fm_scale_and_add(1.0_dp,rxp_psi0(ispin,idir)%matrix,&
                1.0_dp, fm_work1,error=error)

           CALL cp_fm_set_all(fm_work1,0.0_dp,error=error)
           CALL cp_sm_fm_multiply(op_p_ao(ii)%matrix,fm_Rmd_mos(iii)%matrix,&
                fm_work1,ncol=homo,alpha=-1.0_dp,error=error)
           CALL cp_fm_scale_and_add(1.0_dp,rxp_psi0(ispin,idir)%matrix,&
                -1.0_dp, fm_work1,error=error)

       END DO ! idir

       DO idir = 1,3
         CALL cp_fm_release (fm_Rmd_mos(idir)%matrix,error=error)
       END DO
       CALL cp_fm_release (fm_work1,error=error)

    END DO  ! ispin

    DEALLOCATE (first_sgf,last_sgf,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(vecbuf_c0,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO idir = 1,3
      DEALLOCATE(vecbuf_Rmdc0(idir)%array,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END DO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE nmr_operators

! *****************************************************************************
  FUNCTION fac_vecp(a,b,c) RESULT(factor)

    INTEGER                                  :: a, b, c
    REAL(dp)                                 :: factor

     factor = 0.0_dp

     IF ((b.eq.a+1 .OR. b.eq.a-2) .AND. (c.eq.b+1 .OR. c.eq.b-2)) THEN
       factor = 1.0_dp
     ELSEIF ((b.eq.a-1 .OR. b.eq.a+2) .AND. (c.eq.b-1 .OR. c.eq.b+2)) THEN
       factor = -1.0_dp
     END IF

  END FUNCTION fac_vecp

! *****************************************************************************
  FUNCTION ind_m2(ii,iii) RESULT(i)

    INTEGER                                  :: ii, iii, i

    INTEGER                                  :: l(3)

     i = 0
     l(1:3) = 0
     IF(ii==0) THEN
       l(iii) = 1
     ELSEIF(iii==0) THEN
       l(ii) = 1
     ELSEIF(ii==iii) THEN
       l(ii) = 2
       i = coset(l(1),l(2),l(3)) - 1
     ELSE
       l(ii) = 1
       l(iii) = 1
     ENDIF
     i = coset(l(1),l(2),l(3)) - 1
  END FUNCTION ind_m2

! *****************************************************************************
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
! *****************************************************************************
  SUBROUTINE set_vecp_rev(i1,i2,i3)

    INTEGER, INTENT(IN)                      :: i1, i2
    INTEGER, INTENT(OUT)                     :: i3

    IF((i1+i2)==3) THEN
      i3=3
    ELSEIF((i1+i2)==4) THEN
      i3=2
    ELSEIF((i1+i2)==5) THEN
      i3=1
    ELSE
    END IF

  END SUBROUTINE set_vecp_rev

END MODULE qs_linres_nmr_op
