!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines needed for EMD
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

MODULE  rt_propagation_utils
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_binary_read,&
                                             cp_dbcsr_checksum,&
                                             cp_dbcsr_copy,&
                                             cp_dbcsr_distribution,&
                                             cp_dbcsr_p_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: use_restart_wfn,&
                                             use_rt_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE mathconstants,                   ONLY: zero
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_mo_io,                        ONLY: read_mo_set,&
                                             read_rt_mos_from_restart
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: mo_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE rt_propagation_methods,          ONLY: calc_update_rho_sparse
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC:: get_restart_wfn

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_utils'


  CONTAINS

! *****************************************************************************
!> \brief reads the restart file. At the moment only SCF (means only real)
!> \param qs_env ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************


  SUBROUTINE get_restart_wfn(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_restart_wfn', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name
    INTEGER                                  :: i, id_nr, im, ispin, ncol, &
                                                nspin, re, unit_nr
    REAL(KIND=dp)                            :: alpha, cs_pos
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_rmpv, rho_new, rho_old
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new, mos_old
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho_struct
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(section_vals_type), POINTER         :: dft_section, input

    NULLIFY(atomic_kind_set,qs_kind_set,mo_array,particle_set,rho_struct,para_env)

    CALL get_qs_env(qs_env,&
                    qs_kind_set=qs_kind_set,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    mos=mo_array,&
                    input=input,&
                    rtp=rtp,&
                    dft_control=dft_control,&
                    rho=rho_struct,&
                    para_env=para_env,&
                    error=error)
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    id_nr=0
    nspin=SIZE(mo_array)
    p_rmpv => rho_struct%rho_ao
    dft_section =>  section_vals_get_subs_vals(input,"DFT",error=error)
    SELECT CASE(dft_control%rtp_control%initial_wfn)
    CASE(use_restart_wfn)
       CALL read_mo_set(mo_array,atomic_kind_set,qs_kind_set,particle_set,para_env,&
            id_nr=id_nr,multiplicity=dft_control%multiplicity,dft_section=dft_section,&
            error=error)
       DO ispin=1,nspin
          CALL calculate_density_matrix(mo_array(ispin)%mo_set, p_rmpv(ispin)%matrix,error=error)
       ENDDO
       IF(rtp%linear_scaling) THEN
          CALL get_rtp(rtp=rtp,rho_old=rho_old,rho_new=rho_new,error=error)
          DO ispin=1,nspin
             re=2*ispin-1
             im=2*ispin
             CALL cp_fm_get_info(mo_array(ispin)%mo_set%mo_coeff,ncol_global=ncol,error=error)
             alpha=1.0_dp
             IF(SIZE(mo_array)==1) alpha=2*alpha
             CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho_old(re)%matrix,&
                  matrix_v=mo_array(ispin)%mo_set%mo_coeff,matrix_g=mo_array(ispin)%mo_set%mo_coeff,ncol=ncol,&
                  keep_sparsity=.FALSE.,alpha=alpha,error=error)
          END DO
          DO i=1,nspin
             CALL cp_dbcsr_copy(rho_new(i)%matrix,rho_old(i)%matrix,error=error)
          ENDDO
          CALL calc_update_rho_sparse(qs_env,error)
       ELSE
          CALL get_rtp(rtp=rtp,mos_old=mos_old,error=error)
          DO i=1,SIZE(qs_env%mos)
             CALL cp_fm_to_fm(mo_array(i)%mo_set%mo_coeff,mos_old(2*i-1)%matrix,error)
             CALL cp_fm_set_all(mos_old(2*i)%matrix,zero,zero,error)
          END DO
       ENDIF
    CASE(use_rt_restart)
       IF(rtp%linear_scaling) THEN
          CALL get_rtp(rtp=rtp,rho_old=rho_old,rho_new=rho_new,error=error)
          DO ispin=1,nspin
             re=2*ispin-1
             im=2*ispin
             WRITE(file_name,'(A,I0,A)') "LS_DM_SPIN_RE",ispin,"_RESTART.dm"
             CALL cp_dbcsr_binary_read(file_name, distribution=cp_dbcsr_distribution(rho_old(re)%matrix), &
                                       matrix_new=rho_old(re)%matrix, error=error)
             cs_pos = cp_dbcsr_checksum (rho_old(re)%matrix, pos=.TRUE., error=error)
             IF (unit_nr>0) THEN
                WRITE(unit_nr,'(T2,A,E20.8)') "Read restart DM "//TRIM(file_name)//" with checksum: ",cs_pos
             ENDIF
             WRITE(file_name,'(A,I0,A)') "LS_DM_SPIN_IM",ispin,"_RESTART.dm"
             CALL cp_dbcsr_binary_read(file_name, distribution=cp_dbcsr_distribution(rho_old(im)%matrix), &
                                       matrix_new=rho_old(im)%matrix, error=error)
             cs_pos = cp_dbcsr_checksum (rho_old(im)%matrix, pos=.TRUE., error=error)
             IF (unit_nr>0) THEN
                WRITE(unit_nr,'(T2,A,E20.8)') "Read restart DM "//TRIM(file_name)//" with checksum: ",cs_pos
             ENDIF
          ENDDO
          DO i=1,SIZE(rho_new)
             CALL cp_dbcsr_copy(rho_new(i)%matrix,rho_old(i)%matrix,error=error)
          ENDDO
          CALL calc_update_rho_sparse(qs_env,error)
       ELSE
          CALL get_rtp(rtp=rtp,mos_old=mos_old,mos_new=mos_new,error=error)
          CALL read_rt_mos_from_restart(mo_array,mos_old,atomic_kind_set,qs_kind_set,particle_set,para_env,&
               id_nr,dft_control%multiplicity,dft_section, error)
          DO ispin=1,nspin
             CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
                  p_rmpv(ispin)%matrix,error=error)
          ENDDO
       ENDIF
    END SELECT

  END SUBROUTINE get_restart_wfn

END MODULE rt_propagation_utils
