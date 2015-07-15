!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines needed for EMD
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

MODULE  rt_propagation_utils
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             rtp_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_binary_read, cp_dbcsr_checksum, cp_dbcsr_copy, &
       cp_dbcsr_copy_into_existing, cp_dbcsr_create, &
       cp_dbcsr_deallocate_matrix, cp_dbcsr_deallocate_matrix_set, &
       cp_dbcsr_desymmetrize, cp_dbcsr_distribution, cp_dbcsr_filter, &
       cp_dbcsr_init, cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_p_type, &
       cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: use_orb_basis_set,&
                                             use_restart_wfn,&
                                             use_rt_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE mathconstants,                   ONLY: zero
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_dftb_matrices,                ONLY: build_dftb_overlap
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_ks_types,                     ONLY: qs_ks_did_change,&
                                             qs_ks_env_type
  USE qs_mo_io,                        ONLY: read_mo_set,&
                                             read_rt_mos_from_restart
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_overlap,                      ONLY: build_overlap_matrix
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_set,&
                                             qs_rho_type
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC:: get_restart_wfn,&
           calc_S_derivs,&
           calc_update_rho,&
           calc_update_rho_sparse,&
           calculate_P_imaginary

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_utils'


  CONTAINS

! *****************************************************************************
!> \brief Calculates dS/dR respectily the velocity weighted derivatves
!>        only needed for ehrenfest MD.
!>
!> \param qs_env the qs environment
!> \param error ...
!>
!> \par History
!>      02.2009 created [Manuel Guidon]
!>      02.2014 switched to dbcsr matrices [Samuel Andermatt]
!> \author Florian Schiffmann
! *****************************************************************************
  SUBROUTINE calc_S_derivs(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_S_derivs', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: col_atom, handle, i, j, m, &
                                                maxder, n, nder, &
                                                neighbor_list_id, row_atom
    INTEGER, DIMENSION(6, 2)                 :: c_map_mat
    LOGICAL                                  :: failure, return_s_derivatives
    REAL(dp), DIMENSION(:), POINTER          :: block_values
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: C_mat, S_der, s_derivs
    TYPE(cp_dbcsr_type), POINTER             :: B_mat, tmp_mat, tmp_mat2
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(rt_prop_type), POINTER              :: rtp

    failure = .FALSE.
    CALL timeset(routineN,handle)

    return_s_derivatives=.TRUE.

    NULLIFY (atomic_kind_set)
    NULLIFY (force)
    NULLIFY (particle_set)
    NULLIFY (rtp)
    NULLIFY (s_derivs)
    NULLIFY (dft_control)
    NULLIFY (ks_env)

    CALL get_qs_env(qs_env=qs_env,&
                    rtp=rtp,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    neighbor_list_id=neighbor_list_id,&
                    sab_orb=sab_orb,&
                    force=force,&
                    dft_control=dft_control,&
                    ks_env=ks_env,&
                    error=error)

    CALL get_rtp(rtp=rtp,B_mat=B_mat,C_mat=C_mat,S_der=S_der,error=error)

    nder = 2
    maxder = ncoset(nder)

    NULLIFY(tmp_mat)
    ALLOCATE(tmp_mat)
    CALL cp_dbcsr_init(tmp_mat,error=error)
    CALL cp_dbcsr_create(tmp_mat,template=S_der(1)%matrix,matrix_type="N",error=error)

    IF(rtp%iter<2) THEN
       ! calculate the overlap derivative matrices
       IF(dft_control%qs_control%dftb)THEN
          CALL build_dftb_overlap(qs_env,nder,s_derivs,error)
       ELSE
          CALL build_overlap_matrix(ks_env,nderivative=nder,matrix_s=s_derivs,&
               basis_set_id_a=use_orb_basis_set,&
               basis_set_id_b=use_orb_basis_set,sab_nl=sab_orb,error=error)
       END IF        

       NULLIFY(tmp_mat2)
       ALLOCATE(tmp_mat2)
       CALL cp_dbcsr_init(tmp_mat2,error=error)
       CALL cp_dbcsr_create(tmp_mat2,template=S_der(1)%matrix,matrix_type="S",error=error)
       DO m=1,9
          CALL cp_dbcsr_copy(tmp_mat2,s_derivs(m+1)%matrix,error=error)
          CALL cp_dbcsr_desymmetrize(tmp_mat2,S_der(m)%matrix,error=error)
          CALL cp_dbcsr_scale(S_der(m)%matrix,-one,error=error) 
          CALL cp_dbcsr_filter(S_der(m)%matrix,rtp%filter_eps,error=error)
          !The diagonal should be zero
          CALL cp_dbcsr_iterator_start(iter, S_der(m)%matrix)
          DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
             CALL cp_dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
             IF(row_atom==col_atom) block_values=0.0_dp
          END DO
          CALL cp_dbcsr_iterator_stop (iter) 
       END DO
       CALL cp_dbcsr_deallocate_matrix_set(s_derivs,error=error)
       CALL cp_dbcsr_deallocate_matrix(tmp_mat2,error=error)
    END IF

    !calculate scalar product v(Rb)*<alpha|d/dRb beta> (B_mat), and store the first derivatives

     CALL cp_dbcsr_set(B_mat,zero,error=error)
     DO m=1,3
       CALL cp_dbcsr_copy(tmp_mat,S_der(m)%matrix,error=error)
       CALL cp_dbcsr_iterator_start(iter, tmp_mat)
       DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
          CALL cp_dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
          IF(row_atom==col_atom) block_values=0.0_dp
          block_values=block_values*particle_set(col_atom)%v(m)
       END DO
       CALL cp_dbcsr_iterator_stop (iter) 
       CALL cp_dbcsr_add(B_mat,tmp_mat,one,one,error=error)
    END DO
    CALL cp_dbcsr_filter(B_mat,rtp%filter_eps,error=error)
    !calculate C matrix: v(Rb)*<d/dRa alpha| d/dRb beta>

    c_map_mat=0
    n=0
    DO j=1,3
       DO m=j,3
          n=n+1
          c_map_mat(n,1)=j
          IF(m==j)CYCLE
          c_map_mat(n,2)=m
       END DO
    END DO


    DO i=1,3
       CALL cp_dbcsr_set(C_mat(i)%matrix,zero,error=error)
    END DO
    DO m=1,6
       CALL cp_dbcsr_copy(tmp_mat,S_der(m+3)%matrix,error=error)
       DO j=1,2
          IF(c_map_mat(m,j)==0)CYCLE
          CALL cp_dbcsr_add(C_mat(c_map_mat(m,j))%matrix,tmp_mat,one,one,error)
       END DO
    END DO

    DO m=1,3
       CALL cp_dbcsr_iterator_start(iter, C_mat(m)%matrix)
       DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
          CALL cp_dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
          block_values=block_values*particle_set(row_atom)%v(m)
       END DO
       CALL cp_dbcsr_iterator_stop (iter)
       CALL cp_dbcsr_filter(C_mat(m)%matrix,rtp%filter_eps,error=error)
    END DO


    CALL cp_dbcsr_deallocate_matrix(tmp_mat,error=error)
    CALL timestop(handle)
  END SUBROUTINE

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

    CHARACTER(LEN=default_path_length)       :: file_name, project_name
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
    CALL qs_rho_get(rho_struct, rho_ao=p_rmpv, error=error)
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
          project_name = logger%iter_info%project_name
          DO ispin=1,nspin
             re=2*ispin-1
             im=2*ispin
             WRITE(file_name,'(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_RE",ispin,"_RESTART.dm"
             CALL cp_dbcsr_binary_read(file_name, distribution=cp_dbcsr_distribution(rho_old(re)%matrix), &
                                       matrix_new=rho_old(re)%matrix, error=error)
             cs_pos = cp_dbcsr_checksum (rho_old(re)%matrix, pos=.TRUE., error=error)
             IF (unit_nr>0) THEN
                WRITE(unit_nr,'(T2,A,E20.8)') "Read restart DM "//TRIM(file_name)//" with checksum: ",cs_pos
             ENDIF
             WRITE(file_name,'(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_IM",ispin,"_RESTART.dm"
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
  
! *****************************************************************************
!> \brief calculates the density from the complex MOs and passes the density to
!>        qs_env.
!> \param qs_env ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE calc_update_rho(qs_env,error)


    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_update_rho', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, im, ncol, re
    REAL(KIND=dp)                            :: alpha
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao, rho_ao_im
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL timeset(routineN,handle)

    NULLIFY (rho,ks_env,mos,rtp)
    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    rho=rho,&
                    rtp=rtp,&
                    error=error)
    CALL get_rtp(rtp=rtp,mos_new=mos,error=error)
    CALL qs_rho_get(rho_struct=rho,rho_ao=rho_ao,error=error)
    DO i=1,SIZE(mos)/2
       re=2*i-1 ; im =2*i
       alpha=3*one-REAL(SIZE(mos)/2,dp)
       CALL cp_dbcsr_set(rho_ao(i)%matrix,zero,error=error)
       CALL cp_fm_get_info(mos(re)%matrix,ncol_global=ncol,error=error)
       CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho_ao(i)%matrix,&
                               matrix_v=mos(re)%matrix,&
                               ncol=ncol,&
                               alpha=alpha,error=error)
       ! It is actually complex conjugate but i*i=-1 therfore it must be added
       CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho_ao(i)%matrix,&
                        matrix_v=mos(im)%matrix,&
                        ncol=ncol,&
                        alpha=alpha,error=error)
    END DO

    CALL qs_rho_update_rho(rho, qs_env, error=error)

    IF(rtp%do_hfx)THEN
       CALL qs_rho_get(rho_struct=rho,rho_ao_im=rho_ao_im,error=error)
       CALL calculate_P_imaginary(rtp, rho_ao_im, error)
       CALL qs_rho_set(rho, rho_ao_im=rho_ao_im, error=error)
    END IF

    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)

    CALL timestop(handle)

  END SUBROUTINE calc_update_rho


! *****************************************************************************
!> \brief Copies the density matrix back into the qs_env%rho%rho_ao
!> \param qs_env ...
!> \param error ...
!> \author Samuel Andermatt (3.14)
! *****************************************************************************

  SUBROUTINE calc_update_rho_sparse(qs_env,error)


    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_update_rho_sparse', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, im, ispin, re
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao, rho_ao_im, rho_new
    TYPE(cp_dbcsr_type), POINTER             :: S_minus_half, tmp, tmp2
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control

    NULLIFY(rho,ks_env,rtp,tmp,tmp2,dft_control)
    CALL timeset(routineN,handle)
    CALL get_qs_env(qs_env,&
                   ks_env=ks_env,&
                   rho=rho,&
                   rtp=rtp,&
                   dft_control=dft_control,&
                   error=error)
    rtp_control=>dft_control%rtp_control
    CALL get_rtp(rtp=rtp,rho_new=rho_new,error=error)
    CALL qs_rho_get(rho_struct=rho,rho_ao=rho_ao,error=error)
    IF(rtp%do_hfx) CALL qs_rho_get(rho_struct=rho,rho_ao_im=rho_ao_im,error=error)
    IF(rtp_control%orthonormal) THEN
       CALL get_rtp(rtp=rtp,S_minus_half=S_minus_half,error=error)
       ALLOCATE(tmp)
       CALL cp_dbcsr_init(tmp,error=error)
       CALL cp_dbcsr_create(tmp,template=rho_new(1)%matrix,error=error)
       ALLOCATE(tmp2)
       CALL cp_dbcsr_init(tmp2,error=error)
       CALL cp_dbcsr_create(tmp2,template=rho_new(1)%matrix,error=error)
       DO ispin=1,SIZE(rho_ao)
          re = 2*ispin-1
          CALL cp_dbcsr_multiply("N","N",one,S_minus_half,rho_new(re)%matrix,zero,tmp,filter_eps=rtp%filter_eps,error=error)
          CALL cp_dbcsr_multiply("N","N",one,tmp,S_minus_half,zero,tmp2,filter_eps=rtp%filter_eps,error=error)
          CALL cp_dbcsr_set(rho_ao(ispin)%matrix,zero,error=error)
          CALL cp_dbcsr_copy_into_existing(rho_ao(ispin)%matrix,tmp2,error=error)
       END DO
       IF(rtp%do_hfx) THEN
          DO ispin=1,SIZE(rho_ao_im)
             im = 2*ispin
             CALL cp_dbcsr_multiply("N","N",one,S_minus_half,rho_new(im)%matrix,zero,tmp,filter_eps=rtp%filter_eps,error=error)
             CALL cp_dbcsr_multiply("N","N",one,tmp,S_minus_half,zero,tmp2,filter_eps=rtp%filter_eps,error=error)
             CALL cp_dbcsr_set(rho_ao_im(ispin)%matrix,zero,error=error)
             CALL cp_dbcsr_copy_into_existing(rho_ao_im(ispin)%matrix,tmp2,error=error)
          END DO
       ENDIF
       CALL cp_dbcsr_deallocate_matrix(tmp,error=error)
       CALL cp_dbcsr_deallocate_matrix(tmp2,error=error)
    ELSE
       DO ispin=1,SIZE(rho_ao)
          CALL cp_dbcsr_set(rho_ao(ispin)%matrix,zero,error=error)
          CALL cp_dbcsr_copy_into_existing(rho_ao(ispin)%matrix,rho_new(ispin*2-1)%matrix,error=error)
          IF(rtp%do_hfx) CALL cp_dbcsr_copy_into_existing(rho_ao_im(ispin)%matrix,rho_new(ispin*2)%matrix,error=error)
       END DO
    ENDIF

    CALL qs_rho_update_rho(rho, qs_env, error=error)
    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)


    CALL timestop(handle)

  END SUBROUTINE calc_update_rho_sparse
  
! *****************************************************************************
!> \brief ...
!> \param rtp ...
!> \param matrix_p_im ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_P_imaginary(rtp,matrix_p_im,error)
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p_im
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_P_imaginary', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, im, ncol, re
    REAL(KIND=dp)                            :: alpha
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos

    CALL get_rtp(rtp=rtp,mos_new=mos,error=error)

    DO i=1,SIZE(mos)/2
       re=2*i-1 ; im =2*i
       alpha=3.0_dp-REAL(SIZE(matrix_p_im),dp)
       CALL cp_dbcsr_set(matrix_p_im(i)%matrix,0.0_dp,error=error)
       CALL cp_fm_get_info(mos(re)%matrix,ncol_global=ncol,error=error)
       CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=matrix_p_im(i)%matrix,&
                               matrix_v=mos(im)%matrix,&
                               matrix_g=mos(re)%matrix,&
                               ncol=ncol,&
                               alpha=alpha,error=error)
       ! It is actually complex conjugate not only transposed
       alpha=-alpha
       CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=matrix_p_im(i)%matrix,&
                        matrix_v=mos(re)%matrix,&
                        matrix_g=mos(im)%matrix,&
                        ncol=ncol,&
                        alpha=alpha,error=error)
    END DO
  
  END SUBROUTINE calculate_P_imaginary

END MODULE rt_propagation_utils
