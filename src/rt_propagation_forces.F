!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines needed for EMD
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

MODULE  rt_propagation_forces
  USE admm_types,                      ONLY: admm_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_control_types,                ONLY: dft_control_type,&
                                             rtp_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_col_block_sizes, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_deallocate_matrix, cp_dbcsr_distribution, &
       cp_dbcsr_get_block_p, cp_dbcsr_get_data_size, cp_dbcsr_get_data_type, &
       cp_dbcsr_get_num_blocks, cp_dbcsr_init, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_multiply, &
       cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_type, dbcsr_create_dist_r_unrot, dbcsr_distribution_obj, &
       dbcsr_distribution_release, dbcsr_type_no_symmetry
  USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_struct,                    ONLY: cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE input_constants,                 ONLY: use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: one,&
                                             zero
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: add_qs_force,&
                                             qs_force_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_overlap,                      ONLY: build_overlap_force
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_c_mat_force, &
            rt_admm_force

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_forces'


  CONTAINS



! *****************************************************************************
!> \brief calculates the three additional force contributions needed in EMD
!>        P_imag*C , P_imag*B*S^-1*S_der , P*S^-1*H*S_der
!>        driver routine
!> \param qs_env ...
!> \param error ...
!> \par History
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE calc_c_mat_force(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_c_mat_force', &
      routineP = moduleN//':'//routineN

    IF(qs_env%rtp%linear_scaling)THEN
       CALL calc_c_mat_force_ls(qs_env,error)
    ELSE 
       CALL calc_c_mat_force_fm(qs_env,error)
    END IF

  END SUBROUTINE calc_c_mat_force

! *****************************************************************************
!> \brief standard treatment for fm MO based calculations
!>        P_imag*C , P_imag*B*S^-1*S_der , P*S^-1*H*S_der
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calc_c_mat_force_fm(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_c_mat_force_fm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, im, ispin, nao, &
                                                natom, nmo, re, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of
    INTEGER, DIMENSION(:), POINTER           :: col_blk_size
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: C_mat, S_der, SinvB, SinvH
    TYPE(cp_dbcsr_type)                      :: db_mo_tmp1, db_mo_tmp2, &
                                                db_mos_im, db_mos_re
    TYPE(cp_dbcsr_type), POINTER             :: rho_im_sparse, tmp_dbcsr
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(rt_prop_type), POINTER              :: rtp

    failure=.FALSE.

    CALL timeset(routineN,handle)

    NULLIFY(rtp,particle_set,atomic_kind_set,mos)
    NULLIFY(tmp_dbcsr,rho_im_sparse)
    CALL get_qs_env(qs_env=qs_env,rtp=rtp,particle_set=particle_set,&
                 atomic_kind_set=atomic_kind_set,mos=mos,force=force,error=error)

    CALL get_rtp(rtp=rtp,C_mat=C_mat,S_der=S_der,&
                 SinvH=SinvH,SinvB=SinvB,mos_new=mos_new,error=error)

    ALLOCATE(tmp_dbcsr,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(rho_im_sparse,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(tmp_dbcsr,error=error)
    CALL cp_dbcsr_create(tmp_dbcsr,template=SinvH(1)%matrix, error=error)
    CALL cp_dbcsr_init(rho_im_sparse,error=error)
    CALL cp_dbcsr_create(rho_im_sparse,template=SinvH(1)%matrix, error=error)

    natom = SIZE(particle_set)
    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (kind_of(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind,kind_of=kind_of)

    DO ispin=1,SIZE(SinvH)
       re=2*ispin-1
       im=2*ispin
       alpha=mos(ispin)%mo_set%maxocc

       CALL get_mo_set(mos(ispin)%mo_set, nao=nao,nmo=nmo)
       
       NULLIFY(col_blk_size)
       CALL dbcsr_create_dist_r_unrot (dist, cp_dbcsr_distribution(SinvB(ispin)%matrix), nmo, col_blk_size)

       CALL cp_dbcsr_init(db_mos_re, error)
       CALL cp_dbcsr_create(db_mos_re, "D", dist, dbcsr_type_no_symmetry, &
            cp_dbcsr_row_block_sizes(SinvB(ispin)%matrix), col_blk_size,&
            0, 0, error=error)
       CALL cp_dbcsr_init(db_mos_im, error)
       CALL cp_dbcsr_create(db_mos_im,template=db_mos_re, error=error)
       CALL cp_dbcsr_init(db_mo_tmp1, error)
       CALL cp_dbcsr_create(db_mo_tmp1,template=db_mos_re, error=error)
       CALL cp_dbcsr_init(db_mo_tmp2, error)
       CALL cp_dbcsr_create(db_mo_tmp2,template=db_mos_re, error=error)

       CALL copy_fm_to_dbcsr(mos_new(im)%matrix,db_mos_im,error=error)
       CALL copy_fm_to_dbcsr(mos_new(re)%matrix,db_mos_re,error=error)

       CALL cp_dbcsr_multiply("N","N",alpha, SinvB(ispin)%matrix,db_mos_im,0.0_dp,db_mo_tmp1,error=error)       
       CALL cp_dbcsr_multiply("N","N",alpha, SinvH(ispin)%matrix,db_mos_re,1.0_dp,db_mo_tmp1,error=error)       
       CALL cp_dbcsr_multiply("N","N",-alpha,SinvB(ispin)%matrix,db_mos_re,0.0_dp,db_mo_tmp2,error=error)       
       CALL cp_dbcsr_multiply("N","N",alpha, SinvH(ispin)%matrix,db_mos_im,1.0_dp,db_mo_tmp2,error=error)
       CALL cp_dbcsr_multiply("N","T",1.0_dp,db_mo_tmp1,db_mos_re,0.0_dp,tmp_dbcsr,error=error)
       CALL cp_dbcsr_multiply("N","T",1.0_dp,db_mo_tmp2,db_mos_im,1.0_dp,tmp_dbcsr,error=error)
       
       CALL cp_dbcsr_multiply("N","T",alpha,db_mos_re,db_mos_im,0.0_dp,rho_im_sparse,error=error)
       CALL cp_dbcsr_multiply("N","T",-alpha,db_mos_im,db_mos_re,1.0_dp,rho_im_sparse,error=error)
              
       CALL compute_forces(force,tmp_dbcsr,S_der,rho_im_sparse,C_mat,kind_of,atom_of_kind,error) 

       CALL cp_dbcsr_release(db_mos_re,error)     
       CALL cp_dbcsr_release(db_mos_im,error)     
       CALL cp_dbcsr_release(db_mo_tmp1,error)     
       CALL cp_dbcsr_release(db_mo_tmp2,error)     
       
       DEALLOCATE(col_blk_size)
       CALL dbcsr_distribution_release(dist) 

    END DO     
 
    DO i=1,SIZE(force)
       force(i)%ehrenfest(:,:)=- force(i)%ehrenfest(:,:)
    END DO
 
    CALL cp_dbcsr_release(tmp_dbcsr,error)     
    CALL cp_dbcsr_release(rho_im_sparse,error)     
    DEALLOCATE (atom_of_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tmp_dbcsr,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(rho_im_sparse,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     
    CALL timestop(handle)


  END SUBROUTINE calc_c_mat_force_fm

! *****************************************************************************
!> \brief special treatment ofr linear scaling
!>        P_imag*C , P_imag*B*S^-1*S_der , P*S^-1*H*S_der
!> \param qs_env ...
!> \param error ...
!> \par History
!>      02.2014 switched to dbcsr matrices [Samuel Andermatt]
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE calc_c_mat_force_ls(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_c_mat_force_ls', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp, zero = 0.0_dp

    INTEGER                                  :: handle, i, im, ispin, natom, &
                                                re, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: C_mat, rho_new, S_der, SinvB, &
                                                SinvH
    TYPE(cp_dbcsr_type), POINTER             :: S_inv, S_minus_half, tmp, &
                                                tmp2, tmp3
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(rtp,particle_set,atomic_kind_set,dft_control)

    CALL get_qs_env(qs_env,&
                    rtp=rtp,&
                    particle_set=particle_set,&
                    atomic_kind_set=atomic_kind_set,&
                    force=force,&
                    dft_control=dft_control,&
                    error=error)

    rtp_control=>dft_control%rtp_control
    CALL get_rtp(rtp=rtp,C_mat=C_mat,S_der=S_der,S_inv=S_inv,&
                 SinvH=SinvH,SinvB=SinvB,error=error)

    natom = SIZE(particle_set)
    ALLOCATE (atom_of_kind(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (kind_of(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,atom_of_kind=atom_of_kind,kind_of=kind_of)


    NULLIFY(tmp)
    ALLOCATE(tmp,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(tmp,error=error)
    CALL cp_dbcsr_create(tmp,template=SinvB(1)%matrix,error=error)
    NULLIFY(tmp2)
    ALLOCATE(tmp2,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(tmp2,error=error)
    CALL cp_dbcsr_create(tmp2,template=SinvB(1)%matrix,error=error)
    NULLIFY(tmp3)
    ALLOCATE(tmp3,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_init(tmp3,error=error)
    CALL cp_dbcsr_create(tmp3,template=SinvB(1)%matrix,error=error)

    CALL get_rtp(rtp=rtp,rho_new=rho_new,S_minus_half=S_minus_half,error=error)

    DO ispin=1,SIZE(SinvH)
       re=2*ispin-1
       im=2*ispin
       IF(rtp_control%orthonormal) THEN
          CALL cp_dbcsr_multiply("N","N",one,S_minus_half,SinvH(ispin)%matrix,zero,tmp2,&
               filter_eps=rtp%filter_eps,error=error)
          CALL cp_dbcsr_multiply("N","N",one,tmp2,rho_new(re)%matrix,zero,tmp3,&
               filter_eps=rtp%filter_eps,error=error)
          CALL cp_dbcsr_multiply("N","N",one,tmp3,S_minus_half,zero,tmp,&
               filter_eps=rtp%filter_eps,error=error)
          CALL cp_dbcsr_multiply("N","N",one,S_minus_half,SinvB(ispin)%matrix,zero,tmp2,&
               filter_eps=rtp%filter_eps,error=error)
          CALL cp_dbcsr_multiply("N","N",one,tmp2,rho_new(im)%matrix,zero,tmp3,&
               filter_eps=rtp%filter_eps,error=error)
          CALL cp_dbcsr_multiply("N","N",one,tmp3,S_minus_half,one,tmp,&
               filter_eps=rtp%filter_eps,error=error)
       ELSE
          CALL cp_dbcsr_multiply("N","N",one,SinvH(ispin)%matrix,rho_new(re)%matrix,zero,tmp,&
               filter_eps=rtp%filter_eps,error=error)
          CALL cp_dbcsr_multiply("N","N",-one,SinvB(ispin)%matrix,rho_new(im)%matrix,one,tmp,&
               filter_eps=rtp%filter_eps,error=error)
       ENDIF

       CALL compute_forces(force,tmp,S_der,rho_new(im)%matrix,C_mat,kind_of,atom_of_kind,error)

    END DO

    ! recall QS forces, at this point have the other sign.
    DO i=1,SIZE(force)
       force(i)%ehrenfest(:,:)=- force(i)%ehrenfest(:,:)
    END DO

    CALL cp_dbcsr_deallocate_matrix(tmp,error=error)
    CALL cp_dbcsr_deallocate_matrix(tmp2,error=error)
    CALL cp_dbcsr_deallocate_matrix(tmp3,error=error)

    DEALLOCATE (atom_of_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param force ...
!> \param tmp ...
!> \param S_der ...
!> \param rho_im ...
!> \param C_mat ...
!> \param kind_of ...
!> \param atom_of_kind ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE compute_forces(force,tmp,S_der,rho_im,C_mat,kind_of,atom_of_kind,error)
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(cp_dbcsr_type), POINTER             :: tmp
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: S_der
    TYPE(cp_dbcsr_type), POINTER             :: rho_im
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: C_mat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of, atom_of_kind
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: col_atom, i, ikind, &
                                                kind_atom, row_atom
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:), POINTER          :: block_values, block_values2
    TYPE(cp_dbcsr_iterator)                  :: iter

    DO i=1,3
       !Calculate the sum over the hadmard product
       !S_der part

       CALL cp_dbcsr_iterator_start(iter, tmp)
       DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
          CALL cp_dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
          CALL cp_dbcsr_get_block_p(S_der(i)%matrix, row_atom, col_atom, block_values2, found=found)
          IF(found) THEN
             ikind=kind_of(col_atom)
             kind_atom=atom_of_kind(col_atom)
             !The block_values are in a vector format,
             ! so the dot_product is the sum over all elements of the hamand product, that I need
             force(ikind)%ehrenfest(i,kind_atom)=force(ikind)%ehrenfest(i,kind_atom)+&
             2.0_dp*DOT_PRODUCT(block_values,block_values2)
          ENDIF
       END DO
       CALL cp_dbcsr_iterator_stop (iter)

       !C_mat part

       CALL cp_dbcsr_iterator_start(iter, rho_im)
       DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
          CALL cp_dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
          CALL cp_dbcsr_get_block_p(C_mat(i)%matrix,row_atom,col_atom,block_values2,found=found)
          IF(found) THEN
             ikind=kind_of(col_atom)
             kind_atom=atom_of_kind(col_atom)
             !The block_values are in a vector format, so the dot_product is
             ! the sum over all elements of the hamand product, that I need
             force(ikind)%ehrenfest(i,kind_atom)=force(ikind)%ehrenfest(i,kind_atom)+&
             2.0_dp*DOT_PRODUCT(block_values,block_values2)
          ENDIF
       END DO
       CALL cp_dbcsr_iterator_stop (iter)
    END DO

 END SUBROUTINE compute_forces

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rt_admm_force(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_admm_force', &
      routineP = moduleN//':'//routineN

    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: KS_aux_im, KS_aux_re, &
                                                matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_admm
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL get_qs_env(qs_env,&
                    admm_env=admm_env,&
                    rtp=rtp,&
                    matrix_ks_aux_fit=KS_aux_re,&
                    matrix_ks_aux_fit_im=KS_aux_im,&
                    matrix_s_aux_fit=matrix_s_aux_fit,&
                    matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb,&
                    error=error)

    CALL get_rtp(rtp=rtp,mos_new=mos,admm_mos=mos_admm,error=error)

    ! currently only none option
    CALL rt_admm_forces_none(qs_env, admm_env, KS_aux_re, KS_aux_im, &
                             matrix_s_aux_fit, matrix_s_aux_fit_vs_orb , mos_admm,mos,error)

  END SUBROUTINE rt_admm_force

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param admm_env ...
!> \param KS_aux_re ...
!> \param KS_aux_im ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_aux_fit_vs_orb ...
!> \param mos_admm ...
!> \param mos ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rt_admm_forces_none(qs_env,admm_env,KS_aux_re,KS_aux_im,matrix_s_aux_fit, matrix_s_aux_fit_vs_orb,mos_admm,mos,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: KS_aux_re, KS_aux_im, &
                                                matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_admm, mos
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_admm_forces_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: im, ispin, istat, nao, natom, &
                                                naux, nmo, re
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: admm_force
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_type), POINTER             :: matrix_w_q, matrix_w_s
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: tmp_aux_aux, tmp_aux_mo, &
                                                tmp_aux_mo1, tmp_aux_nao
    TYPE(cp_fm_struct_type), POINTER         :: mstruct
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_aux_fit_asymm, &
                                                sab_aux_fit_vs_orb
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    NULLIFY(sab_aux_fit_asymm, sab_aux_fit_vs_orb, ks_env)
    failure = .FALSE.
!   CALL cp_fm_create(tmp_aux_aux,admm_env%fm_struct_tmp,name="fm matrix",error=error)

    CALL get_qs_env(qs_env,&
                    sab_aux_fit_asymm=sab_aux_fit_asymm,&
                    sab_aux_fit_vs_orb=sab_aux_fit_vs_orb,&
                    ks_env=ks_env,&
                    error=error)

    ALLOCATE(matrix_w_s)
    CALL cp_dbcsr_init (matrix_w_s, error)
    CALL cp_dbcsr_create(matrix_w_s, 'W MATRIX AUX S', &
         cp_dbcsr_distribution(matrix_s_aux_fit(1)%matrix), dbcsr_type_no_symmetry, &
         cp_dbcsr_row_block_sizes(matrix_s_aux_fit(1)%matrix),&
         cp_dbcsr_col_block_sizes(matrix_s_aux_fit(1)%matrix), &
         cp_dbcsr_get_num_blocks(matrix_s_aux_fit(1)%matrix), &
         cp_dbcsr_get_data_size(matrix_s_aux_fit(1)%matrix),&
         cp_dbcsr_get_data_type(matrix_s_aux_fit(1)%matrix), &
         error=error)
    CALL cp_dbcsr_alloc_block_from_nbl(matrix_w_s,sab_aux_fit_asymm,error=error)

    ALLOCATE(matrix_w_q)
    CALL cp_dbcsr_init(matrix_w_q, error=error)
    CALL cp_dbcsr_copy(matrix_w_q,matrix_s_aux_fit_vs_orb(1)%matrix,&
                    "W MATRIX AUX Q",error=error)

    DO ispin=1,SIZE(KS_aux_re)
       re=2*ispin-1; im=2*ispin
       naux=admm_env%nao_aux_fit; nmo=admm_env%nmo(ispin); nao=admm_env%nao_orb
       
       ALLOCATE(tmp_aux_aux(2),tmp_aux_nao(2),tmp_aux_mo(2),tmp_aux_mo1(2))
       CALL cp_fm_create(tmp_aux_aux(1)%matrix,admm_env%work_aux_aux%matrix_struct,name="taa",error=error)
       CALL cp_fm_create(tmp_aux_aux(2)%matrix,admm_env%work_aux_aux%matrix_struct,name="taa",error=error)
       CALL cp_fm_create(tmp_aux_nao(1)%matrix,admm_env%work_aux_orb%matrix_struct,name="tao",error=error)
       CALL cp_fm_create(tmp_aux_nao(2)%matrix,admm_env%work_aux_orb%matrix_struct,name="tao",error=error)
       mstruct => admm_env%work_aux_nmo(ispin)%matrix%matrix_struct
       CALL cp_fm_create(tmp_aux_mo(1)%matrix,mstruct,name="tam",error=error)
       CALL cp_fm_create(tmp_aux_mo(2)%matrix,mstruct,name="tam",error=error)
       CALL cp_fm_create(tmp_aux_mo1(1)%matrix,mstruct,name="tam",error=error)
       CALL cp_fm_create(tmp_aux_mo1(2)%matrix,mstruct,name="tam",error=error)

! First calculate H=KS_aux*C~, real part ends on work_aux_aux2, imaginary part ends at work_aux_aux3
       CALL cp_dbcsr_sm_fm_multiply(KS_aux_re(ispin)%matrix,mos_admm(re)%matrix,tmp_aux_mo(re)%matrix,nmo,4.0_dp,0.0_dp,error)
       CALL cp_dbcsr_sm_fm_multiply(KS_aux_re(ispin)%matrix,mos_admm(im)%matrix,tmp_aux_mo(im)%matrix,nmo,4.0_dp,0.0_dp,error)
       CALL cp_dbcsr_sm_fm_multiply(KS_aux_im(ispin)%matrix,mos_admm(im)%matrix,tmp_aux_mo(re)%matrix,nmo,-4.0_dp,1.0_dp,error)
       CALL cp_dbcsr_sm_fm_multiply(KS_aux_im(ispin)%matrix,mos_admm(re)%matrix,tmp_aux_mo(im)%matrix,nmo,4.0_dp,1.0_dp,error)
 
! Next step compute S-1*H
       CALL cp_gemm('N','N',naux,nmo,naux,1.0_dp,admm_env%S_inv,tmp_aux_mo(re)%matrix,0.0_dp,tmp_aux_mo1(re)%matrix,error)
       CALL cp_gemm('N','N',naux,nmo,naux,1.0_dp,admm_env%S_inv,tmp_aux_mo(im)%matrix,0.0_dp,tmp_aux_mo1(im)%matrix,error)

! Here we go on with Ws=S-1*H * C^H (take care of sign of the imaginary part!!!)

       CALL cp_gemm("N","T",naux,nao,nmo,-1.0_dp,tmp_aux_mo1(re)%matrix, mos(re)%matrix, 0.0_dp,&
                       tmp_aux_nao(re)%matrix, error)
       CALL cp_gemm("N","T",naux,nao,nmo,-1.0_dp,tmp_aux_mo1(im)%matrix, mos(im)%matrix, 1.0_dp,&
                       tmp_aux_nao(re)%matrix, error)
       CALL cp_gemm("N","T",naux,nao,nmo,1.0_dp,tmp_aux_mo1(re)%matrix, mos(im)%matrix, 0.0_dp,&
                       tmp_aux_nao(im)%matrix, error)
       CALL cp_gemm("N","T",naux,nao,nmo,-1.0_dp,tmp_aux_mo1(im)%matrix, mos(re)%matrix, 1.0_dp,&
                       tmp_aux_nao(im)%matrix, error)

! Let's do the final bit  Wq=S-1*H * C^H * A^T
       CALL cp_gemm('N','T',naux,naux,nao,-1.0_dp,tmp_aux_nao(re)%matrix,admm_env%A,0.0_dp,tmp_aux_aux(re)%matrix,error)     
       CALL cp_gemm('N','T',naux,naux,nao,-1.0_dp,tmp_aux_nao(im)%matrix,admm_env%A,0.0_dp,tmp_aux_aux(im)%matrix,error)     

       ! *** copy to sparse matrix
       CALL copy_fm_to_dbcsr(tmp_aux_nao(re)%matrix, matrix_w_q,keep_sparsity=.TRUE.,&
            error=error)


       ! *** copy to sparse matrix
       CALL copy_fm_to_dbcsr(tmp_aux_aux(re)%matrix, matrix_w_s,keep_sparsity=.TRUE.,&
            error=error)

! *** This can be done in one call w_total = w_alpha + w_beta
       ! allocate force vector
       CALL get_qs_env(qs_env=qs_env,natom=natom,error=error)
       ALLOCATE(admm_force(3,natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       admm_force = 0.0_dp
       CALL build_overlap_force(ks_env, admm_force,&
            basis_set_id_a=use_aux_fit_basis_set, basis_set_id_b=use_aux_fit_basis_set, &
            sab_nl=sab_aux_fit_asymm, matrix_p=matrix_w_s, error=error)
       CALL build_overlap_force(ks_env, admm_force,&
            basis_set_id_a=use_aux_fit_basis_set, basis_set_id_b=use_orb_basis_set, &
            sab_nl=sab_aux_fit_vs_orb, matrix_p=matrix_w_q, error=error)
       ! add forces
       CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,&
                       force=force,error=error)
       CALL add_qs_force(admm_force, force, "overlap_admm", atomic_kind_set, error)
       DEALLOCATE(admm_force,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       
       ! *** Deallocated weighted density matrices
       CALL cp_dbcsr_deallocate_matrix(matrix_w_s,error)
       CALL cp_dbcsr_deallocate_matrix(matrix_w_q,error)
       CALL cp_fm_vect_dealloc(tmp_aux_aux,error)
       CALL cp_fm_vect_dealloc(tmp_aux_nao,error)
       CALL cp_fm_vect_dealloc(tmp_aux_mo,error)
       CALL cp_fm_vect_dealloc(tmp_aux_mo1,error)
     END DO
    
  END SUBROUTINE rt_admm_forces_none

END MODULE rt_propagation_forces
