!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines needed for EMD
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

MODULE  rt_propagation_utils


  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       array_i1d_obj, array_release, cp_dbcsr_add, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix, &
       cp_dbcsr_deallocate_matrix_set, cp_dbcsr_desymmetrize, &
       cp_dbcsr_filter, cp_dbcsr_get_block, cp_dbcsr_init, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_multiply, &
       cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_type, &
       dbcsr_create_dist_r_unrot, dbcsr_distribution_obj, &
       dbcsr_distribution_release, dbcsr_type_no_symmetry
  USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                             cp_dbcsr_plus_fm_fm_t
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE input_constants,                 ONLY: use_orb_basis_set,&
                                             use_restart_wfn,&
                                             use_rt_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants
  USE orbital_pointers,                ONLY: ncoset
  USE particle_types,                  ONLY: particle_type
  USE qs_dftb_matrices,                ONLY: build_dftb_overlap
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type,&
                                             read_mo_set,&
                                             read_rt_mos_from_restart
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_overlap,                      ONLY: build_overlap_matrix
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC:: calc_S_derivs,&
           get_restart_wfn,&
           calculate_dp_imag,&
           calc_c_mat_force

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_utils'


  CONTAINS

! *****************************************************************************
!> \brief Calculates dS/dR respectily the velocity weighted derivatves
!>        only needed for ehrenfest MD.
!>
!> \param qs_env the qs environment
!> \param error
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

    INTEGER                                  :: col_atom, handle, i, j, m, &
                                                maxder, n, nder, &
                                                neighbor_list_id, row_atom
    INTEGER, DIMENSION(6, 2)                 :: c_map_mat
    LOGICAL                                  :: failure, return_s_derivatives
    REAL(dp), DIMENSION(:), POINTER          :: block_values
    REAL(KIND=dp)                            :: one = 1.0_dp, zero = 0.0_dp
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: C_mat, S_der, s_derivs
    TYPE(cp_dbcsr_type), POINTER             :: B_mat, tmp_mat, tmp_mat2
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(rt_prop_type), POINTER              :: rtp

    failure = .FALSE.
    CALL timeset(routineN,handle)

    return_s_derivatives=.TRUE.

    NULLIFY (atomic_kind_set)
    NULLIFY (force)
    NULLIFY (particle_set)
    NULLIFY (rtp)
    NULLIFY (s_derivs)

    CALL get_qs_env(qs_env=qs_env,&
                    rtp=rtp,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    neighbor_list_id=neighbor_list_id,&
                    sab_orb=sab_orb,&
                    force=force,&
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
       IF(qs_env%dft_control%qs_control%dftb)THEN
          CALL build_dftb_overlap(qs_env,nder,s_derivs,error)
       ELSE
          CALL build_overlap_matrix(qs_env,nderivative=nder,matrix_s=s_derivs,&
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
             IF(row_atom==col_atom) block_values=0
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
          IF(row_atom==col_atom) block_values=0
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
!> \author Florian Schiffmann (02.09)
! *****************************************************************************


  SUBROUTINE get_restart_wfn(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_restart_wfn', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, id_nr, im, ispin, ncol, &
                                                nspin, re
    REAL(KIND=dp)                            :: alpha
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_rmpv, rho_old
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new, mos_old
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(section_vals_type), POINTER         :: dft_section, input

    NULLIFY(atomic_kind_set,mo_array,particle_set)

    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,&
         particle_set=particle_set,mos=mo_array, &
         input=input,rtp=rtp,dft_control=dft_control,&
         error=error)

    id_nr=0
    nspin=SIZE(mo_array)
    p_rmpv => qs_env%rho%rho_ao
    dft_section =>  section_vals_get_subs_vals(input,"DFT",error=error)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)


    SELECT CASE(dft_control%rtp_control%initial_wfn)
    CASE(use_restart_wfn)
       CALL read_mo_set(mo_array,atomic_kind_set,particle_set,qs_env%para_env,&
            id_nr=id_nr,multiplicity=dft_control%multiplicity,dft_section=dft_section,&
            error=error)
       DO ispin=1,nspin
          CALL calculate_density_matrix(mo_array(ispin)%mo_set, p_rmpv(ispin)%matrix,error=error)
       ENDDO
       IF(rtp%linear_scaling) THEN
          CALL get_rtp(rtp=rtp,rho_old=rho_old,error=error)
          DO ispin=1,SIZE(mo_array)
             re=2*ispin-1
             im=2*ispin
             CALL cp_fm_get_info(mo_array(ispin)%mo_set%mo_coeff,ncol_global=ncol,error=error)
             alpha=1.0_dp
             IF(SIZE(mo_array)==1) alpha=2*alpha
             CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho_old(re)%matrix,&
                  matrix_v=mo_array(ispin)%mo_set%mo_coeff,matrix_g=mo_array(ispin)%mo_set%mo_coeff,ncol=ncol,&
                  keep_sparsity=.FALSE.,alpha=alpha,error=error)
          END DO
       ELSE
          CALL get_rtp(rtp=rtp,mos_old=mos_old,error=error)
          DO i=1,SIZE(qs_env%mos)
             CALL cp_fm_to_fm(mo_array(i)%mo_set%mo_coeff,mos_old(2*i-1)%matrix,error)
             CALL cp_fm_set_all(mos_old(2*i)%matrix,zero,zero,error)
          END DO
       ENDIF
    CASE(use_rt_restart)
       CALL get_rtp(rtp=rtp,mos_old=mos_old,mos_new=mos_new,error=error)
       CALL read_rt_mos_from_restart(mo_array,mos_old,atomic_kind_set,particle_set,qs_env%para_env,&
            id_nr,dft_control%multiplicity,dft_section, error)
       DO ispin=1,nspin
          CALL calculate_density_matrix(mo_array(ispin)%mo_set,&
               p_rmpv(ispin)%matrix,error=error)
       ENDDO
    END SELECT




  END SUBROUTINE get_restart_wfn

! *****************************************************************************
!> \brief calculates the imaginary part of the density matrix
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE calculate_dp_imag(dp_imag,mos_new,mo_set,error)
    TYPE(cp_fm_type), POINTER                :: dp_imag
    TYPE(cp_fm_p_type), DIMENSION(2)         :: mos_new
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dp_imag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao, nmo
    REAL(KIND=dp)                            :: alpha
    TYPE(cp_fm_type), POINTER                :: fm_tmp

    CALL timeset(routineN,handle)

    IF ( .NOT. mo_set%uniform_occupation ) THEN ! not all orbitals 1..homo are equally occupied
       NULLIFY(fm_tmp)
       CALL cp_fm_create(fm_tmp,mo_set%mo_coeff%matrix_struct,error=error)
       CALL cp_fm_to_fm(mos_new(2)%matrix,fm_tmp,error=error)
       CALL cp_fm_column_scale(fm_tmp,mo_set%occupation_numbers(1:mo_set%homo))
       alpha=1.0_dp
    ELSE
       alpha=mo_set%maxocc
       CALL cp_fm_create(fm_tmp,mo_set%mo_coeff%matrix_struct,error=error)
       CALL cp_fm_to_fm(mos_new(2)%matrix,fm_tmp,error=error)
    ENDIF
    CALL get_mo_set(mo_set, nao=nao,nmo=nmo)

    CALL cp_gemm("N","T",nao,nao,nmo,alpha,mos_new(1)%matrix,fm_tmp,0.0_dp,dp_imag,error)

    CALL cp_gemm("N","T",nao,nao,nmo,-alpha,fm_tmp,mos_new(1)%matrix,1.0_dp,dp_imag,error)

    CALL cp_fm_release(fm_tmp,error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_dp_imag

! *****************************************************************************
!> \brief calculates the three additional force contributions needed in EMD
!>        P_imag*C , P_imag*B*S^-1*S_der , P*S^-1*H*S_der
!>        driver routine
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
! *****************************************************************************
  SUBROUTINE calc_c_mat_force_fm(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_c_mat_force_fm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, im, ispin, nao, &
                                                natom, nmo, re, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha
    TYPE(array_i1d_obj)                      :: col_blk_size
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
       
       CALL dbcsr_create_dist_r_unrot (dist, SinvB(ispin)%matrix%matrix%m%dist, nmo, col_blk_size)
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
       
       CALL array_release(col_blk_size)
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
!> \par History
!>      02.2014 switched to dbcsr matrices [Samuel Andermatt]
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE calc_c_mat_force_ls(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_c_mat_force_ls', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, im, ispin, n, &
                                                natom, re, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: one = 1.0_dp, zero = 0.0_dp
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: C_mat, rho_new, S_der, SinvB, &
                                                SinvH
    TYPE(cp_dbcsr_type), POINTER             :: S_inv, tmp
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(rt_prop_type), POINTER              :: rtp

    failure=.FALSE.

    CALL timeset(routineN,handle)

    NULLIFY(rtp,particle_set,atomic_kind_set)
    CALL get_qs_env(qs_env=qs_env,rtp=rtp,particle_set=particle_set,&
                 atomic_kind_set=atomic_kind_set,force=force,error=error)

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

    CALL get_rtp(rtp=rtp,rho_new=rho_new,error=error)

    DO ispin=1,SIZE(SinvH)
       re=2*ispin-1
       im=2*ispin

       CALL cp_dbcsr_multiply("N","N",one,SinvH(ispin)%matrix,rho_new(re)%matrix,zero,tmp,&
            filter_eps=rtp%filter_eps,error=error)
       CALL cp_dbcsr_multiply("N","N",-one,SinvB(ispin)%matrix,rho_new(im)%matrix,one,tmp,filter_eps=rtp%filter_eps,error=error)

       CALL compute_forces(force,tmp,S_der,rho_new(im)%matrix,C_mat,kind_of,atom_of_kind,error)

    END DO

    ! recall QS forces, at this point have the other sign.
    DO i=1,SIZE(force)
       force(i)%ehrenfest(:,:)=- force(i)%ehrenfest(:,:)
    END DO

    CALL cp_dbcsr_deallocate_matrix(tmp,error=error)

    DEALLOCATE (atom_of_kind,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (kind_of,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE

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
          ALLOCATE(block_values2(SIZE(block_values)))
          CALL cp_dbcsr_get_block(S_der(i)%matrix, row_atom, col_atom, block_values2, found=found)
          IF(found) THEN
             ikind=kind_of(col_atom)
             kind_atom=atom_of_kind(col_atom)
             !The block_values are in a vector format, so the dot_product is the sum over all elements of the hamand product, that I need
             force(ikind)%ehrenfest(i,kind_atom)=force(ikind)%ehrenfest(i,kind_atom)+&
             2.0_dp*DOT_PRODUCT(block_values,block_values2)
          ENDIF
          DEALLOCATE(block_values2)
       END DO
       CALL cp_dbcsr_iterator_stop (iter)

       !C_mat part

       CALL cp_dbcsr_iterator_start(iter, rho_im)
       DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
          CALL cp_dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
          ALLOCATE(block_values2(SIZE(block_values)))
          CALL cp_dbcsr_get_block(C_mat(i)%matrix,row_atom,col_atom,block_values2,found=found)
          IF(found) THEN
             ikind=kind_of(col_atom)
             kind_atom=atom_of_kind(col_atom)
             !The block_values are in a vector format, so the dot_product is the sum over all elements of the hamand product, that I need
             force(ikind)%ehrenfest(i,kind_atom)=force(ikind)%ehrenfest(i,kind_atom)+&
             2.0_dp*DOT_PRODUCT(block_values,block_values2)
          ENDIF
          DEALLOCATE(block_values2)
       END DO
       CALL cp_dbcsr_iterator_stop (iter)
    END DO

 END SUBROUTINE compute_forces

END MODULE rt_propagation_utils
