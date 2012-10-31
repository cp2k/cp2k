!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines to calculate RI-RPA energy
!> \par History
!>      06.2012 created [Mauro Del Ben]
! *****************************************************************************
MODULE rpa_ri_gpw
! USE f77_blas
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale_and_add
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush,&
                                             m_memory
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_comm_split_direct,&
                                             mp_min,&
                                             mp_sendrecv,&
                                             mp_sum
  USE mp2_types,                       ONLY: mp2_type
  USE qs_environment_types,            ONLY: qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_ri_gpw'

  PUBLIC :: rpa_ri_compute_en


  CONTAINS

  SUBROUTINE rpa_ri_compute_en(Erpa,mp2_env,BIb_C,qs_env,para_env,para_env_sub,color_sub,&
                               ends_array,ends_B_virtual,sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&
                               Eigenval,nmo,homo,dimen_RI,unit_nr,error,error_sub,&
                               BIb_C_beta,homo_beta,Eigenval_beta,&
                               ends_B_virtual_beta,sizes_B_virtual_beta, starts_B_virtual_beta)
    REAL(KIND=dp)                            :: Erpa
    TYPE(mp2_type), POINTER                  :: mp2_env
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb_C
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub
    INTEGER                                  :: color_sub
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ends_array, ends_B_virtual, &
      sizes_array, sizes_B_virtual, starts_array, starts_B_virtual
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    INTEGER                                  :: nmo, homo, dimen_RI, unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error, error_sub
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :), OPTIONAL           :: BIb_C_beta
    INTEGER, OPTIONAL                        :: homo_beta
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: Eigenval_beta
    INTEGER, ALLOCATABLE, DIMENSION(:), &
      OPTIONAL                               :: ends_B_virtual_beta, &
                                                sizes_B_virtual_beta, &
                                                starts_B_virtual_beta

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rpa_ri_compute_en', &
      routineP = moduleN//':'//routineN

    INTEGER :: a, best_integ_group_size, best_num_integ_point, &
      color_rpa_group, dimen_ia, dimen_ia_beta, handle, i, ia_global, iiB, &
      integ_group_size, iproc, itmp(2), jjB, min_integ_group_size, my_B_size, &
      my_B_size_beta, my_B_virtual_end, my_B_virtual_end_beta, &
      my_B_virtual_start, my_B_virtual_start_beta, my_group_L_end, &
      my_group_L_size, my_group_L_start, my_homo_beta, my_ia_end, &
      my_ia_end_beta, my_ia_size, my_ia_size_beta, my_ia_start, &
      my_ia_start_beta, ngroup, num_integ_group, num_integ_points, &
      pos_integ_group, proc_receive, proc_send, proc_shift, rec_B_size, &
      rec_B_virtual_end, rec_B_virtual_start, stat, virtual, virtual_beta
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ends_ia, ends_ia_beta, &
                                                sizes_ia, sizes_ia_beta, &
                                                starts_ia, starts_ia_beta, &
                                                sub_proc_map
    LOGICAL                                  :: failure, my_open_shell
    REAL(KIND=dp) :: avail_mem, mem_for_iaK, mem_for_QK, mem_min, &
      mem_per_group, mem_real, needed_mem, rpa_memory
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: BIb_C_2D, BIb_C_2D_beta
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb_C_rec
    TYPE(cp_fm_type), POINTER                :: fm_mat_G, fm_mat_G_beta, &
                                                fm_mat_Q, fm_mat_Q_beta, &
                                                fm_mat_S, fm_mat_S_beta
    TYPE(cp_para_env_type), POINTER          :: para_env_RPA

    CALL timeset(routineN,handle)
    failure=.FALSE.

    num_integ_points=mp2_env%ri_rpa%rpa_num_quad_points
    rpa_memory=mp2_env%mp2_memory

    my_open_shell=.FALSE.
    IF(PRESENT(BIb_C_beta).AND.&
       PRESENT(ends_B_virtual_beta).AND.&
       PRESENT(sizes_B_virtual_beta).AND.&
       PRESENT(starts_B_virtual_beta).AND.&
       PRESENT(homo_beta).AND.&
       PRESENT(Eigenval_beta)) my_open_shell=.TRUE.

    virtual=nmo-homo

    ngroup=para_env%num_pe/para_env_sub%num_pe

    ! Calculate available memory and create integral group according to that
    ! mem_for_iaK is the memory needed for storing the 3 centre integrals
    mem_for_iaK=REAL(homo,KIND=dp)*virtual*dimen_RI*8.0_dp/(1024_dp**2)
    mem_for_QK=REAL(dimen_RI,KIND=dp)*dimen_RI*8.0_dp/(1024_dp**2)

    IF(my_open_shell) THEN
      virtual_beta=nmo-homo_beta
      mem_for_iaK=mem_for_iaK+REAL(homo_beta,KIND=dp)*virtual_beta*dimen_RI*8.0_dp/(1024_dp**2)
      mem_for_QK=mem_for_QK*2.0_dp
    END IF

    mem_real=m_memory()
    mem_real=(mem_real+1024*1024-1)/(1024*1024)
    ! mp_min .... a hack.. it should be mp_max, but as it turns out, on some processes the previously freed memory (hfx)
    ! has not been given back to the OS yet. 
    CALL mp_min(mem_real,para_env%group)

    mem_min=REAL(homo,KIND=dp)*MAXVAL(sizes_B_virtual)*MAXVAL(sizes_array)*8.0_dp/(1024**2)
 
    IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T69,F9.2,A3)') 'RI_INFO| Minimum required memory per MPI process:',&
                                                        mem_min, ' MB'
   
    mem_real=rpa_memory-mem_real
    mem_real=MAX(mem_real,mem_min)

    IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T69,F9.2,A3)') 'RI_INFO| Available memory per MPI process:',&
                                                        mem_real, ' MB'
    
    mem_per_group=mem_real*para_env_sub%num_pe

    needed_mem=mem_for_iaK*2.0_dp+mem_for_QK

    ! here we try to find the best rpa group size
    best_integ_group_size=ngroup
    best_num_integ_point=num_integ_points

    min_integ_group_size=MAX(1,ngroup/num_integ_points)

    integ_group_size=min_integ_group_size-1
    DO iiB=min_integ_group_size+1, ngroup
      integ_group_size=integ_group_size+1

      ! check that the ngroup is a multiple of  integ_group_size
      IF(MOD(ngroup,integ_group_size)/=0) CYCLE

      ! check for memory 
      avail_mem=integ_group_size*mem_per_group
      IF(avail_mem<needed_mem) CYCLE
      
      ! check the number of integration points is a multiple of the  number of integ_group
      num_integ_group=ngroup/integ_group_size
      IF(num_integ_points<num_integ_group) CYCLE
      IF(MOD(num_integ_points,num_integ_group)/=0) CYCLE

      ! if all the test passed then decide
      IF((num_integ_points/num_integ_group)<best_num_integ_point) THEN 
        best_num_integ_point=num_integ_points/num_integ_group
        best_integ_group_size=integ_group_size
      END IF

    END DO

    integ_group_size=best_integ_group_size

    IF (unit_nr>0) THEN
       WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                               "RI_INFO| Number of group for (ia|K) integral replication:", integ_group_size
       WRITE (UNIT=unit_nr,FMT="(T3,A)")&
                               "INTEG_INFO| Clenshaw-Curtius quadrature"
       WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                               "INTEG_INFO| Number of integration points:", num_integ_points
       WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                               "INTEG_INFO| Number of integration points per RPA group:", best_num_integ_point
       CALL m_flush(unit_nr)
    END IF
    
    num_integ_group=ngroup/integ_group_size

    pos_integ_group=MOD(color_sub,integ_group_size)
    color_rpa_group=color_sub/integ_group_size

    my_group_L_size=sizes_array(color_sub)
    my_group_L_start=starts_array(color_sub)
    my_group_L_end=ends_array(color_sub)

    my_B_virtual_start=starts_B_virtual(para_env_sub%mepos)
    my_B_virtual_end=ends_B_virtual(para_env_sub%mepos)
    my_B_size=sizes_B_virtual(para_env_sub%mepos)

    IF(my_open_shell) THEN
      my_B_virtual_start_beta=starts_B_virtual_beta(para_env_sub%mepos)
      my_B_virtual_end_beta=ends_B_virtual_beta(para_env_sub%mepos)
      my_B_size_beta=sizes_B_virtual_beta(para_env_sub%mepos)
      my_homo_beta=homo_beta
    ELSE
      my_B_virtual_start_beta=my_B_virtual_start
      my_B_virtual_end_beta=my_B_virtual_end
      my_B_size_beta=my_B_size
      my_homo_beta=homo
    END IF

    ! reorder the local data in such a way to help the next stage of matrix creation
    ! now the data inside the group are divided into a ia x K matrix
    dimen_ia=homo*virtual
    ALLOCATE(sizes_ia(0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    sizes_ia=0
    ALLOCATE(starts_ia(0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    starts_ia=0
    ALLOCATE(ends_ia(0:para_env_sub%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ends_ia=0

    DO iproc=0, para_env_sub%num_pe-1
      itmp=get_limit(dimen_ia,para_env_sub%num_pe,iproc)
      starts_ia(iproc)=itmp(1)
      ends_ia(iproc)=itmp(2)
      sizes_ia(iproc)=itmp(2)-itmp(1)+1
    END DO

    my_ia_start=starts_ia(para_env_sub%mepos)
    my_ia_end=ends_ia(para_env_sub%mepos)
    my_ia_size=sizes_ia(para_env_sub%mepos)

    ! create the sub_proc_map
    ALLOCATE(sub_proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe-1))
    DO i=0,para_env_sub%num_pe-1
      sub_proc_map(i)=i
      sub_proc_map(-i-1)=para_env_sub%num_pe-i-1
      sub_proc_map(para_env_sub%num_pe+i)=i
    END DO

    ! reorder data
    ALLOCATE(BIb_C_2D(my_ia_size,my_group_L_size),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    BIb_C_2D=0.0_dp

    DO iiB=1, homo
      DO jjB=1, my_B_size
        ia_global=(iiB-1)*virtual+my_B_virtual_start+jjB-1
        IF(ia_global>=my_ia_start.AND.ia_global<=my_ia_end) THEN
          BIb_C_2D(ia_global-my_ia_start+1,1:my_group_L_size)=BIb_C(1:my_group_L_size,jjB,iiB)
        END IF
      END DO
    END DO

    DO proc_shift=1, para_env_sub%num_pe-1
      proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
      proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)

      rec_B_size=sizes_B_virtual(proc_receive)
      rec_B_virtual_end=ends_B_virtual(proc_receive)
      rec_B_virtual_start=starts_B_virtual(proc_receive)

      ALLOCATE(BIb_C_rec(my_group_L_size,rec_B_size,homo),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      BIb_C_rec=0.0_dp
      
      CALL  mp_sendrecv(BIb_C,proc_send,&
                        BIb_C_rec,proc_receive,&
                        para_env_sub%group)

      DO iiB=1, homo
        DO jjB=1, rec_B_size
          ia_global=(iiB-1)*virtual+rec_B_virtual_start+jjB-1
          IF(ia_global>=my_ia_start.AND.ia_global<=my_ia_end) THEN
            BIb_C_2D(ia_global-my_ia_start+1,1:my_group_L_size)=BIb_C_rec(1:my_group_L_size,jjB,iiB)
          END IF
        END DO
      END DO

      DEALLOCATE(BIb_C_rec)
    END DO

    DEALLOCATE(BIb_C)
    DEALLOCATE(starts_B_virtual)
    DEALLOCATE(ends_B_virtual)
    DEALLOCATE(sizes_B_virtual)

    IF(my_open_shell) THEN
      ! do the same for the beta component
      dimen_ia_beta=homo_beta*virtual_beta
      ALLOCATE(sizes_ia_beta(0:para_env_sub%num_pe-1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      sizes_ia_beta=0
      ALLOCATE(starts_ia_beta(0:para_env_sub%num_pe-1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      starts_ia_beta=0
      ALLOCATE(ends_ia_beta(0:para_env_sub%num_pe-1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ends_ia_beta=0

      DO iproc=0, para_env_sub%num_pe-1
        itmp=get_limit(dimen_ia_beta,para_env_sub%num_pe,iproc)
        starts_ia_beta(iproc)=itmp(1)
        ends_ia_beta(iproc)=itmp(2)
        sizes_ia_beta(iproc)=itmp(2)-itmp(1)+1
      END DO

      my_ia_start_beta=starts_ia_beta(para_env_sub%mepos)
      my_ia_end_beta=ends_ia_beta(para_env_sub%mepos)
      my_ia_size_beta=sizes_ia_beta(para_env_sub%mepos)

      ! reorder data
      ALLOCATE(BIb_C_2D_beta(my_ia_size_beta,my_group_L_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      BIb_C_2D_beta=0.0_dp

      DO iiB=1, homo_beta
        DO jjB=1, my_B_size_beta
          ia_global=(iiB-1)*virtual_beta+my_B_virtual_start_beta+jjB-1
          IF(ia_global>=my_ia_start_beta.AND.ia_global<=my_ia_end_beta) THEN
            BIb_C_2D_beta(ia_global-my_ia_start_beta+1,1:my_group_L_size)=BIb_C_beta(1:my_group_L_size,jjB,iiB)
          END IF
        END DO
      END DO

      DO proc_shift=1, para_env_sub%num_pe-1
        proc_send=sub_proc_map(para_env_sub%mepos+proc_shift)
        proc_receive=sub_proc_map(para_env_sub%mepos-proc_shift)

        rec_B_size=sizes_B_virtual_beta(proc_receive)
        rec_B_virtual_end=ends_B_virtual_beta(proc_receive)
        rec_B_virtual_start=starts_B_virtual_beta(proc_receive)

        ALLOCATE(BIb_C_rec(my_group_L_size,rec_B_size,homo_beta),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
        BIb_C_rec=0.0_dp

        CALL  mp_sendrecv(BIb_C_beta,proc_send,&
                          BIb_C_rec,proc_receive,&
                          para_env_sub%group)

        DO iiB=1, homo_beta
          DO jjB=1, rec_B_size
            ia_global=(iiB-1)*virtual_beta+rec_B_virtual_start+jjB-1
            IF(ia_global>=my_ia_start_beta.AND.ia_global<=my_ia_end_beta) THEN
              BIb_C_2D_beta(ia_global-my_ia_start_beta+1,1:my_group_L_size)=BIb_C_rec(1:my_group_L_size,jjB,iiB)
            END IF
          END DO
        END DO

        DEALLOCATE(BIb_C_rec)
      END DO

      DEALLOCATE(BIb_C_beta)
      DEALLOCATE(starts_B_virtual_beta)
      DEALLOCATE(ends_B_virtual_beta)
      DEALLOCATE(sizes_B_virtual_beta)
    END IF

    ! now create the matrices needed for the calculation, Q, S and G
    ! Q and G will have omega dependence
    CALL create_integ_mat(BIb_C_2D,para_env,ngroup,integ_group_size,num_integ_group,&
                          dimen_RI,homo,virtual,dimen_ia,color_rpa_group,&
                          sizes_ia,starts_ia,ends_ia,my_ia_size,my_ia_start,my_ia_end,&
                          sizes_array,starts_array,ends_array,my_group_L_size,my_group_L_start,my_group_L_end,&
                          para_env_RPA,fm_mat_S,fm_mat_G,fm_mat_Q,error,error_sub)

    IF(my_open_shell) THEN
      ! the same for beta
      CALL create_integ_mat(BIb_C_2D_beta,para_env,ngroup,integ_group_size,num_integ_group,&
                            dimen_RI,homo_beta,virtual_beta,dimen_ia_beta,color_rpa_group,&
                            sizes_ia_beta,starts_ia_beta,ends_ia_beta,my_ia_size_beta,my_ia_start_beta,my_ia_end_beta,&
                            sizes_array,starts_array,ends_array,my_group_L_size,my_group_L_start,my_group_L_end,&
                            para_env_RPA,fm_mat_S_beta,fm_mat_G_beta,fm_mat_Q_beta,error,error_sub,&
                            .TRUE.,fm_mat_Q%matrix_struct%context)

      ! go with clenshaw curtius quadrature
      CALL clenshaw_curtius_quad(Erpa,para_env,para_env_RPA,unit_nr,homo,virtual,dimen_RI,dimen_ia,Eigenval,&
                                 num_integ_points,num_integ_group,color_rpa_group,&
                                 fm_mat_S,fm_mat_G,fm_mat_Q,error,&
                                 homo_beta,virtual_beta,dimen_ia_beta,Eigenval_beta,fm_mat_S_beta,fm_mat_G_beta,fm_mat_Q_beta)

    ELSE
      CALL clenshaw_curtius_quad(Erpa,para_env,para_env_RPA,unit_nr,homo,virtual,dimen_RI,dimen_ia,Eigenval,&
                                 num_integ_points,num_integ_group,color_rpa_group,&
                                 fm_mat_S,fm_mat_G,fm_mat_Q,error)
    END IF

    DEALLOCATE(sizes_ia)
    DEALLOCATE(starts_ia)
    DEALLOCATE(ends_ia)

    DEALLOCATE(sub_proc_map)

    DEALLOCATE(sizes_array)
    DEALLOCATE(starts_array)
    DEALLOCATE(ends_array)

    CALL cp_para_env_release(para_env_RPA,error=error_sub)

    CALL cp_fm_release(fm_mat_S,error=error)
    CALL cp_fm_release(fm_mat_G,error=error)
    CALL cp_fm_release(fm_mat_Q,error=error)
    IF(my_open_shell) THEN
      DEALLOCATE(sizes_ia_beta)
      DEALLOCATE(starts_ia_beta)
      DEALLOCATE(ends_ia_beta)
      CALL cp_fm_release(fm_mat_S_beta,error=error)
      CALL cp_fm_release(fm_mat_G_beta,error=error)
      CALL cp_fm_release(fm_mat_Q_beta,error=error)
    END IF

    CALL timestop(handle)

    END SUBROUTINE rpa_ri_compute_en

    SUBROUTINE create_integ_mat(BIb_C_2D,para_env,ngroup,integ_group_size,num_integ_group,&
                                dimen_RI,homo,virtual,dimen_ia,color_rpa_group,&
                                sizes_ia,starts_ia,ends_ia,my_ia_size,my_ia_start,my_ia_end,&
                                sizes_array,starts_array,ends_array,my_group_L_size,my_group_L_start,my_group_L_end,&
                                para_env_RPA,fm_mat_S,fm_mat_G,fm_mat_Q,error,error_sub,beta_case,blacs_env_ext)

    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: BIb_C_2D
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER :: ngroup, integ_group_size, num_integ_group, dimen_RI, homo, &
      virtual, dimen_ia, color_rpa_group
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sizes_ia, starts_ia, ends_ia
    INTEGER                                  :: my_ia_size, my_ia_start, &
                                                my_ia_end
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sizes_array, starts_array, &
                                                ends_array
    INTEGER                                  :: my_group_L_size, &
                                                my_group_L_start, &
                                                my_group_L_end
    TYPE(cp_para_env_type), POINTER          :: para_env_RPA
    TYPE(cp_fm_type), POINTER                :: fm_mat_S, fm_mat_G, fm_mat_Q
    TYPE(cp_error_type), INTENT(inout)       :: error, error_sub
    LOGICAL, OPTIONAL                        :: beta_case
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env_ext

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_integ_mat', &
      routineP = moduleN//':'//routineN

    INTEGER :: comm_exchange, comm_rpa, handle, i, i_global, iiB, j_global, &
      jjB, ncol_local, nrow_local, proc_receive, proc_receive_static, &
      proc_send, proc_send_static, proc_shift, rec_ia_end, rec_ia_size, &
      rec_ia_start, rec_L_end, rec_L_size, rec_L_start, stat, sub_sub_color
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: RPA_proc_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: RPA_info_end, RPA_info_size, &
                                                RPA_info_start
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure, my_beta_case, &
                                                my_blacs_ext
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: BIb_C_rec
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_para_env_type), POINTER          :: para_env_exchange

    CALL timeset(routineN,handle)
    failure=.FALSE.

    my_beta_case=.FALSE.
    IF(PRESENT(beta_case)) my_beta_case=beta_case

    my_blacs_ext=.FALSE.
    IF(PRESENT(blacs_env_ext)) my_blacs_ext=.TRUE.

    ! create the RPA para_env
    IF(.NOT.my_beta_case) THEN
      CALL mp_comm_split_direct(para_env%group,comm_rpa,color_rpa_group)
      NULLIFY(para_env_RPA)
      CALL cp_para_env_create(para_env_RPA,comm_rpa,error=error)
    END IF

    ! create the RPA blacs env 
    NULLIFY(blacs_env)
    IF(my_blacs_ext) THEN 
      blacs_env=>blacs_env_ext
    ELSE
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env_RPA, error=error)
    END IF

    ! create the S full matrix, that is the (ia|K) matrix with K colomns
    ! and homo*virtual rows
    NULLIFY(fm_mat_S)
    NULLIFY(fm_struct)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=dimen_ia,&
                             ncol_global=dimen_RI,para_env=para_env_RPA,error=error)
    CALL cp_fm_create(fm_mat_S,fm_struct,name="fm_mat_S",error=error)

    CALL cp_fm_set_all(matrix=fm_mat_S,alpha=0.0_dp,error=error)

    ! fill the matrix
    CALL cp_fm_get_info(matrix=fm_mat_S,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        error=error)

    ! create the RPA proc_map
    ALLOCATE(RPA_proc_map(-para_env_RPA%num_pe:2*para_env_RPA%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    RPA_proc_map=0
    DO i=0,para_env_RPA%num_pe-1
      RPA_proc_map(i)=i
      RPA_proc_map(-i-1)=para_env_RPA%num_pe-i-1
      RPA_proc_map(para_env_RPA%num_pe+i)=i
    END DO

    ! create the info array, first index: 1-> L info, 2-> ia info
    ALLOCATE(RPA_info_start(2,0:para_env_RPA%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    RPA_info_start=0
    ALLOCATE(RPA_info_end(2,0:para_env_RPA%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    RPA_info_end=0
    ALLOCATE(RPA_info_size(2,0:para_env_RPA%num_pe-1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    RPA_info_size=0

    RPA_info_start(1,para_env_RPA%mepos)=my_group_L_start
    RPA_info_start(2,para_env_RPA%mepos)=my_ia_start

    RPA_info_end(1,para_env_RPA%mepos)=my_group_L_end
    RPA_info_end(2,para_env_RPA%mepos)=my_ia_end

    RPA_info_size(1,para_env_RPA%mepos)=my_group_L_size
    RPA_info_size(2,para_env_RPA%mepos)=my_ia_size

    CALL mp_sum(RPA_info_start,para_env_RPA%group)
    CALL mp_sum(RPA_info_end,para_env_RPA%group)
    CALL mp_sum(RPA_info_size,para_env_RPA%group)

    ! local data
    DO jjB=1, ncol_local
      j_global=col_indices(jjB)
      IF(j_global>=my_group_L_start.AND.j_global<=my_group_L_end) THEN
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          IF(i_global>=my_ia_start.AND.i_global<=my_ia_end) THEN
            fm_mat_S%local_data(iiB,jjB)=BIb_C_2D(i_global-my_ia_start+1,j_global-my_group_L_start+1)
          END IF
        END DO
      END IF
    END DO

    proc_send_static=RPA_proc_map(para_env_RPA%mepos+1)
    proc_receive_static=RPA_proc_map(para_env_RPA%mepos-1)

    ! start collect data from other proc in the RPA group
    DO proc_shift=1, para_env_RPA%num_pe-1
      proc_send=RPA_proc_map(para_env_RPA%mepos+proc_shift)
      proc_receive=RPA_proc_map(para_env_RPA%mepos-proc_shift)

      rec_L_size=RPA_info_size(1,proc_receive)
      rec_L_start=RPA_info_start(1,proc_receive)
      rec_L_end=RPA_info_end(1,proc_receive)

      rec_ia_size=RPA_info_size(2,proc_receive)
      rec_ia_start=RPA_info_start(2,proc_receive)
      rec_ia_end=RPA_info_end(2,proc_receive)

      ALLOCATE(BIb_C_rec(rec_ia_size,rec_L_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      BIb_C_rec=0.0_dp

      CALL  mp_sendrecv(BIb_C_2D,proc_send_static,&
                        BIb_C_rec,proc_receive_static,&
                        para_env_RPA%group)

      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        IF(j_global>=rec_L_start.AND.j_global<=rec_L_end) THEN
          DO iiB=1, nrow_local
            i_global=row_indices(iiB)
            IF(i_global>=rec_ia_start.AND.i_global<=rec_ia_end) THEN
              fm_mat_S%local_data(iiB,jjB)=BIb_C_rec(i_global-rec_ia_start+1,j_global-rec_L_start+1)
            END IF
          END DO
        END IF
      END DO

      DEALLOCATE(BIb_C_2D)
      ALLOCATE(BIb_C_2D(rec_ia_size,rec_L_size),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      BIb_C_2D=BIb_C_rec

      DEALLOCATE(BIb_C_rec)
    END DO

    DEALLOCATE(BIb_C_2D)

    ! deallocaete the info array
    DEALLOCATE(RPA_info_start)
    DEALLOCATE(RPA_info_end)
    DEALLOCATE(RPA_info_size)

    ! mp_sum the local data across processes belonging to different RPA group.
    ! first create the para_env then mp_sum
    sub_sub_color=para_env_RPA%mepos
    CALL mp_comm_split_direct(para_env%group,comm_exchange,sub_sub_color)
    NULLIFY(para_env_exchange)
    CALL cp_para_env_create(para_env_exchange,comm_exchange,error=error_sub)

    CALL mp_sum(fm_mat_S%local_data,para_env_exchange%group)

    CALL cp_para_env_release(para_env_exchange,error=error_sub)

    ! create the twin matrix for the mat-mat-mul (mat_Q)
    CALL cp_fm_create(fm_mat_G,fm_struct,name="fm_mat_G",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)
    
    CALL cp_fm_set_all(matrix=fm_mat_G,alpha=0.0_dp,error=error)

    ! create the Q matrix dimen_RIxdimen_RI where the result of the mat-mat-mult will be stored
    NULLIFY(fm_mat_Q)
    NULLIFY(fm_struct)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=dimen_RI,&
                             ncol_global=dimen_RI,para_env=para_env_RPA,error=error)
    CALL cp_fm_create(fm_mat_Q,fm_struct,name="fm_mat_Q",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL cp_fm_set_all(matrix=fm_mat_Q,alpha=0.0_dp,error=error)


    ! release blacs_env
    IF(.NOT.my_blacs_ext) CALL cp_blacs_env_release(blacs_env,error=error)

    CALL timestop(handle)

    END SUBROUTINE create_integ_mat

    SUBROUTINE clenshaw_curtius_quad(Erpa,para_env,para_env_RPA,unit_nr,homo,virtual,dimen_RI,dimen_ia,Eigenval,&
                                     num_integ_points,num_integ_group,color_rpa_group,&
                                     fm_mat_S,fm_mat_G,fm_mat_Q,error,&
                                     homo_beta,virtual_beta,dimen_ia_beta,Eigenval_beta,fm_mat_S_beta,fm_mat_G_beta,fm_mat_Q_beta)
    REAL(KIND=dp)                            :: Erpa
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_RPA
    INTEGER                                  :: unit_nr, homo, virtual, &
                                                dimen_RI, dimen_ia
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    INTEGER                                  :: num_integ_points, &
                                                num_integ_group, &
                                                color_rpa_group
    TYPE(cp_fm_type), POINTER                :: fm_mat_S, fm_mat_G, fm_mat_Q
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, OPTIONAL                        :: homo_beta, virtual_beta, &
                                                dimen_ia_beta
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: Eigenval_beta
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: fm_mat_S_beta, fm_mat_G_beta, &
                                                fm_mat_Q_beta

    CHARACTER(LEN=*), PARAMETER :: routineN = 'clenshaw_curtius_quad', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: avirt, handle, i_global, iiB, &
                                                info_chol, iocc, j_global, &
                                                jjB, jquad, ncol_local, &
                                                nrow_local, stat
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure, my_open_shell
    REAL(KIND=dp)                            :: a_scaling, alpha, eigen_diff, &
                                                FComega, omega
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Q_log, tj, trace_Qomega, wj

    CALL timeset(routineN,handle)
    failure=.FALSE.
    
    my_open_shell=.FALSE.
    IF(PRESENT(homo_beta).AND.&
       PRESENT(virtual_beta).AND.&
       PRESENT(dimen_ia_beta).AND.&
       PRESENT(Eigenval_beta).AND.&
       PRESENT(fm_mat_S_beta).AND.&
       PRESENT(fm_mat_G_beta).AND.&
       PRESENT(fm_mat_Q_beta)) my_open_shell=.TRUE.

    ALLOCATE(tj(num_integ_points),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    tj=0.0_dp
    ALLOCATE(wj(num_integ_points),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    wj=0.0_dp
    ALLOCATE(Q_log(dimen_RI),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(trace_Qomega(dimen_RI),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO jquad=1, num_integ_points-1
      tj(jquad)=jquad*pi/(2.0_dp*num_integ_points)
      wj(jquad)=pi/(num_integ_points*SIN(tj(jquad))**2)
    END DO
    tj(num_integ_points)=pi/2.0_dp
    wj(num_integ_points)=pi/(2.0_dp*num_integ_points*SIN(tj(num_integ_points))**2)

    a_scaling=1.0_dp
    IF(my_open_shell) THEN
      CALL calc_scaling_factor(a_scaling,para_env,para_env_RPA,homo,virtual,dimen_RI,dimen_ia,Eigenval,&
                               num_integ_points,num_integ_group,color_rpa_group,&
                               tj,wj,fm_mat_S,fm_mat_G,error,&
                               homo_beta,virtual_beta,dimen_ia_beta,Eigenval_beta,fm_mat_S_beta,fm_mat_G_beta)
      alpha=2.0_dp
    ELSE
      CALL calc_scaling_factor(a_scaling,para_env,para_env_RPA,homo,virtual,dimen_RI,dimen_ia,Eigenval,&
                               num_integ_points,num_integ_group,color_rpa_group,&
                               tj,wj,fm_mat_S,fm_mat_G,error)
      alpha=4.0_dp
    END IF

    IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.5)') 'INTEG_INFO| Scaling parameter:', a_scaling

    Erpa=0.0_dp
    wj=wj*a_scaling
    DO jquad=1, num_integ_points
      IF(MODULO(jquad,num_integ_group)/=color_rpa_group) CYCLE
      
      omega=a_scaling/TAN(tj(jquad))

      ! copy fm_mat_S into fm_mat_G
      CALL cp_fm_to_fm(source=fm_mat_S,destination=fm_mat_G,error=error)

      ! get info of fm_mat_G
      CALL cp_fm_get_info(matrix=fm_mat_G,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices,&
                          error=error)

      ! update G matrix with the new value of omega
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)

          iocc=MAX(1,i_global-1)/virtual+1
          avirt=i_global-(iocc-1)*virtual
          eigen_diff=Eigenval(avirt+homo)-Eigenval(iocc)

          fm_mat_G%local_data(iiB,jjB)=fm_mat_G%local_data(iiB,jjB)*eigen_diff/(eigen_diff**2+omega**2)

        END DO
      END DO

      ! alpha=4.0 is valid only for closed shell systems (alpha is 2.0 in the spin-orbital basis)
      CALL cp_fm_gemm(transa="T",transb="N",m=dimen_RI,n=dimen_RI,k=dimen_ia,alpha=alpha,&
                      matrix_a=fm_mat_S,matrix_b=fm_mat_G,beta=0.0_dp,&
                      matrix_c=fm_mat_Q,error=error)

      IF(my_open_shell) THEN
        ! the same for the beta spin
        ! copy fm_mat_S into fm_mat_G
        CALL cp_fm_to_fm(source=fm_mat_S_beta,destination=fm_mat_G_beta,error=error)
        ! get info of fm_mat_G_beta
        CALL cp_fm_get_info(matrix=fm_mat_G_beta,&
                            nrow_local=nrow_local,&
                            ncol_local=ncol_local,&
                            row_indices=row_indices,&
                            col_indices=col_indices,&
                            error=error)
        ! update G matrix with the new value of omega
        DO jjB=1, ncol_local
          j_global=col_indices(jjB)
          DO iiB=1, nrow_local
            i_global=row_indices(iiB)

            iocc=MAX(1,i_global-1)/virtual_beta+1
            avirt=i_global-(iocc-1)*virtual_beta
            eigen_diff=Eigenval_beta(avirt+homo_beta)-Eigenval_beta(iocc)

            fm_mat_G_beta%local_data(iiB,jjB)=fm_mat_G_beta%local_data(iiB,jjB)*eigen_diff/(eigen_diff**2+omega**2)

          END DO
        END DO

        CALL cp_fm_gemm(transa="T",transb="N",m=dimen_RI,n=dimen_RI,k=dimen_ia_beta,alpha=alpha,&
                        matrix_a=fm_mat_S_beta,matrix_b=fm_mat_G_beta,beta=0.0_dp,&
                        matrix_c=fm_mat_Q_beta,error=error)

        CALL cp_fm_scale_and_add(alpha=1.0_dp,matrix_a=fm_mat_Q,beta=1.0_dp,matrix_b=fm_mat_Q_beta,error=error)
        ! fm_mat_Q%local_data=fm_mat_Q%local_data+fm_mat_Q_beta%local_data        

      END IF

      ! get info of fm_mat_Q
      CALL cp_fm_get_info(matrix=fm_mat_Q,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices,&
                          error=error)

      ! calcualte the trace of Q and add 1 on the diagonal
      trace_Qomega=0.0_dp
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          IF(j_global==i_global.AND.i_global<=dimen_RI) THEN
            trace_Qomega(i_global)=fm_mat_Q%local_data(iiB,jjB)
            fm_mat_Q%local_data(iiB,jjB)=fm_mat_Q%local_data(iiB,jjB)+1.0_dp
          END IF
        END DO
      END DO
      CALL mp_sum(trace_Qomega,para_env_RPA%group)

      ! calculate Trace(Log(Matrix)) as Log(DET(Matrix)) via cholesky decomposition
      CALL cp_fm_cholesky_decompose(matrix=fm_mat_Q, n=dimen_RI, info_out=info_chol, error=error)
      CPPostcondition(info_chol==0,cp_failure_level,routineP,error,failure)

      ! get info of cholesky_decomposed(fm_mat_Q)
      CALL cp_fm_get_info(matrix=fm_mat_Q,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices,&
                          error=error)

      Q_log=0.0_dp
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          IF(j_global==i_global.AND.i_global<=dimen_RI) THEN
            Q_log(i_global)=2.0_dp*LOG(fm_mat_Q%local_data(iiB,jjB))
          END IF
        END DO
      END DO
      CALL mp_sum(Q_log,para_env_RPA%group)

      FComega=0.0_dp
      DO iiB=1, dimen_RI
        IF(MODULO(iiB,para_env_RPA%num_pe)/=para_env_RPA%mepos) CYCLE 
        ! FComega=FComega+(LOG(Q_log(iiB))-trace_Qomega(iiB))/2.0_dp
        FComega=FComega+(Q_log(iiB)-trace_Qomega(iiB))/2.0_dp
      END DO

      Erpa=Erpa+FComega*wj(jquad)
    END DO

    CALL mp_sum(Erpa,para_env%group)
    Erpa=Erpa/(pi*2.0_dp)

    DEALLOCATE(tj)
    DEALLOCATE(wj)
    DEALLOCATE(Q_log)
    DEALLOCATE(trace_Qomega)

    CALL timestop(handle)

    END SUBROUTINE clenshaw_curtius_quad

    SUBROUTINE  calc_scaling_factor(a_scaling_ext,para_env,para_env_RPA,homo,virtual,dimen_RI,dimen_ia,Eigenval,&
                                    num_integ_points,num_integ_group,color_rpa_group,&
                                    tj_ext,wj_ext,fm_mat_S,fm_mat_G,error,&
                                    homo_beta,virtual_beta,dimen_ia_beta,Eigenval_beta,fm_mat_S_beta,fm_mat_G_beta)
    REAL(KIND=dp)                            :: a_scaling_ext
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_RPA
    INTEGER                                  :: homo, virtual, dimen_RI, &
                                                dimen_ia
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    INTEGER                                  :: num_integ_points, &
                                                num_integ_group, &
                                                color_rpa_group
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: tj_ext, wj_ext
    TYPE(cp_fm_type), POINTER                :: fm_mat_S, fm_mat_G
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, OPTIONAL                        :: homo_beta, virtual_beta, &
                                                dimen_ia_beta
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: Eigenval_beta
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: fm_mat_S_beta, fm_mat_G_beta

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_scaling_factor', &
      routineP = moduleN//':'//routineN

    INTEGER :: avirt, color_col, color_col_beta, color_row, color_row_beta, &
      comm_col, comm_col_beta, comm_row, comm_row_beta, handle, i_global, &
      icycle, iiB, iocc, jjB, jquad, ncol_local, ncol_local_beta, nrow_local, &
      nrow_local_beta, stat
    INTEGER, DIMENSION(:), POINTER           :: col_indices, &
                                                col_indices_beta, &
                                                row_indices, row_indices_beta
    LOGICAL                                  :: failure, my_open_shell
    REAL(KIND=dp) :: a_high, a_low, a_scaling, conv_param, eigen_diff, eps, &
      first_deriv, four, left_term, one, pig, right_term, right_term_ref, &
      right_term_ref_beta, step, two, zero
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cottj, D_ia, D_ia_beta, &
      iaia_RI, iaia_RI_beta, iaia_RI_dp, iaia_RI_dp_beta, M_ia, M_ia_beta, &
      tj, wj
    TYPE(cp_para_env_type), POINTER          :: para_env_col, &
                                                para_env_col_beta, &
                                                para_env_row, &
                                                para_env_row_beta

    CALL timeset(routineN,handle)
    failure=.FALSE.

    my_open_shell=.FALSE.
    IF(PRESENT(homo_beta).AND.&
       PRESENT(virtual_beta).AND.&
       PRESENT(dimen_ia_beta).AND.&
       PRESENT(Eigenval_beta).AND.&
       PRESENT(fm_mat_S_beta).AND.&
       PRESENT(fm_mat_G_beta)) my_open_shell=.TRUE.

    ZERO=0.0_dp
    ONE=1.0_dp
    TWO=2.0_dp
    FOUR=4.0_dp
    PIG=pi
    eps=1.0E-10_dp

    ALLOCATE(cottj(num_integ_points),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(tj(num_integ_points),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(wj(num_integ_points),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! calculate the cotangent of the abscissa tj
    DO jquad=1, num_integ_points
      tj(jquad)=tj_ext(jquad)
      wj(jquad)=wj_ext(jquad)
      cottj(jquad)=ONE/TAN(tj(jquad))
    END DO

    ! calculate the (ia|ia) RI integrals 
    ! ----------------------------------
    ! 1) copy fm_mat_S into fm_mat_G (working array)
    CALL cp_fm_to_fm(source=fm_mat_S,destination=fm_mat_G,error=error)
    CALL cp_fm_get_info(matrix=fm_mat_G,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        error=error)

    ! allocate the local buffer of iaia_RI integrals (dp kind)
    ALLOCATE(iaia_RI_dp(nrow_local),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    iaia_RI_dp=0.0_dp

    ! 2) perform the local multiplication SUM_K (ia|K)*(ia|K)
    DO jjB=1, ncol_local
      DO iiB=1, nrow_local
       iaia_RI_dp(iiB)=iaia_RI_dp(iiB)+fm_mat_G%local_data(iiB,jjB)*fm_mat_G%local_data(iiB,jjB)
      END DO
    END DO

    ! 3) sum the result with the processes of the RPA_group having the same rows
    !          _______K_______               _ 
    !         |   |   |   |   |             | |
    !     --> | 1 | 5 | 9 | 13|   SUM -->   | |
    !         |___|__ |___|___|             |_|
    !         |   |   |   |   |             | |
    !     --> | 2 | 6 | 10| 14|   SUM -->   | |
    !      ia |___|___|___|___|             |_|   (ia|ia)_RI
    !         |   |   |   |   |             | |
    !     --> | 3 | 7 | 11| 15|   SUM -->   | |
    !         |___|___|___|___|             |_|
    !         |   |   |   |   |             | |
    !     --> | 4 | 8 | 12| 16|   SUM -->   | |
    !         |___|___|___|___|             |_|
    !
   
    color_row=fm_mat_G%matrix_struct%context%mepos(1)
    CALL mp_comm_split_direct(para_env_RPA%group,comm_row,color_row)
    NULLIFY(para_env_row)
    CALL cp_para_env_create(para_env_row,comm_row,error=error)
    
    CALL mp_sum(iaia_RI_dp,para_env_row%group)

    ! convert the iaia_RI_dp into double-double precision
    ALLOCATE(iaia_RI(nrow_local),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO iiB=1, nrow_local
      iaia_RI(iiB)=iaia_RI_dp(iiB)
    END DO
    DEALLOCATE(iaia_RI_dp)


    ! 4) calculate the right hand term, D_ia is the matrix containing the 
    ! orbital energy differences, M_ia is the diagonal of the full RPA 'excitation'
    ! matrix 
    ALLOCATE(D_ia(nrow_local),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(M_ia(nrow_local),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO iiB=1, nrow_local
      i_global=row_indices(iiB)

      iocc=MAX(1,i_global-1)/virtual+1
      avirt=i_global-(iocc-1)*virtual
      eigen_diff=Eigenval(avirt+homo)-Eigenval(iocc)

      D_ia(iiB)=eigen_diff
    END DO

    DO iiB=1, nrow_local
      M_ia(iiB)=D_ia(iiB)*D_ia(iiB)+TWO*D_ia(iiB)*iaia_RI(iiB)
    END DO
 
    right_term_ref=ZERO
    DO iiB=1, nrow_local
      right_term_ref=right_term_ref+(SQRT(M_ia(iiB))-D_ia(iiB)-iaia_RI(iiB))
    END DO 
    right_term_ref=right_term_ref/TWO
    ! right_term_ref=accurate_sum((SQRT(M_ia)-D_ia-iaia_RI))/2.0_dp

    ! sum the result with the processes of the RPA_group having the same col
    color_col=fm_mat_G%matrix_struct%context%mepos(2)
    CALL mp_comm_split_direct(para_env_RPA%group,comm_col,color_col)
    NULLIFY(para_env_col)
    CALL cp_para_env_create(para_env_col,comm_col,error=error)

    ! allocate communication array for columns
    CALL mp_sum(right_term_ref,para_env_col%group)    

    ! In the open shell case do point 1-2-3 for the beta spin
    IF(my_open_shell) THEN
      CALL cp_fm_to_fm(source=fm_mat_S_beta,destination=fm_mat_G_beta,error=error)
      CALL cp_fm_get_info(matrix=fm_mat_G_beta,&
                          nrow_local=nrow_local_beta,&
                          ncol_local=ncol_local_beta,&
                          row_indices=row_indices_beta,&
                          col_indices=col_indices_beta,&
                          error=error)

      ALLOCATE(iaia_RI_dp_beta(nrow_local_beta),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      iaia_RI_dp_beta=0.0_dp

      DO jjB=1, ncol_local_beta
        DO iiB=1, nrow_local_beta
         iaia_RI_dp_beta(iiB)=iaia_RI_dp_beta(iiB)+fm_mat_G_beta%local_data(iiB,jjB)*fm_mat_G_beta%local_data(iiB,jjB)
        END DO
      END DO

      color_row_beta=fm_mat_G_beta%matrix_struct%context%mepos(1)
      CALL mp_comm_split_direct(para_env_RPA%group,comm_row_beta,color_row_beta)
      NULLIFY(para_env_row_beta)
      CALL cp_para_env_create(para_env_row_beta,comm_row_beta,error=error)

      CALL mp_sum(iaia_RI_dp_beta,para_env_row_beta%group)

      ALLOCATE(iaia_RI_beta(nrow_local_beta),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DO iiB=1, nrow_local_beta
        iaia_RI_beta(iiB)=iaia_RI_dp_beta(iiB)
      END DO
      DEALLOCATE(iaia_RI_dp_beta)

      ALLOCATE(D_ia_beta(nrow_local_beta),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      ALLOCATE(M_ia_beta(nrow_local_beta),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO iiB=1, nrow_local_beta
        i_global=row_indices_beta(iiB)

        iocc=MAX(1,i_global-1)/virtual_beta+1
        avirt=i_global-(iocc-1)*virtual_beta
        eigen_diff=Eigenval_beta(avirt+homo_beta)-Eigenval_beta(iocc)

        D_ia_beta(iiB)=eigen_diff
      END DO

      DO iiB=1, nrow_local_beta
        M_ia_beta(iiB)=D_ia_beta(iiB)*D_ia_beta(iiB)+TWO*D_ia_beta(iiB)*iaia_RI_beta(iiB)
      END DO

      right_term_ref_beta=ZERO
      DO iiB=1, nrow_local_beta
        right_term_ref_beta=right_term_ref_beta+(SQRT(M_ia_beta(iiB))-D_ia_beta(iiB)-iaia_RI_beta(iiB))
      END DO
      right_term_ref_beta=right_term_ref_beta/TWO

      ! sum the result with the processes of the RPA_group having the same col
      color_col_beta=fm_mat_G_beta%matrix_struct%context%mepos(2)
      CALL mp_comm_split_direct(para_env_RPA%group,comm_col_beta,color_col_beta)
      NULLIFY(para_env_col_beta)
      CALL cp_para_env_create(para_env_col_beta,comm_col_beta,error=error)

      CALL mp_sum(right_term_ref_beta,para_env_col_beta%group)

      right_term_ref=right_term_ref+right_term_ref_beta
    END IF

    ! bcast the result
    IF(para_env%mepos==0) THEN
      CALL mp_bcast(right_term_ref,0,para_env%group)
    ELSE
      right_term_ref=0.0_dp
      CALL mp_bcast(right_term_ref,0,para_env%group)
    END IF

    ! 5) start iteration for solving the non-linear equation by bisection
    ! find limit, here step=0.5 seems a good compromise
    conv_param=100.0_dp*EPSILON(right_term_ref)
    step=0.5_dp
    a_low=0.0_dp
    a_high=step
    right_term=-right_term_ref
    DO icycle=1, num_integ_points*2
      a_scaling=a_high
 
      CALL calculate_objfunc(a_scaling,left_term,first_deriv,num_integ_points,my_open_shell,&
                             ZERO,ONE,FOUR,PIG,M_ia,cottj,wj,D_ia,D_ia_beta,M_ia_beta,&
                             nrow_local,nrow_local_beta,num_integ_group,color_rpa_group,&
                             para_env,para_env_row,para_env_row_beta)
      left_term=left_term/FOUR/PIG*a_scaling

      IF(ABS(left_term)>ABS(right_term).OR.ABS(left_term+right_term)<=conv_param) EXIT
      a_low=a_high
      a_high=a_high+step

    END DO

    IF(ABS(left_term+right_term)>=conv_param) THEN
      IF(a_scaling>=2*num_integ_points*step) THEN
        a_scaling=1.0_dp
      ELSE

        DO icycle=1, num_integ_points*2
          a_scaling=(a_low+a_high)/2.0_dp

          CALL calculate_objfunc(a_scaling,left_term,first_deriv,num_integ_points,my_open_shell,&
                                 ZERO,ONE,FOUR,PIG,M_ia,cottj,wj,D_ia,D_ia_beta,M_ia_beta,&
                                 nrow_local,nrow_local_beta,num_integ_group,color_rpa_group,&
                                 para_env,para_env_row,para_env_row_beta)
          left_term=left_term/FOUR/PIG*a_scaling

          IF(ABS(left_term)>ABS(right_term)) THEN
            a_high=a_scaling
          ELSE
            a_low=a_scaling
          END IF

          ! IF(para_env%mepos==0) THEN
          !   WRITE(*,*) a_scaling, a_high-a_low, left_term+right_term
          ! END IF

          IF(ABS(a_high-a_low)<1.0D-5) EXIT

        END DO

      END IF
    END IF

    a_scaling_ext=a_scaling   
    CALL mp_bcast(a_scaling_ext,0,para_env%group)

    DEALLOCATE(cottj)
    DEALLOCATE(tj)
    DEALLOCATE(wj)
    DEALLOCATE(iaia_RI)
    DEALLOCATE(D_ia)
    DEALLOCATE(M_ia)
    CALL cp_para_env_release(para_env_row,error=error)   
    CALL cp_para_env_release(para_env_col,error=error)

    IF(my_open_shell) THEN
      DEALLOCATE(iaia_RI_beta)
      DEALLOCATE(D_ia_beta)
      DEALLOCATE(M_ia_beta)
      CALL cp_para_env_release(para_env_row_beta,error=error)
      CALL cp_para_env_release(para_env_col_beta,error=error)
    END IF

    CALL timestop(handle)

    END SUBROUTINE

    SUBROUTINE calculate_objfunc(a_scaling,left_term,first_deriv,num_integ_points,my_open_shell,&
                                 ZERO,ONE,FOUR,PIG,M_ia,cottj,wj,D_ia,D_ia_beta,M_ia_beta,&
                                 nrow_local,nrow_local_beta,num_integ_group,color_rpa_group,&
                                 para_env,para_env_row,para_env_row_beta)
    REAL(KIND=dp)                            :: a_scaling, left_term, &
                                                first_deriv
    INTEGER                                  :: num_integ_points
    LOGICAL                                  :: my_open_shell
    REAL(KIND=dp)                            :: zero, one, four, pig
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: M_ia, cottj, wj, D_ia, &
                                                D_ia_beta, M_ia_beta
    INTEGER                                  :: nrow_local, nrow_local_beta, &
                                                num_integ_group, &
                                                color_rpa_group
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_row, &
                                                para_env_row_beta

    INTEGER                                  :: iiB, jquad
    REAL(KIND=dp)                            :: first_deriv_beta, &
                                                left_term_beta, omega

      left_term=ZERO
      first_deriv=ZERO
      left_term_beta=ZERO
      first_deriv_beta=ZERO
      DO jquad=1, num_integ_points
        ! parallelize over integration points
        IF(MODULO(jquad,num_integ_group)/=color_rpa_group) CYCLE
        omega=a_scaling*cottj(jquad)

        DO iiB=1, nrow_local
          ! parallelize over ia elements in the para_env_row group
          IF(MODULO(iiB,para_env_row%num_pe)/=para_env_row%mepos) CYCLE
          ! calculate left_term
          left_term=left_term+wj(jquad)*&
                    (LOG(ONE+(M_ia(iiB)-D_ia(iiB)**2)/(omega**2+D_ia(iiB)**2))-&
                    (M_ia(iiB)-D_ia(iiB)**2)/(omega**2+D_ia(iiB)**2))
          first_deriv=first_deriv+wj(jquad)*cottj(jquad)**2*&
                      ((-M_ia(iiB)+D_ia(iiB)**2)**2/((omega**2+D_ia(iiB)**2)**2*(omega**2+M_ia(iiB))))
        END DO

        IF(my_open_shell) THEN
          DO iiB=1, nrow_local_beta
            ! parallelize over ia elements in the para_env_row group
            IF(MODULO(iiB,para_env_row_beta%num_pe)/=para_env_row_beta%mepos) CYCLE
            ! calculate left_term
            left_term_beta=left_term_beta+wj(jquad)*&
                      (LOG(ONE+(M_ia_beta(iiB)-D_ia_beta(iiB)**2)/(omega**2+D_ia_beta(iiB)**2))-&
                      (M_ia_beta(iiB)-D_ia_beta(iiB)**2)/(omega**2+D_ia_beta(iiB)**2))
            first_deriv_beta=first_deriv_beta+wj(jquad)*cottj(jquad)**2*&
                             ((-M_ia_beta(iiB)+D_ia_beta(iiB)**2)**2/((omega**2+D_ia_beta(iiB)**2)**2*(omega**2+M_ia_beta(iiB))))
          END DO
        END IF

      END DO
    
      ! sum the contribution from all proc, starting form the row group
      CALL mp_sum(left_term,para_env%group)
      CALL mp_sum(first_deriv,para_env%group)

      IF(my_open_shell) THEN
        CALL mp_sum(left_term_beta,para_env%group)
        CALL mp_sum(first_deriv_beta,para_env%group)

        left_term=left_term+left_term_beta
        first_deriv=first_deriv+first_deriv_beta
      END IF

    END SUBROUTINE calculate_objfunc

END MODULE rpa_ri_gpw
