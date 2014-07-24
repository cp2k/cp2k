!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines to calculate RI-RPA energy
!> \par History
!>      06.2012 created [Mauro Del Ben]
! *****************************************************************************
MODULE rpa_ri_gpw
! USE f77_blas
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release,&
                                             cp_blacs_env_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add,&
                                             cp_fm_syrk
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_indxg2l,&
                                             cp_fm_indxg2p,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: wfc_mm_style_gemm,&
                                             wfc_mm_style_syrk
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush,&
                                             m_memory,&
                                             m_walltime
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: &
       mp_bcast, mp_comm_split_direct, mp_irecv, mp_isend, mp_min, &
       mp_sendrecv, mp_sum, mp_wait, mp_waitall
  USE minimax,                         ONLY: check_range
  USE minimax_rpa,                     ONLY: get_rpa_minimax_coeff
  USE mp2_laplace,                     ONLY: laplace_minimax_approx
  USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                             mp2_type
  USE qs_environment_types,            ONLY: qs_environment_type
  USE rpa_communication,               ONLY: fm_redistribute,&
                                             initialize_buffer,&
                                             release_buffer
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_ri_gpw'

  PUBLIC :: rpa_ri_compute_en


  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param Erpa ...
!> \param mp2_env ...
!> \param BIb_C ...
!> \param qs_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param color_sub ...
!> \param ends_array ...
!> \param ends_B_virtual ...
!> \param sizes_array ...
!> \param sizes_B_virtual ...
!> \param starts_array ...
!> \param starts_B_virtual ...
!> \param Eigenval ...
!> \param nmo ...
!> \param homo ...
!> \param dimen_RI ...
!> \param unit_nr ...
!> \param error ...
!> \param error_sub ...
!> \param do_ri_sos_laplace_mp2 ...
!> \param BIb_C_beta ...
!> \param homo_beta ...
!> \param Eigenval_beta ...
!> \param ends_B_virtual_beta ...
!> \param sizes_B_virtual_beta ...
!> \param starts_B_virtual_beta ...
! *****************************************************************************
  SUBROUTINE rpa_ri_compute_en(Erpa,mp2_env,BIb_C,qs_env,para_env,para_env_sub,color_sub,&
                               ends_array,ends_B_virtual,sizes_array,sizes_B_virtual,starts_array,starts_B_virtual,&
                               Eigenval,nmo,homo,dimen_RI,unit_nr,error,error_sub,do_ri_sos_laplace_mp2,&
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
    LOGICAL                                  :: do_ri_sos_laplace_mp2
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
    INTEGER, PARAMETER                       :: occ_chunk = 128

    INTEGER :: best_integ_group_size, best_num_integ_point, color_rpa_group, &
      dimen_ia, dimen_ia_beta, handle, handle2, i, ia_global, ierr, iiB, &
      input_integ_group_size, integ_group_size, iproc, itmp(2), jjB, &
      min_integ_group_size, my_B_size, my_B_size_beta, my_B_virtual_end, &
      my_B_virtual_end_beta, my_B_virtual_start, my_B_virtual_start_beta, &
      my_group_L_end, my_group_L_size, my_group_L_start, my_homo_beta, &
      my_ia_end, my_ia_end_beta, my_ia_size, my_ia_size_beta, my_ia_start, &
      my_ia_start_beta, ngroup, num_integ_group, num_integ_points, occ_high, &
      occ_low, pos_integ_group, proc_receive, proc_send
    INTEGER :: proc_shift, rec_B_size, rec_B_virtual_end, &
      rec_B_virtual_start, stat, virtual, virtual_beta
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ends_ia, ends_ia_beta, &
                                                sizes_ia, sizes_ia_beta, &
                                                starts_ia, starts_ia_beta, &
                                                sub_proc_map
    LOGICAL                                  :: do_minimax_quad, failure, &
                                                my_open_shell, &
                                                skip_integ_group_opt
    REAL(KIND=dp) :: allowed_memory, avail_mem, E_Range, Emax, Emax_beta, &
      Emin, Emin_beta, mem_for_iaK, mem_for_QK, mem_min, mem_per_group, &
      mem_real, needed_mem
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: BIb_C_2D, BIb_C_2D_beta
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: BIb_C_rec
    TYPE(cp_fm_type), POINTER                :: fm_mat_Q, fm_mat_Q_beta, &
                                                fm_mat_Q_gemm, &
                                                fm_mat_Q_gemm_beta, fm_mat_S, &
                                                fm_mat_S_beta
    TYPE(cp_para_env_type), POINTER          :: para_env_RPA

    CALL timeset(routineN,handle)
    failure=.FALSE.

    my_open_shell=.FALSE.
    IF(PRESENT(BIb_C_beta).AND.&
       PRESENT(ends_B_virtual_beta).AND.&
       PRESENT(sizes_B_virtual_beta).AND.&
       PRESENT(starts_B_virtual_beta).AND.&
       PRESENT(homo_beta).AND.&
       PRESENT(Eigenval_beta)) my_open_shell=.TRUE.

    virtual=nmo-homo

    IF(do_ri_sos_laplace_mp2) THEN
      num_integ_points=mp2_env%ri_laplace%n_quadrature
      input_integ_group_size=mp2_env%ri_laplace%integ_group_size

      ! check the range for the minimax approximation
      Emin=2.0_dp*(Eigenval(homo+1)-Eigenval(homo))
      Emax=2.0_dp*(MAXVAL(Eigenval)-MINVAL(Eigenval))
      IF(my_open_shell.AND.homo_beta>0) THEN
        Emin_beta=2.0_dp*(Eigenval_beta(homo_beta+1)-Eigenval_beta(homo_beta))
        Emax_beta=2.0_dp*(MAXVAL(Eigenval_beta)-MINVAL(Eigenval_beta))
        Emin=MIN(Emin,Emin_beta)
        Emax=MAX(Emax,Emax_beta)
      END IF
      E_Range=Emax/Emin
      IF(E_Range<2.0_dp) E_Range=2.0_dp
      ierr=0
      CALL check_range(num_integ_points,E_Range,ierr)
      IF(ierr/=0) THEN
        jjB=num_integ_points-1
        DO iiB=1, jjB
          num_integ_points=num_integ_points-1
          ierr=0
          CALL check_range(num_integ_points,E_Range,ierr)
          IF(ierr==0) EXIT
        END DO
      END IF
      CPPostcondition(num_integ_points>=1,cp_failure_level,routineP,error_sub,failure)
    ELSE
      num_integ_points=mp2_env%ri_rpa%rpa_num_quad_points
      input_integ_group_size=mp2_env%ri_rpa%rpa_integ_group_size
      do_minimax_quad=mp2_env%ri_rpa%minimax_quad
      IF(do_minimax_quad.AND.num_integ_points>20) THEN
        CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                       "The required number of quadrature point exceeds the maximum possible in the "//&
                       "Minimax quadrature scheme. The number of quadrature point has been reset to 20.",&
                       only_ionode=.TRUE.)
        num_integ_points=20
      END IF
    END IF
    allowed_memory=mp2_env%mp2_memory

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

    mem_min=2.0_dp*REAL(homo,KIND=dp)*MAXVAL(sizes_B_virtual)*MAXVAL(sizes_array)*8.0_dp/(1024**2)
 
    IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T69,F9.2,A3)') 'RI_INFO| Minimum required memory per MPI process:',&
                                                        mem_min, ' MB'
   
    mem_real=allowed_memory-mem_real
    mem_real=MAX(mem_real,mem_min)

    IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T69,F9.2,A3)') 'RI_INFO| Available memory per MPI process:',&
                                                        mem_real, ' MB'
    
    mem_per_group=mem_real*para_env_sub%num_pe

    needed_mem=mem_for_iaK*2.0_dp+mem_for_QK*3.0_dp

    ! here we try to find the best rpa/lapleca group size
    skip_integ_group_opt=.FALSE.
     
    IF(input_integ_group_size>0) THEN
      IF(MOD(input_integ_group_size,para_env_sub%num_pe)==0) THEN
        best_integ_group_size=input_integ_group_size/para_env_sub%num_pe
        IF(MOD(ngroup,best_integ_group_size)==0) THEN
          num_integ_group=ngroup/best_integ_group_size
          IF((num_integ_points>num_integ_group).AND.MOD(num_integ_points,num_integ_group)==0) THEN
            best_num_integ_point=num_integ_points/num_integ_group
            skip_integ_group_opt=.TRUE.
          ELSE 
            IF (unit_nr>0) WRITE(unit_nr,'(T3,A)') 'NUM_QUAD_POINTS not multiple of the number of INTEG_GROUP'
          END IF
        ELSE
          IF (unit_nr>0) WRITE(unit_nr,'(T3,A)') 'Total number of groups not multiple of SIZE_INTEG_GROUP'
        END IF
      ELSE
        IF (unit_nr>0) WRITE(unit_nr,'(T3,A)') 'SIZE_INTEG_GROUP not multiple of GROUP_SIZE'
      END IF
    END IF

    IF(.NOT.skip_integ_group_opt) THEN
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
    END IF

    integ_group_size=best_integ_group_size

    IF (unit_nr>0) THEN
      IF(do_ri_sos_laplace_mp2) THEN 
         WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                                 "RI_INFO| Group size for laplace numerical integration:", integ_group_size*para_env_sub%num_pe
         WRITE (UNIT=unit_nr,FMT="(T3,A)")&
                                 "INTEG_INFO| MINIMAX approximation"
         WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                                 "INTEG_INFO| Number of integration points:", num_integ_points
         WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                                 "INTEG_INFO| Number of integration points per Laplace group:", best_num_integ_point
      ELSE
         WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                                 "RI_INFO| Group size for frequency integration:", integ_group_size*para_env_sub%num_pe
         IF(do_minimax_quad) THEN
           WRITE (UNIT=unit_nr,FMT="(T3,A)")&
                                   "INTEG_INFO| MINIMAX quadrature"
         ELSE
           WRITE (UNIT=unit_nr,FMT="(T3,A)")&
                                   "INTEG_INFO| Clenshaw-Curtius quadrature"
         END IF
         WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                                 "INTEG_INFO| Number of integration points:", num_integ_points
         WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                                 "INTEG_INFO| Number of integration points per RPA group:", best_num_integ_point
      END IF
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

    CALL timeset(routineN//"_reorder",handle2)

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

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,ia_global) &
    !$OMP          SHARED(homo,my_B_size,virtual,my_B_virtual_start,my_ia_start,my_ia_end,BIb_C,BIb_C_2D,&
    !$OMP          my_group_L_size)
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

      ! do this in chunks to avoid high memory overhead  for both BIb_C_rec and buffers in mp_sendrecv 
      ! TODO: fix this more cleanly with a rewrite sending only needed data etc.
      ! TODO: occ_chunk should presumably be precomputed so that messages are limited to e.g. 100Mb.
      ALLOCATE(BIb_C_rec(my_group_L_size,rec_B_size,MIN(homo,occ_chunk)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)

      DO occ_low=1,homo,occ_chunk
         occ_high=MIN(homo,occ_low+occ_chunk-1)
         CALL  mp_sendrecv(BIb_C(:,:,occ_low:occ_high),proc_send,&
                           BIb_C_rec(:,:,1:occ_high-occ_low+1),proc_receive,&
                           para_env_sub%group)
         !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,ia_global) &
         !$OMP          SHARED(occ_low,occ_high,rec_B_size,virtual,rec_B_virtual_start,my_ia_start,my_ia_end,BIb_C_rec,BIb_C_2D,&
         !$OMP                 my_group_L_size)
         DO iiB=occ_low, occ_high
           DO jjB=1, rec_B_size
             ia_global=(iiB-1)*virtual+rec_B_virtual_start+jjB-1
             IF(ia_global>=my_ia_start.AND.ia_global<=my_ia_end) THEN
               BIb_C_2D(ia_global-my_ia_start+1,1:my_group_L_size)=BIb_C_rec(1:my_group_L_size,jjB,iiB-occ_low+1)
             END IF
           END DO
         END DO
      ENDDO
      
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

     !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,ia_global) &
     !$OMP          SHARED(homo_beta,my_B_size_beta,virtual_beta,my_B_virtual_start_beta,&
     !$OMP                 my_ia_start_beta,my_ia_end_beta,BIb_C_beta,BIb_C_2D_beta,my_group_L_size)
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

        ALLOCATE(BIb_C_rec(my_group_L_size,rec_B_size,MIN(homo_beta,occ_chunk)),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)

        DO occ_low=1,homo_beta,occ_chunk
          occ_high=MIN(homo_beta,occ_low+occ_chunk-1)
          CALL  mp_sendrecv(BIb_C_beta(:,:,occ_low:occ_high),proc_send,&
                            BIb_C_rec(:,:,1:occ_high-occ_low+1),proc_receive,&
                            para_env_sub%group)

         !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,ia_global) &
         !$OMP          SHARED(occ_low,occ_high,rec_B_size,virtual_beta,rec_B_virtual_start,&
         !$OMP                 my_ia_start_beta,my_ia_end_beta,BIb_C_rec,BIb_C_2D_beta,my_group_L_size)
          DO iiB=occ_low, occ_high
            DO jjB=1, rec_B_size
              ia_global=(iiB-1)*virtual_beta+rec_B_virtual_start+jjB-1
              IF(ia_global>=my_ia_start_beta.AND.ia_global<=my_ia_end_beta) THEN
                BIb_C_2D_beta(ia_global-my_ia_start_beta+1,1:my_group_L_size)=BIb_C_rec(1:my_group_L_size,jjB,iiB-occ_low+1)
              END IF
            END DO
          END DO
        ENDDO

        DEALLOCATE(BIb_C_rec)
      END DO

      DEALLOCATE(BIb_C_beta)
      DEALLOCATE(starts_B_virtual_beta)
      DEALLOCATE(ends_B_virtual_beta)
      DEALLOCATE(sizes_B_virtual_beta)
    END IF
    CALL timestop(handle2)

    ! now create the matrices needed for the calculation, Q, S and G
    ! Q and G will have omega dependence
    CALL create_integ_mat(BIb_C_2D,para_env,para_env_sub,color_sub,ngroup,integ_group_size,num_integ_group,&
                          dimen_RI,homo,virtual,dimen_ia,color_rpa_group,&
                          mp2_env%block_size_row,mp2_env%block_size_col,unit_nr,&
                          sizes_ia,starts_ia,ends_ia,my_ia_size,my_ia_start,my_ia_end,&
                          sizes_array,starts_array,ends_array,my_group_L_size,my_group_L_start,my_group_L_end,&
                          para_env_RPA,fm_mat_S,fm_mat_Q_gemm,fm_mat_Q,error,error_sub)

    IF(my_open_shell) THEN
      ! the same for beta
      CALL create_integ_mat(BIb_C_2D_beta,para_env,para_env_sub,color_sub,ngroup,integ_group_size,num_integ_group,&
                            dimen_RI,homo_beta,virtual_beta,dimen_ia_beta,color_rpa_group,&
                            mp2_env%block_size_row,mp2_env%block_size_col,unit_nr,&
                            sizes_ia_beta,starts_ia_beta,ends_ia_beta,my_ia_size_beta,my_ia_start_beta,my_ia_end_beta,&
                            sizes_array,starts_array,ends_array,my_group_L_size,my_group_L_start,my_group_L_end,&
                            para_env_RPA,fm_mat_S_beta,fm_mat_Q_gemm_beta,fm_mat_Q_beta,error,error_sub,&
                            .TRUE.,fm_mat_Q%matrix_struct%context)

      IF(do_ri_sos_laplace_mp2) THEN
        ! go with laplace MINIMAX MP2
        CALL laplace_minimax_approx(Erpa,para_env,para_env_RPA,unit_nr,homo,virtual,dimen_RI,dimen_ia,Eigenval,&
                                    num_integ_points,num_integ_group,color_rpa_group,&
                                    fm_mat_S,fm_mat_Q_gemm,fm_mat_Q,error,&
                                    homo_beta,virtual_beta,dimen_ia_beta,Eigenval_beta,fm_mat_S_beta,&
                                    fm_mat_Q_gemm_beta,fm_mat_Q_beta)
      ELSE
        ! go with clenshaw-curtius/minimax quadrature
        CALL rpa_numerical_integ(Erpa,para_env,para_env_RPA,unit_nr,homo,virtual,dimen_RI,dimen_ia,Eigenval,&
                                 num_integ_points,num_integ_group,color_rpa_group,&
                                 fm_mat_S,fm_mat_Q_gemm,fm_mat_Q,mp2_env%ri_rpa%mm_style,do_minimax_quad,error,&
                                 homo_beta,virtual_beta,dimen_ia_beta,Eigenval_beta,fm_mat_S_beta,&
                                 fm_mat_Q_gemm_beta,fm_mat_Q_beta)
      END IF
    ELSE
      IF(do_ri_sos_laplace_mp2) THEN
        ! go with laplace MINIMAX MP2
        CALL laplace_minimax_approx(Erpa,para_env,para_env_RPA,unit_nr,homo,virtual,dimen_RI,dimen_ia,Eigenval,&
                                    num_integ_points,num_integ_group,color_rpa_group,&
                                    fm_mat_S,fm_mat_Q_gemm,fm_mat_Q,error)
      ELSE
        ! go with clenshaw-curtius/minimax quadrature
        CALL rpa_numerical_integ(Erpa,para_env,para_env_RPA,unit_nr,homo,virtual,dimen_RI,dimen_ia,Eigenval,&
                                 num_integ_points,num_integ_group,color_rpa_group,&
                                 fm_mat_S,fm_mat_Q_gemm,fm_mat_Q,mp2_env%ri_rpa%mm_style,do_minimax_quad,error)
      END IF
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
    !XXX CALL cp_fm_release(fm_mat_G,error=error)
    CALL cp_fm_release(fm_mat_Q_gemm,error=error)
    CALL cp_fm_release(fm_mat_Q,error=error)
    IF(my_open_shell) THEN
      DEALLOCATE(sizes_ia_beta)
      DEALLOCATE(starts_ia_beta)
      DEALLOCATE(ends_ia_beta)
      CALL cp_fm_release(fm_mat_S_beta,error=error)
      !XXX CALL cp_fm_release(fm_mat_G_beta,error=error)
      CALL cp_fm_release(fm_mat_Q_gemm_beta,error=error)
      CALL cp_fm_release(fm_mat_Q_beta,error=error)
    END IF

    CALL timestop(handle)

    END SUBROUTINE rpa_ri_compute_en

! *****************************************************************************
!> \brief ...
!> \param BIb_C_2D ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param color_sub ...
!> \param ngroup ...
!> \param integ_group_size ...
!> \param num_integ_group ...
!> \param dimen_RI ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_ia ...
!> \param color_rpa_group ...
!> \param ext_row_block_size ...
!> \param ext_col_block_size ...
!> \param unit_nr ...
!> \param sizes_ia ...
!> \param starts_ia ...
!> \param ends_ia ...
!> \param my_ia_size ...
!> \param my_ia_start ...
!> \param my_ia_end ...
!> \param sizes_array ...
!> \param starts_array ...
!> \param ends_array ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param para_env_RPA ...
!> \param fm_mat_S ...
!> \param fm_mat_Q_gemm ...
!> \param fm_mat_Q ...
!> \param error ...
!> \param error_sub ...
!> \param beta_case ...
!> \param blacs_env_ext ...
! *****************************************************************************
    SUBROUTINE create_integ_mat(BIb_C_2D,para_env,para_env_sub,color_sub,ngroup,integ_group_size,num_integ_group,&
                                dimen_RI,homo,virtual,dimen_ia,color_rpa_group,&
                                ext_row_block_size,ext_col_block_size,unit_nr,&
                                sizes_ia,starts_ia,ends_ia,my_ia_size,my_ia_start,my_ia_end,&
                                sizes_array,starts_array,ends_array,my_group_L_size,my_group_L_start,my_group_L_end,&
                                para_env_RPA,fm_mat_S,fm_mat_Q_gemm,fm_mat_Q,error,error_sub,beta_case,&
                                blacs_env_ext)

    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: BIb_C_2D
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_sub
    INTEGER :: color_sub, ngroup, integ_group_size, num_integ_group, &
      dimen_RI, homo, virtual, dimen_ia, color_rpa_group, ext_row_block_size, &
      ext_col_block_size, unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sizes_ia, starts_ia, ends_ia
    INTEGER                                  :: my_ia_size, my_ia_start, &
                                                my_ia_end
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sizes_array, starts_array, &
                                                ends_array
    INTEGER                                  :: my_group_L_size, &
                                                my_group_L_start, &
                                                my_group_L_end
    TYPE(cp_para_env_type), POINTER          :: para_env_RPA
    TYPE(cp_fm_type), POINTER                :: fm_mat_S, fm_mat_Q_gemm, &
                                                fm_mat_Q
    TYPE(cp_error_type), INTENT(inout)       :: error, error_sub
    LOGICAL, OPTIONAL                        :: beta_case
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env_ext

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_integ_mat', &
      routineP = moduleN//':'//routineN

    INTEGER :: comm_exchange, comm_rpa, dummy_proc, end_col_block, &
      end_row_block, grid_2D(2), handle, handle2, handle3, i, i_global, &
      i_local, i_sub, iiB, iii, iproc, iproc_col, iproc_row, itmp(2), &
      j_global, j_local, j_sub, jjB, mepos_in_RPA_group, my_num_col_blocks, &
      my_num_row_blocks, mypcol, myprow, ncol_block, ncol_block_mat, &
      ncol_local, npcol, nprow, nrow_block, nrow_block_mat, nrow_local, &
      num_rec_cols, number_of_rec, number_of_send, proc_receive, &
      proc_receive_static, proc_send, proc_send_static, proc_shift, &
      rec_counter, rec_ia_end, rec_ia_size, rec_ia_start, rec_L_end, &
      rec_L_size
    INTEGER :: rec_L_start, rec_pcol, rec_prow, ref_send_pcol, ref_send_prow, &
      row_col_proc_ratio, send_counter, send_pcol, send_prow, &
      size_rec_buffer, size_send_buffer, start_col_block, start_row_block, &
      stat, sub_sub_color
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: iii_vet, index_col_rec, &
                                                map_rec_size, map_send_size, &
                                                req_send, RPA_proc_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :) :: blocks_ranges_col, &
      blocks_ranges_row, grid_2_mepos, grid_ref_2_send_pos, &
      group_grid_2_mepos, mepos_2_grid, RPA_info_end, RPA_info_size, &
      RPA_info_start
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure, my_beta_case, &
                                                my_blacs_ext
    REAL(KIND=dp)                            :: part_ia, part_RI
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: BIb_C_rec
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env, blacs_env_Q
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct
    TYPE(cp_para_env_type), POINTER          :: para_env_exchange
    TYPE(integ_mat_buffer_type), &
      ALLOCATABLE, DIMENSION(:)              :: buffer_rec, buffer_send

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
    IF(para_env_RPA%num_pe>1) THEN
      row_col_proc_ratio=dimen_ia/dimen_RI
      row_col_proc_ratio=MAX(1,row_col_proc_ratio)
      IF(.FALSE.) THEN
        iproc_row=MAX(INT(SQRT(REAL(para_env_RPA%num_pe*row_col_proc_ratio,KIND=dp))),1)-1
        DO iproc=1, para_env_RPA%num_pe
          iproc_row=iproc_row+1
          IF(MOD(para_env_RPA%num_pe,iproc_row)==0) EXIT
        END DO
      ELSE
        iproc_row=MIN(MAX(INT(SQRT(REAL(para_env_RPA%num_pe*row_col_proc_ratio,KIND=dp))),1),para_env_RPA%num_pe)+1
        DO iproc=1, para_env_RPA%num_pe
          iproc_row=iproc_row-1
          IF(MOD(para_env_RPA%num_pe,iproc_row)==0) EXIT
        END DO
      END IF
      iproc_col=para_env_RPA%num_pe/iproc_row
      grid_2D(1)=iproc_row
      grid_2D(2)=iproc_col
    ELSE
      grid_2D=1
    END IF
    CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env_RPA, grid_2d=grid_2d, error=error)
    IF (unit_nr>0) THEN
      WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                                 "MATRIX_INFO| Number row processes:", grid_2D(1)
      WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                                 "MATRIX_INFO| Number column processes:", grid_2D(2)
    END IF
    ! CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env_RPA, error=error)

    ! define the block_size for the row
    IF(ext_row_block_size>0) THEN 
      nrow_block_mat=ext_row_block_size
    ELSE
      nrow_block_mat=dimen_ia/grid_2D(1)/2
      nrow_block_mat=MAX(nrow_block_mat,1)
    END IF

    ! define the block_size for the column
    IF(ext_col_block_size>0) THEN
      ncol_block_mat=ext_col_block_size
    ELSE
      ncol_block_mat=dimen_RI/grid_2D(2)/2
      ncol_block_mat=MAX(ncol_block_mat,1)
    END IF

    IF (unit_nr>0) THEN
      WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                                 "MATRIX_INFO| Row block size:", nrow_block_mat
      WRITE (UNIT=unit_nr,FMT="(T3,A,T75,i6)")&
                                 "MATRIX_INFO| Column block size:", ncol_block_mat
    END IF

    ! create the S full matrix, that is the (ia|K) matrix with K colomns
    ! and homo*virtual rows
    NULLIFY(fm_mat_S)
    NULLIFY(fm_struct)
    ! CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=dimen_ia,&
    !                          ncol_global=dimen_RI,para_env=para_env_RPA,error=error)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=dimen_ia,&
                             ncol_global=dimen_RI,para_env=para_env_RPA,&
                             nrow_block=nrow_block_mat, ncol_block=ncol_block_mat, force_block=.TRUE.,error=error)
    CALL cp_fm_create(fm_mat_S,fm_struct,name="fm_mat_S",error=error)

    CALL cp_fm_set_all(matrix=fm_mat_S,alpha=0.0_dp,error=error)

    ! fill the matrix
    CALL cp_fm_get_info(matrix=fm_mat_S,&
                        nrow_local=nrow_local,&
                        ncol_local=ncol_local,&
                        row_indices=row_indices,&
                        col_indices=col_indices,&
                        nrow_block=nrow_block,&
                        ncol_block=ncol_block,&
                        error=error)
    myprow=fm_mat_S%matrix_struct%context%mepos(1)
    mypcol=fm_mat_S%matrix_struct%context%mepos(2)
    nprow =fm_mat_S%matrix_struct%context%num_pe(1)
    npcol =fm_mat_S%matrix_struct%context%num_pe(2)

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

    IF(.TRUE.) THEN
      CALL timeset(routineN//"_comm",handle2)
      ! new replication scheme
      ! 0) create array with processes positions
      CALL timeset(routineN//"_info",handle3)
      ALLOCATE(grid_2_mepos(0:nprow-1,0:npcol-1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      grid_2_mepos=0
      ALLOCATE(mepos_2_grid(0:para_env_RPA%num_pe-1,2),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      mepos_2_grid=0

      grid_2_mepos(myprow,mypcol)=para_env_RPA%mepos
      mepos_2_grid(para_env_RPA%mepos,1)=myprow
      mepos_2_grid(para_env_RPA%mepos,2)=mypcol
 
      CALL mp_sum(grid_2_mepos,para_env_RPA%group)
      CALL mp_sum(mepos_2_grid,para_env_RPA%group)

      ! 1) loop over my local data and define a map for the proc to send data
      ALLOCATE(map_send_size(0:para_env_RPA%num_pe-1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      map_send_size=0
      DO jjB=my_group_L_start, my_group_L_end
        send_pcol=cp_fm_indxg2p(jjB,ncol_block,dummy_proc,&
                          fm_mat_S%matrix_struct%first_p_pos(2),npcol)
        DO iiB=my_ia_start, my_ia_end
          send_prow=cp_fm_indxg2p(iiB,nrow_block,dummy_proc,&
                            fm_mat_S%matrix_struct%first_p_pos(1),nprow)
          proc_send=grid_2_mepos(send_prow,send_pcol)
          map_send_size(proc_send)=map_send_size(proc_send)+1
        END DO
      END DO  

      ! 2) loop over my local data of fm_mat_S and define a map for the proc from which rec data
      ALLOCATE(map_rec_size(0:para_env_RPA%num_pe-1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      map_rec_size=0
      mepos_in_RPA_group=MOD(color_sub,integ_group_size)
      ALLOCATE(group_grid_2_mepos(0:para_env_sub%num_pe-1,0:integ_group_size-1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      group_grid_2_mepos=0
      group_grid_2_mepos(para_env_sub%mepos,mepos_in_RPA_group)=para_env_RPA%mepos
      CALL mp_sum(group_grid_2_mepos,para_env_RPA%group)
      part_ia=FLOAT(dimen_ia)/FLOAT(para_env_sub%num_pe)
      part_RI=FLOAT(dimen_RI)/FLOAT(ngroup)

      DO jjB=1, ncol_local
        j_global=col_indices(jjB) 
        ! check the group holding this element
        ! dirty way, if someone has a better idea ...
        rec_pcol=INT(FLOAT(j_global-1)/part_RI)
        rec_pcol=MAX(0,rec_pcol)
        rec_pcol=MIN(rec_pcol,ngroup-1)
        DO
          itmp=get_limit(dimen_RI,ngroup,rec_pcol)
          IF(j_global>=itmp(1).AND.j_global<=itmp(2)) EXIT
          IF(j_global<itmp(1)) rec_pcol=rec_pcol-1
          IF(j_global>itmp(2)) rec_pcol=rec_pcol+1
        END DO
        ! if the group is not in the same RPA group cycle
        IF((rec_pcol/integ_group_size).NE.color_rpa_group) CYCLE
        ! convert global position to position into the RPA group
        rec_pcol=MOD(rec_pcol,integ_group_size)

        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          ! check the process in the group holding this element
          rec_prow=INT(FLOAT(i_global-1)/part_ia)
          rec_prow=MAX(0,rec_prow)
          rec_prow=MIN(rec_prow,para_env_sub%num_pe-1)
          DO
            itmp=get_limit(dimen_ia,para_env_sub%num_pe,rec_prow)
            IF(i_global>=itmp(1).AND.i_global<=itmp(2)) EXIT
            IF(i_global<itmp(1)) rec_prow=rec_prow-1
            IF(i_global>itmp(2)) rec_prow=rec_prow+1
          END DO

          proc_receive=group_grid_2_mepos(rec_prow,rec_pcol)

          map_rec_size(proc_receive)=map_rec_size(proc_receive)+1

        END DO  ! i_global
      END DO  ! j_global

      ! 3) check if the local data has to be stored in the new fm_mat_S
      IF(map_rec_size(para_env_RPA%mepos)>0) THEN
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
      END IF
      CALL timestop(handle3)

      ! 4) reorder data in the send_buffer
      CALL timeset(routineN//"_buffer_s",handle3)
      number_of_send=0
      DO proc_shift=1, para_env_RPA%num_pe-1
        proc_send=RPA_proc_map(para_env_RPA%mepos+proc_shift)
        IF(map_send_size(proc_send)>0) THEN
          number_of_send=number_of_send+1
        END IF
      END DO

      ALLOCATE(buffer_send(number_of_send),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)

      ! this array given a pair (ref_send_prow,ref_send_pcol) returns
      ! the position in the buffer_send associated to that process
      ALLOCATE(grid_ref_2_send_pos(0:nprow-1,0:npcol-1),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      grid_ref_2_send_pos=0

      send_counter=0
      DO proc_shift=1, para_env_RPA%num_pe-1
        proc_send=RPA_proc_map(para_env_RPA%mepos+proc_shift)
        size_send_buffer=map_send_size(proc_send)

        IF(map_send_size(proc_send)>0) THEN
          send_counter=send_counter+1
          ! prepare the sending buffer
          ALLOCATE(buffer_send(send_counter)%msg(size_send_buffer),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
          buffer_send(send_counter)%msg=0.0_dp
          buffer_send(send_counter)%proc=proc_send

          ref_send_prow=mepos_2_grid(proc_send,1)
          ref_send_pcol=mepos_2_grid(proc_send,2)

          grid_ref_2_send_pos(ref_send_prow,ref_send_pcol)=send_counter

          ! iii=0
          ! DO iiB=my_ia_start, my_ia_end
          !   send_prow=cp_fm_indxg2p(iiB,nrow_block,dummy_proc,&
          !                           fm_mat_S%matrix_struct%first_p_pos(1),nprow)
          !   IF(send_prow.NE.ref_send_prow) CYCLE
          !   DO jjB=my_group_L_start, my_group_L_end
          !     send_pcol=cp_fm_indxg2p(jjB,ncol_block,dummy_proc,&
          !                             fm_mat_S%matrix_struct%first_p_pos(2),npcol)
          !     IF(send_pcol.NE.ref_send_pcol) CYCLE
          !     iii=iii+1
          !     buffer_send(send_counter)%msg(iii)=BIb_C_2D(iiB-my_ia_start+1,jjB-my_group_L_start+1)
          !   END DO
          ! END DO

        END IF

      END DO

      ! loop over the locally held data and fill the buffer_send
      ! for doing that we need an array index 
      ALLOCATE(iii_vet(number_of_send),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
    
      iii_vet=0
      DO iiB=my_ia_start, my_ia_end
        send_prow=cp_fm_indxg2p(iiB,nrow_block,dummy_proc,&
                                fm_mat_S%matrix_struct%first_p_pos(1),nprow)
        DO jjB=my_group_L_start, my_group_L_end
          send_pcol=cp_fm_indxg2p(jjB,ncol_block,dummy_proc,&
                                  fm_mat_S%matrix_struct%first_p_pos(2),npcol)
          ! we don't need to send to ourselves
          IF(grid_2_mepos(send_prow,send_pcol)==para_env_RPA%mepos) CYCLE

          send_counter=grid_ref_2_send_pos(send_prow,send_pcol)
          iii_vet(send_counter)=iii_vet(send_counter)+1
          iii=iii_vet(send_counter)
          buffer_send(send_counter)%msg(iii)=BIb_C_2D(iiB-my_ia_start+1,jjB-my_group_L_start+1)
        END DO
      END DO

      DEALLOCATE(iii_vet)
      DEALLOCATE(grid_ref_2_send_pos)
      DEALLOCATE(BIb_C_2D)
      CALL timestop(handle3)

      ! 5) create the buffer for receive and post the irecv
      CALL timeset(routineN//"_isendrecv",handle3)
      number_of_rec=0
      DO proc_shift=1, para_env_RPA%num_pe-1
        proc_receive=RPA_proc_map(para_env_RPA%mepos-proc_shift)
        IF(map_rec_size(proc_receive)>0) THEN
          number_of_rec=number_of_rec+1
        END IF
      END DO
  
      ALLOCATE(buffer_rec(number_of_rec),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)

      rec_counter=0
      DO proc_shift=1, para_env_RPA%num_pe-1
        proc_receive=RPA_proc_map(para_env_RPA%mepos-proc_shift)
        size_rec_buffer=map_rec_size(proc_receive)

        IF(map_rec_size(proc_receive)>0) THEN
          rec_counter=rec_counter+1
          ! prepare the buffer for receive
          ALLOCATE(buffer_rec(rec_counter)%msg(size_rec_buffer),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
          buffer_rec(rec_counter)%msg=0.0_dp
          buffer_rec(rec_counter)%proc=proc_receive

          ! post the receiving message
          CALL mp_irecv(buffer_rec(rec_counter)%msg,proc_receive,para_env_RPA%group,buffer_rec(rec_counter)%msg_req)
        END IF
      END DO

      ! 6) send data
      ALLOCATE(req_send(number_of_send),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      send_counter=0
      DO proc_shift=1, para_env_RPA%num_pe-1
        proc_send=RPA_proc_map(para_env_RPA%mepos+proc_shift)
        IF(map_send_size(proc_send)>0) THEN
          send_counter=send_counter+1
          CALL mp_isend(buffer_send(send_counter)%msg,proc_send,para_env_RPA%group,buffer_send(send_counter)%msg_req)
          req_send(send_counter)=buffer_send(send_counter)%msg_req
        END IF
      END DO
      CALL timestop(handle3)

      ! 8) fill the fm_mat_S matrix
      CALL timeset(routineN//"_fill",handle3)
      ! In order to perform this step efficiently first we have to know the  
      ! ranges of the blocks over which a given process hold its local data.
      ! Start with the rows ...
      my_num_row_blocks=1
      DO iiB=1, nrow_local-1
        IF(ABS(row_indices(iiB+1)-row_indices(iiB))==1) CYCLE
        my_num_row_blocks=my_num_row_blocks+1
      END DO
      ALLOCATE(blocks_ranges_row(2,my_num_row_blocks),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      blocks_ranges_row=0
      blocks_ranges_row(1,1)=row_indices(1)
      iii=1
      DO iiB=1, nrow_local-1
        IF(ABS(row_indices(iiB+1)-row_indices(iiB))==1) CYCLE
        iii=iii+1
        blocks_ranges_row(2,iii-1)=row_indices(iiB)
        blocks_ranges_row(1,iii)=row_indices(iiB+1)
      END DO
      blocks_ranges_row(2,my_num_row_blocks)=row_indices(MAX(nrow_local,1))

      ! and columns
      my_num_col_blocks=1
      DO jjB=1, ncol_local-1
        IF(ABS(col_indices(jjB+1)-col_indices(jjB))==1) CYCLE
        my_num_col_blocks=my_num_col_blocks+1
      END DO
      ALLOCATE(blocks_ranges_col(2,my_num_col_blocks),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
      blocks_ranges_col=0
      blocks_ranges_col(1,1)=col_indices(1)
      iii=1
      DO jjB=1, ncol_local-1
        IF(ABS(col_indices(jjB+1)-col_indices(jjB))==1) CYCLE
        iii=iii+1
        blocks_ranges_col(2,iii-1)=col_indices(jjB)
        blocks_ranges_col(1,iii)=col_indices(jjB+1)
      END DO
      blocks_ranges_col(2,my_num_col_blocks)=col_indices(MAX(ncol_local,1))

      rec_counter=0
      DO proc_shift=1, para_env_RPA%num_pe-1
        proc_receive=RPA_proc_map(para_env_RPA%mepos-proc_shift)
        size_rec_buffer=map_rec_size(proc_receive)

        IF(map_rec_size(proc_receive)>0) THEN
          rec_counter=rec_counter+1

          rec_L_size=RPA_info_size(1,proc_receive)
          rec_L_start=RPA_info_start(1,proc_receive)
          rec_L_end=RPA_info_end(1,proc_receive)
          ! precompute the number of received columns and relative index
          num_rec_cols=0
          DO jjB=1, my_num_col_blocks
            start_col_block=MAX(blocks_ranges_col(1,jjB),rec_L_start)
            end_col_block=MIN(blocks_ranges_col(2,jjB),rec_L_end)
            DO j_sub=start_col_block, end_col_block
              num_rec_cols=num_rec_cols+1
            END DO
          END DO
          ALLOCATE(index_col_rec(num_rec_cols),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error_sub,failure)
          index_col_rec=0
          iii=0
          DO jjB=1, my_num_col_blocks
            start_col_block=MAX(blocks_ranges_col(1,jjB),rec_L_start)
            end_col_block=MIN(blocks_ranges_col(2,jjB),rec_L_end)
            DO j_sub=start_col_block, end_col_block
              iii=iii+1
              j_local=cp_fm_indxg2l(j_sub,ncol_block,dummy_proc,&
                                    fm_mat_S%matrix_struct%first_p_pos(2),npcol)
              index_col_rec(iii)=j_local
            END DO
          END DO

          rec_ia_size=RPA_info_size(2,proc_receive)
          rec_ia_start=RPA_info_start(2,proc_receive)
          rec_ia_end=RPA_info_end(2,proc_receive)

          ! wait for the message 
          CALL mp_wait(buffer_rec(rec_counter)%msg_req)

          iii=0
          DO iiB=1, my_num_row_blocks
            start_row_block=MAX(blocks_ranges_row(1,iiB),rec_ia_start)
            end_row_block=MIN(blocks_ranges_row(2,iiB),rec_ia_end)
            DO i_sub=start_row_block, end_row_block
              i_local=cp_fm_indxg2l(i_sub,nrow_block,dummy_proc,&
                                    fm_mat_S%matrix_struct%first_p_pos(1),nprow)
              ! DO jjB=1, my_num_col_blocks
              !   start_col_block=MAX(blocks_ranges_col(1,jjB),rec_L_start)
              !   end_col_block=MIN(blocks_ranges_col(2,jjB),rec_L_end)
              !   DO j_sub=start_col_block, end_col_block
              !     j_local=cp_fm_indxg2l(j_sub,ncol_block,dummy_proc,&
              !                           fm_mat_S%matrix_struct%first_p_pos(2),npcol)
              !     iii=iii+1
              !     fm_mat_S%local_data(i_local,j_local)=buffer_rec(rec_counter)%msg(iii)
              !   END DO
              ! END DO
              DO jjB=1, num_rec_cols 
                iii=iii+1
                j_local=index_col_rec(jjB)
                fm_mat_S%local_data(i_local,j_local)=buffer_rec(rec_counter)%msg(iii)
              END DO
            END DO
          END DO

          ! iii=0
          ! DO iiB=rec_ia_start, rec_ia_end
          !   rec_prow=cp_fm_indxg2p(iiB,nrow_block,dummy_proc,&
          !                    fm_mat_S%matrix_struct%first_p_pos(1),nprow)
          !   IF(rec_prow.NE.myprow) CYCLE
          !   i_local=cp_fm_indxg2l(iiB,nrow_block,dummy_proc,&
          !                         fm_mat_S%matrix_struct%first_p_pos(1),nprow)
          !   DO jjB=rec_L_start, rec_L_end
          !     rec_pcol=cp_fm_indxg2p(jjB,ncol_block,dummy_proc,&
          !                      fm_mat_S%matrix_struct%first_p_pos(2),npcol)
          !     IF(rec_pcol.NE.mypcol) CYCLE
          !     j_local=cp_fm_indxg2l(jjB,ncol_block,dummy_proc,&
          !                     fm_mat_S%matrix_struct%first_p_pos(2),npcol)
          !     iii=iii+1
          !     fm_mat_S%local_data(i_local,j_local)=buffer_rec(rec_counter)%msg(iii)
          !   END DO
          ! END DO

          DEALLOCATE(buffer_rec(rec_counter)%msg) 
          DEALLOCATE(index_col_rec)
          
        END IF
      END DO
      DEALLOCATE(buffer_rec)
      
      DEALLOCATE(blocks_ranges_row)
      DEALLOCATE(blocks_ranges_col)
      
      CALL timestop(handle3)
      
      ! wait for all messeges to be sent
      CALL timeset(routineN//"_waitall",handle3)
      CALL mp_waitall(req_send(:))
      DO send_counter=1, number_of_send
        DEALLOCATE(buffer_send(send_counter)%msg)
      END DO
      DEALLOCATE(buffer_send)
      CALL timestop(handle3)

      DEALLOCATE(group_grid_2_mepos)
      DEALLOCATE(map_send_size)
      DEALLOCATE(map_rec_size)
      DEALLOCATE(grid_2_mepos)
      DEALLOCATE(mepos_2_grid)

      CALL timestop(handle2)
    ELSE   
      ! old inefficient replication scheme
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
        BIb_C_2D(:,:)=BIb_C_rec

        DEALLOCATE(BIb_C_rec)
      END DO

      DEALLOCATE(BIb_C_2D)
    END IF

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

    CALL timeset(routineN//"_sum",handle2)
    CALL mp_sum(fm_mat_S%local_data,para_env_exchange%group)
    CALL timestop(handle2)

    CALL cp_para_env_release(para_env_exchange,error=error_sub)

    ! create the twin matrix for the mat-mat-mul (mat_Q)
    !XXX CALL cp_fm_create(fm_mat_G,fm_struct,name="fm_mat_G",error=error)
    !XXX CALL cp_fm_set_all(matrix=fm_mat_G,alpha=0.0_dp,error=error)

    CALL cp_fm_struct_release(fm_struct,error=error)
    

    ! create the Q matrix dimen_RIxdimen_RI where the result of the mat-mat-mult will be stored
    NULLIFY(fm_mat_Q_gemm)
    NULLIFY(fm_struct)
    !CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=dimen_RI,&
    !                         ncol_global=dimen_RI,para_env=para_env_RPA,error=error)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=dimen_RI,&
                             ncol_global=dimen_RI,para_env=para_env_RPA,&
                             nrow_block=nrow_block_mat, ncol_block=ncol_block_mat, force_block=.TRUE.,error=error)
    CALL cp_fm_create(fm_mat_Q_gemm,fm_struct,name="fm_mat_Q_gemm",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL cp_fm_set_all(matrix=fm_mat_Q_gemm,alpha=0.0_dp,error=error)

    ! create the Q matrix with a different blacs env
    NULLIFY(blacs_env_Q)
    IF(my_blacs_ext) THEN 
      blacs_env_Q=>blacs_env_ext
    ELSE
      CALL cp_blacs_env_create(blacs_env=blacs_env_Q, para_env=para_env_RPA, error=error)
    END IF

    NULLIFY(fm_mat_Q)
    NULLIFY(fm_struct)
    CALL cp_fm_struct_create(fm_struct,context=blacs_env_Q,nrow_global=dimen_RI,&
                             ncol_global=dimen_RI,para_env=para_env_RPA,error=error)
    CALL cp_fm_create(fm_mat_Q,fm_struct,name="fm_mat_Q",error=error)
    CALL cp_fm_struct_release(fm_struct,error=error)

    CALL cp_fm_set_all(matrix=fm_mat_Q,alpha=0.0_dp,error=error)

    ! release blacs_env
    IF(.NOT.my_blacs_ext) CALL cp_blacs_env_release(blacs_env_Q,error=error)
    CALL cp_blacs_env_release(blacs_env,error=error)

    CALL timestop(handle)

    END SUBROUTINE create_integ_mat

! *****************************************************************************
!> \brief ...
!> \param Erpa ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param unit_nr ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param dimen_ia ...
!> \param Eigenval ...
!> \param num_integ_points ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param fm_mat_S ...
!> \param fm_mat_Q_gemm ...
!> \param fm_mat_Q ...
!> \param mm_style ...
!> \param do_minimax_quad ...
!> \param error ...
!> \param homo_beta ...
!> \param virtual_beta ...
!> \param dimen_ia_beta ...
!> \param Eigenval_beta ...
!> \param fm_mat_S_beta ...
!> \param fm_mat_Q_gemm_beta ...
!> \param fm_mat_Q_beta ...
! *****************************************************************************
    SUBROUTINE rpa_numerical_integ(Erpa,para_env,para_env_RPA,unit_nr,homo,virtual,dimen_RI,dimen_ia,Eigenval,&
                                   num_integ_points,num_integ_group,color_rpa_group,&
                                   fm_mat_S,fm_mat_Q_gemm,fm_mat_Q,mm_style,do_minimax_quad,error,&
                                   homo_beta,virtual_beta,dimen_ia_beta,Eigenval_beta,fm_mat_S_beta,&
                                   fm_mat_Q_gemm_beta,fm_mat_Q_beta)
    REAL(KIND=dp)                            :: Erpa
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_RPA
    INTEGER                                  :: unit_nr, homo, virtual, &
                                                dimen_RI, dimen_ia
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    INTEGER                                  :: num_integ_points, &
                                                num_integ_group, &
                                                color_rpa_group
    TYPE(cp_fm_type), POINTER                :: fm_mat_S, fm_mat_Q_gemm, &
                                                fm_mat_Q
    INTEGER                                  :: mm_style
    LOGICAL                                  :: do_minimax_quad
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, OPTIONAL                        :: homo_beta, virtual_beta, &
                                                dimen_ia_beta
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: Eigenval_beta
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: fm_mat_S_beta, &
                                                fm_mat_Q_gemm_beta, &
                                                fm_mat_Q_beta

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rpa_numerical_integ', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: use_gemm = .FALSE.

    INTEGER :: avirt, handle, i_global, ierr, iiB, info_chol, iocc, j_global, &
      jjB, jquad, my_num_dgemm_call, ncol_local, nrow_local, number_of_rec, &
      number_of_rec_beta, number_of_send, number_of_send_beta, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: map_rec_size, &
                                                map_rec_size_beta, &
                                                map_send_size, &
                                                map_send_size_beta, &
                                                RPA_proc_map
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: local_size_source, &
                                                local_size_source_beta
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    LOGICAL                                  :: failure, first_cycle, &
                                                my_open_shell
    REAL(KIND=dp) :: a_scaling, actual_flop_rate, alpha, E_Range, eigen_diff, &
      Emax, Emax_beta, Emin, Emin_beta, FComega, my_flop_rate, omega, &
      omega_old, t_end, t_start
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: Q_log, tj, trace_Qomega, wj, &
                                                x_tw
    TYPE(integ_mat_buffer_type), &
      ALLOCATABLE, DIMENSION(:)              :: buffer_rec, buffer_rec_beta, &
                                                buffer_send, buffer_send_beta

    CALL timeset(routineN,handle)
    failure=.FALSE.
    
    my_open_shell=.FALSE.
    IF(PRESENT(homo_beta).AND.&
       PRESENT(virtual_beta).AND.&
       PRESENT(dimen_ia_beta).AND.&
       PRESENT(Eigenval_beta).AND.&
       PRESENT(fm_mat_S_beta).AND.&
       PRESENT(fm_mat_Q_gemm_beta).AND.&
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

    IF(do_minimax_quad) THEN
      ! MINIMAX quadrature
      ALLOCATE(x_tw(2*num_integ_points),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      x_tw=0.0_dp

      Emin=Eigenval(homo+1)-Eigenval(homo)
      Emax=MAXVAL(Eigenval)-MINVAL(Eigenval)
      IF(my_open_shell.AND.homo_beta>0) THEN
        Emin_beta=Eigenval_beta(homo_beta+1)-Eigenval_beta(homo_beta)
        Emax_beta=MAXVAL(Eigenval_beta)-MINVAL(Eigenval_beta)
        Emin=MIN(Emin,Emin_beta)
        Emax=MAX(Emax,Emax_beta)
      END IF
      E_Range=Emax/Emin

      ierr=0
      CALL get_rpa_minimax_coeff(num_integ_points,E_Range,x_tw,ierr,error)

      DO jquad=1, num_integ_points
        tj(jquad)=x_tw(jquad)
        wj(jquad)=x_tw(jquad+num_integ_points)
      END DO

      DEALLOCATE(x_tw)

      IF (unit_nr>0) THEN
        WRITE (UNIT=unit_nr,FMT="(T3,A,T66,F15.4)")&
                            "INTEG_INFO| Range for the minimax approximation:", E_Range
        WRITE (UNIT=unit_nr,FMT="(T3,A,T54,A,T72,A)") "INTEG_INFO| Minimax parameters:","Weights","Abscissas"
        DO jquad=1, num_integ_points
          WRITE (UNIT=unit_nr,FMT="(T41,F20.10,F20.10)") wj(jquad), tj(jquad)
        END DO
        CALL m_flush(unit_nr)
      END IF

      ! scale the minimax parameters
      tj=tj*Emin
      wj=wj*Emin

    ELSE
      ! Clenshaw-Curtius quadrature
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
                                 tj,wj,fm_mat_S,error,&
                                 homo_beta,virtual_beta,dimen_ia_beta,Eigenval_beta,fm_mat_S_beta)
      ELSE
        CALL calc_scaling_factor(a_scaling,para_env,para_env_RPA,homo,virtual,dimen_RI,dimen_ia,Eigenval,&
                                 num_integ_points,num_integ_group,color_rpa_group,&
                                 tj,wj,fm_mat_S,error)
      END IF

      IF (unit_nr>0) WRITE(unit_nr,'(T3,A,T56,F25.5)') 'INTEG_INFO| Scaling parameter:', a_scaling

      wj=wj*a_scaling

    END IF

    ! initialize buffer for matrix redistribution
    CALL initialize_buffer(fm_mat_Q_gemm,fm_mat_Q,RPA_proc_map,buffer_rec,buffer_send,&
                           number_of_rec,number_of_send,&
                           map_send_size,map_rec_size,local_size_source,para_env_RPA,error)
    IF(my_open_shell) THEN
      CALL initialize_buffer(fm_mat_Q_gemm_beta,fm_mat_Q_beta,RPA_proc_map,buffer_rec_beta,buffer_send_beta,&
                             number_of_rec_beta,number_of_send_beta,&
                             map_send_size_beta,map_rec_size_beta,local_size_source_beta,para_env_RPA,error)
    END IF

    IF(my_open_shell) THEN
      alpha=2.0_dp
    ELSE
      alpha=4.0_dp
    END IF

    Erpa=0.0_dp
    first_cycle=.TRUE.
    omega_old=0.0_dp
    my_num_dgemm_call=0
    my_flop_rate=0.0_dp
    DO jquad=1, num_integ_points
      IF(MODULO(jquad,num_integ_group)/=color_rpa_group) CYCLE
      
      IF(do_minimax_quad) THEN
        omega=tj(jquad)
      ELSE
        omega=a_scaling/TAN(tj(jquad))
      END IF

      !XXX ! copy fm_mat_S into fm_mat_G
      !XXX CALL cp_fm_to_fm(source=fm_mat_S,destination=fm_mat_G,error=error)

      !XXX ! get info of fm_mat_G
      !XXX CALL cp_fm_get_info(matrix=fm_mat_G,&
      !XXX                     nrow_local=nrow_local,&
      !XXX                     ncol_local=ncol_local,&
      !XXX                     row_indices=row_indices,&
      !XXX                     col_indices=col_indices,&
      !XXX                     error=error)
      
      ! get info of fm_mat_S
      CALL cp_fm_get_info(matrix=fm_mat_S,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices,&
                          error=error)

      ! update G matrix with the new value of omega
      IF(first_cycle) THEN
        ! In this case just update the matrix (symmetric form) witrh 
        ! SQRT((epsi_a-epsi_i)/((epsi_a-epsi_i)**2+omega**2))
        !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
        !$OMP             SHARED(ncol_local,nrow_local,col_indices,row_indices,Eigenval,fm_mat_S,virtual,homo,omega)
        DO jjB=1, ncol_local
          j_global=col_indices(jjB)
          DO iiB=1, nrow_local
            i_global=row_indices(iiB)

            iocc=MAX(1,i_global-1)/virtual+1
            avirt=i_global-(iocc-1)*virtual
            eigen_diff=Eigenval(avirt+homo)-Eigenval(iocc)

            fm_mat_S%local_data(iiB,jjB)=fm_mat_S%local_data(iiB,jjB)*&
                                         SQRT(eigen_diff/(eigen_diff**2+omega**2))

          END DO
        END DO
      ELSE
        ! In this case the update has to remove the old omega component thus
        ! SQRT(((epsi_a-epsi_i)**2+omega_old**2)/((epsi_a-epsi_i)**2+omega**2))
        !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
        !$OMP             SHARED(ncol_local,nrow_local,col_indices,row_indices,Eigenval,&
        !$OMP                    fm_mat_S,virtual,homo,omega,omega_old)
        DO jjB=1, ncol_local
          j_global=col_indices(jjB)
          DO iiB=1, nrow_local
            i_global=row_indices(iiB)

            iocc=MAX(1,i_global-1)/virtual+1
            avirt=i_global-(iocc-1)*virtual
            eigen_diff=Eigenval(avirt+homo)-Eigenval(iocc)

            fm_mat_S%local_data(iiB,jjB)=fm_mat_S%local_data(iiB,jjB)*&
                                         SQRT((eigen_diff**2+omega_old**2)/(eigen_diff**2+omega**2))

          END DO
        END DO
      END IF

      ! alpha=4.0 is valid only for closed shell systems (alpha is 2.0 in the spin-orbital basis)
      t_start=m_walltime()
      SELECT CASE(mm_style)
      CASE(wfc_mm_style_gemm)
         ! waste-fully computes the full symmetrix matrix, but maybe faster than cp_fm_syrk for optimized cp_fm_gemm
         CALL cp_gemm(transa="T",transb="N",m=dimen_RI,n=dimen_RI,k=dimen_ia,alpha=alpha,&
                         matrix_a=fm_mat_S,matrix_b=fm_mat_S,beta=0.0_dp,&
                         matrix_c=fm_mat_Q_gemm,error=error)
      CASE(wfc_mm_style_syrk)
         ! will only compute the upper half of the matrix, which is fine, since we only use it for cholesky later
         CALL cp_fm_syrk(uplo='U',trans='T',k=dimen_ia,alpha=alpha,matrix_a=fm_mat_S,&
                         ia=1,ja=1,beta=0.0_dp,matrix_c=fm_mat_Q_gemm,error=error)
      CASE DEFAULT
         CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT
      t_end=m_walltime()

      actual_flop_rate=2.0_dp*REAL(dimen_ia,KIND=dp)*dimen_RI*REAL(dimen_RI,KIND=dp)
      actual_flop_rate=actual_flop_rate/(MAX(0.01_dp,t_end-t_start))
      IF(para_env_RPA%mepos==0) my_flop_rate=my_flop_rate+actual_flop_rate
      my_num_dgemm_call=my_num_dgemm_call+1

      ! copy/redistribute fm_mat_Q_gemm to fm_mat_Q
      CALL cp_fm_set_all(matrix=fm_mat_Q,alpha=0.0_dp,error=error)
      CALL fm_redistribute(fm_mat_Q_gemm,fm_mat_Q,RPA_proc_map,buffer_rec,buffer_send,&
                           number_of_rec,number_of_send,&
                           map_send_size,map_rec_size,local_size_source,para_env_RPA,error)

      IF(my_open_shell) THEN
        ! the same for the beta spin
        !XXX ! copy fm_mat_S into fm_mat_G
        !XXX CALL cp_fm_to_fm(source=fm_mat_S_beta,destination=fm_mat_G_beta,error=error)
        !XXX ! get info of fm_mat_G_beta
        !XXX CALL cp_fm_get_info(matrix=fm_mat_G_beta,&
        !XXX                     nrow_local=nrow_local,&
        !XXX                     ncol_local=ncol_local,&
        !XXX                     row_indices=row_indices,&
        !XXX                     col_indices=col_indices,&
        !XXX                     error=error)

        ! get info of fm_mat_S_beta
        CALL cp_fm_get_info(matrix=fm_mat_S_beta,&
                            nrow_local=nrow_local,&
                            ncol_local=ncol_local,&
                            row_indices=row_indices,&
                            col_indices=col_indices,&
                            error=error)

        ! update G matrix with the new value of omega
        IF(first_cycle) THEN
          !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
          !$OMP             SHARED(ncol_local,nrow_local,col_indices,row_indices,&
          !$OMP                    Eigenval_beta,fm_mat_S_beta,virtual_beta,homo_beta,omega)
          DO jjB=1, ncol_local
            j_global=col_indices(jjB)
            DO iiB=1, nrow_local
              i_global=row_indices(iiB)

              iocc=MAX(1,i_global-1)/virtual_beta+1
              avirt=i_global-(iocc-1)*virtual_beta
              eigen_diff=Eigenval_beta(avirt+homo_beta)-Eigenval_beta(iocc)

              fm_mat_S_beta%local_data(iiB,jjB)=fm_mat_S_beta%local_data(iiB,jjB)*&
                                                SQRT(eigen_diff/(eigen_diff**2+omega**2))

            END DO
          END DO
        ELSE
          !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
          !$OMP             SHARED(ncol_local,nrow_local,col_indices,row_indices,&
          !$OMP                    Eigenval_beta,fm_mat_S_beta,virtual_beta,homo_beta,omega,omega_old)
          DO jjB=1, ncol_local
            j_global=col_indices(jjB)
            DO iiB=1, nrow_local
              i_global=row_indices(iiB)

              iocc=MAX(1,i_global-1)/virtual_beta+1
              avirt=i_global-(iocc-1)*virtual_beta
              eigen_diff=Eigenval_beta(avirt+homo_beta)-Eigenval_beta(iocc)

              fm_mat_S_beta%local_data(iiB,jjB)=fm_mat_S_beta%local_data(iiB,jjB)*&
                                                SQRT((eigen_diff**2+omega_old**2)/(eigen_diff**2+omega**2))

            END DO
          END DO
        END IF

        t_start=m_walltime()
        SELECT CASE(mm_style)
        CASE(wfc_mm_style_gemm)
           CALL cp_gemm(transa="T",transb="N",m=dimen_RI,n=dimen_RI,k=dimen_ia_beta,alpha=alpha,&
                           matrix_a=fm_mat_S_beta,matrix_b=fm_mat_S_beta,beta=0.0_dp,&
                           matrix_c=fm_mat_Q_gemm_beta,error=error)
        CASE(wfc_mm_style_syrk)
           CALL cp_fm_syrk(uplo='U',trans='T',k=dimen_ia_beta,alpha=alpha,matrix_a=fm_mat_S_beta,&
                           ia=1,ja=1,beta=0.0_dp,matrix_c=fm_mat_Q_gemm_beta,error=error)
        CASE DEFAULT
           CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
        END SELECT
        t_end=m_walltime()
        actual_flop_rate=2.0_dp*REAL(dimen_ia_beta,KIND=dp)*dimen_RI*REAL(dimen_RI,KIND=dp)/(MAX(0.01_dp,t_end-t_start))
        IF(para_env_RPA%mepos==0) my_flop_rate=my_flop_rate+actual_flop_rate
        my_num_dgemm_call=my_num_dgemm_call+1

        ! copy/redistribute fm_mat_Q_gemm to fm_mat_Q
        CALL cp_fm_set_all(matrix=fm_mat_Q_beta,alpha=0.0_dp,error=error)
        CALL fm_redistribute(fm_mat_Q_gemm_beta,fm_mat_Q_beta,RPA_proc_map,buffer_rec_beta,buffer_send_beta,&
                             number_of_rec_beta,number_of_send_beta,&
                             map_send_size_beta,map_rec_size_beta,local_size_source_beta,para_env_RPA,error)

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
      !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,i_global,j_global) &
      !$OMP             SHARED(ncol_local,nrow_local,col_indices,row_indices,trace_Qomega,fm_mat_Q,dimen_RI)
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
      !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,i_global,j_global) &
      !$OMP             SHARED(ncol_local,nrow_local,col_indices,row_indices,Q_log,fm_mat_Q,dimen_RI)
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

      ! save omega and reset the first_cycle flag
      first_cycle=.FALSE.
      omega_old=omega
    END DO

    CALL mp_sum(Erpa,para_env%group)
    Erpa=Erpa/(pi*2.0_dp)
    IF(do_minimax_quad) Erpa=Erpa/2.0_dp

    IF(para_env_RPA%mepos==0) my_flop_rate=my_flop_rate/REAL(MAX(my_num_dgemm_call,1),KIND=dp)/1.0E9_dp
    CALL mp_sum(my_flop_rate,para_env%group)
    my_flop_rate=my_flop_rate/para_env%num_pe
    IF (unit_nr>0) WRITE (UNIT=unit_nr,FMT="(T3,A,T66,F15.2)")&
                          "PERFORMANCE| PDGEMM flop rate (Gflops / MPI rank):", my_flop_rate

    ! release buffer
    CALL release_buffer(RPA_proc_map,buffer_rec,buffer_send,&
                        number_of_rec,number_of_send,&
                        map_send_size,map_rec_size,local_size_source,error)

    IF(my_open_shell) THEN
      CALL release_buffer(RPA_proc_map,buffer_rec_beta,buffer_send_beta,&
                          number_of_rec_beta,number_of_send_beta,&
                          map_send_size_beta,map_rec_size_beta,local_size_source_beta,error)
    END IF

    DEALLOCATE(tj)
    DEALLOCATE(wj)
    DEALLOCATE(Q_log)
    DEALLOCATE(trace_Qomega)

    CALL timestop(handle)

    END SUBROUTINE rpa_numerical_integ

! *****************************************************************************
!> \brief ...
!> \param a_scaling_ext ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param dimen_ia ...
!> \param Eigenval ...
!> \param num_integ_points ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param tj_ext ...
!> \param wj_ext ...
!> \param fm_mat_S ...
!> \param error ...
!> \param homo_beta ...
!> \param virtual_beta ...
!> \param dimen_ia_beta ...
!> \param Eigenval_beta ...
!> \param fm_mat_S_beta ...
! *****************************************************************************
    SUBROUTINE  calc_scaling_factor(a_scaling_ext,para_env,para_env_RPA,homo,virtual,dimen_RI,dimen_ia,Eigenval,&
                                    num_integ_points,num_integ_group,color_rpa_group,&
                                    tj_ext,wj_ext,fm_mat_S,error,&
                                    homo_beta,virtual_beta,dimen_ia_beta,Eigenval_beta,fm_mat_S_beta)
    REAL(KIND=dp)                            :: a_scaling_ext
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_RPA
    INTEGER                                  :: homo, virtual, dimen_RI, &
                                                dimen_ia
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    INTEGER                                  :: num_integ_points, &
                                                num_integ_group, &
                                                color_rpa_group
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: tj_ext, wj_ext
    TYPE(cp_fm_type), POINTER                :: fm_mat_S
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, OPTIONAL                        :: homo_beta, virtual_beta, &
                                                dimen_ia_beta
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: Eigenval_beta
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: fm_mat_S_beta

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
       PRESENT(fm_mat_S_beta)) my_open_shell=.TRUE.

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
    ! 1) get info fm_mat_S 
    !XXX CALL cp_fm_to_fm(source=fm_mat_S,destination=fm_mat_G,error=error)
    CALL cp_fm_get_info(matrix=fm_mat_S,&
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
       iaia_RI_dp(iiB)=iaia_RI_dp(iiB)+fm_mat_S%local_data(iiB,jjB)*fm_mat_S%local_data(iiB,jjB)
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
   
    color_row=fm_mat_S%matrix_struct%context%mepos(1)
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
    color_col=fm_mat_S%matrix_struct%context%mepos(2)
    CALL mp_comm_split_direct(para_env_RPA%group,comm_col,color_col)
    NULLIFY(para_env_col)
    CALL cp_para_env_create(para_env_col,comm_col,error=error)

    ! allocate communication array for columns
    CALL mp_sum(right_term_ref,para_env_col%group)    

    ! In the open shell case do point 1-2-3 for the beta spin
    IF(my_open_shell) THEN
      !XXX CALL cp_fm_to_fm(source=fm_mat_S_beta,destination=fm_mat_G_beta,error=error)
      CALL cp_fm_get_info(matrix=fm_mat_S_beta,&
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
         iaia_RI_dp_beta(iiB)=iaia_RI_dp_beta(iiB)+fm_mat_S_beta%local_data(iiB,jjB)*fm_mat_S_beta%local_data(iiB,jjB)
        END DO
      END DO

      color_row_beta=fm_mat_S_beta%matrix_struct%context%mepos(1)
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
      color_col_beta=fm_mat_S_beta%matrix_struct%context%mepos(2)
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

! *****************************************************************************
!> \brief ...
!> \param a_scaling ...
!> \param left_term ...
!> \param first_deriv ...
!> \param num_integ_points ...
!> \param my_open_shell ...
!> \param ZERO ...
!> \param ONE ...
!> \param FOUR ...
!> \param PIG ...
!> \param M_ia ...
!> \param cottj ...
!> \param wj ...
!> \param D_ia ...
!> \param D_ia_beta ...
!> \param M_ia_beta ...
!> \param nrow_local ...
!> \param nrow_local_beta ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param para_env ...
!> \param para_env_row ...
!> \param para_env_row_beta ...
! *****************************************************************************
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
