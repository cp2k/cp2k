!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Rountines to calculate MP2 energy with laplace approach
!> \par History
!>      11.2012 created [Mauro Del Ben]
! *****************************************************************************
MODULE mp2_laplace
! 
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_set_all,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE message_passing,                 ONLY: mp_sum
  USE minimax,                         ONLY: get_minimax_coeff
  USE mp2_types,                       ONLY: integ_mat_buffer_type
  USE rpa_communication,               ONLY: fm_redistribute,&
                                             initialize_buffer,&
                                             release_buffer
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_laplace'

  PUBLIC :: laplace_minimax_approx


  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param Emp2 ...
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
!> \param error ...
!> \param homo_beta ...
!> \param virtual_beta ...
!> \param dimen_ia_beta ...
!> \param Eigenval_beta ...
!> \param fm_mat_S_beta ...
!> \param fm_mat_Q_gemm_beta ...
!> \param fm_mat_Q_beta ...
! *****************************************************************************
  SUBROUTINE laplace_minimax_approx(Emp2,para_env,para_env_RPA,unit_nr,homo,virtual,dimen_RI,dimen_ia,Eigenval,&
                                    num_integ_points,num_integ_group,color_rpa_group,&
                                    fm_mat_S,fm_mat_Q_gemm,fm_mat_Q,error,&
                                    homo_beta,virtual_beta,dimen_ia_beta,Eigenval_beta,fm_mat_S_beta,&
                                    fm_mat_Q_gemm_beta,fm_mat_Q_beta)
    REAL(KIND=dp)                            :: Emp2
    TYPE(cp_para_env_type), POINTER          :: para_env, para_env_RPA
    INTEGER                                  :: unit_nr, homo, virtual, &
                                                dimen_RI, dimen_ia
    REAL(KIND=dp), DIMENSION(:)              :: Eigenval
    INTEGER                                  :: num_integ_points, &
                                                num_integ_group, &
                                                color_rpa_group
    TYPE(cp_fm_type), POINTER                :: fm_mat_S, fm_mat_Q_gemm, &
                                                fm_mat_Q
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, OPTIONAL                        :: homo_beta, virtual_beta, &
                                                dimen_ia_beta
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: Eigenval_beta
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: fm_mat_S_beta, &
                                                fm_mat_Q_gemm_beta, &
                                                fm_mat_Q_beta

    CHARACTER(LEN=*), PARAMETER :: routineN = 'laplace_minimax_approx', &
      routineP = moduleN//':'//routineN

    INTEGER :: avirt, handle, i_global, ierr, iiB, iocc, j_global, jjB, &
      jquad, my_num_dgemm_call, ncol_local, nrow_local, number_of_rec, &
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
    REAL(KIND=dp) :: actual_flop_rate, alpha_old, E_Range, Emax, Emax_beta, &
      Emin, Emin_beta, laplace_transf, my_flop_rate, t_end, t_start, &
      trace_XX, w_old
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: aj, awj, wj
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

    ! Calculate weights and exponents with minimax approximation
    ALLOCATE(awj(2*num_integ_points),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    awj=0.0_dp

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
    CALL get_minimax_coeff(num_integ_points,E_Range,awj,ierr,error)

    ALLOCATE(aj(num_integ_points),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    aj=0.0_dp
    ALLOCATE(wj(num_integ_points),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    wj=0.0_dp

    DO jquad=1, num_integ_points
      aj(jquad)=awj(jquad)
      wj(jquad)=awj(jquad+num_integ_points)
    END DO

    DEALLOCATE(awj)
  
    IF (unit_nr>0) THEN
      WRITE (UNIT=unit_nr,FMT="(T3,A,T66,F15.4)")&
                          "MINIMAX| Range of the energy denominator:", E_Range
      WRITE (UNIT=unit_nr,FMT="(T3,A,T54,A,T72,A)") "MINIMAX| Parameters:","Weights","Exponents"
      DO jquad=1, num_integ_points
        WRITE (UNIT=unit_nr,FMT="(T41,F20.10,F20.10)") wj(jquad), aj(jquad)
      END DO
      CALL m_flush(unit_nr)
    END IF

    ! scale the minimax parameters
    aj=aj/Emin
    wj=wj/Emin

    ! initialize buffer for matrix redistribution
    CALL initialize_buffer(fm_mat_Q_gemm,fm_mat_Q,RPA_proc_map,buffer_rec,buffer_send,&
                           number_of_rec,number_of_send,&
                           map_send_size,map_rec_size,local_size_source,para_env_RPA,error)
    IF(my_open_shell) THEN
      CALL initialize_buffer(fm_mat_Q_gemm_beta,fm_mat_Q_beta,RPA_proc_map,buffer_rec_beta,buffer_send_beta,&
                             number_of_rec_beta,number_of_send_beta,&
                             map_send_size_beta,map_rec_size_beta,local_size_source_beta,para_env_RPA,error)
    END IF    

    Emp2=0.0_dp
    first_cycle=.TRUE.
    alpha_old=0.0_dp
    w_old=0.0_dp
    my_num_dgemm_call=0
    my_flop_rate=0.0_dp
    DO jquad=1, num_integ_points
      IF(MODULO(jquad,num_integ_group)/=color_rpa_group) CYCLE
      
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

      ! update G matrix with the new value of w and a
      IF(first_cycle) THEN
        DO jjB=1, ncol_local
          j_global=col_indices(jjB)
          DO iiB=1, nrow_local
            i_global=row_indices(iiB)

            iocc=MAX(1,i_global-1)/virtual+1
            avirt=i_global-(iocc-1)*virtual

            laplace_transf=SQRT(wj(jquad))*EXP(Eigenval(iocc)*aj(jquad))*EXP(-Eigenval(avirt+homo)*aj(jquad))
            laplace_transf=SQRT(laplace_transf)

            fm_mat_S%local_data(iiB,jjB)=fm_mat_S%local_data(iiB,jjB)*laplace_transf

          END DO
        END DO
      ELSE
        DO jjB=1, ncol_local
          j_global=col_indices(jjB)
          DO iiB=1, nrow_local
            i_global=row_indices(iiB)

            iocc=MAX(1,i_global-1)/virtual+1
            avirt=i_global-(iocc-1)*virtual

            laplace_transf=SQRT(wj(jquad)/w_old)*&
                           EXP(Eigenval(iocc)*(aj(jquad)-alpha_old))*EXP(-Eigenval(avirt+homo)*(aj(jquad)-alpha_old))
            laplace_transf=SQRT(laplace_transf)

            fm_mat_S%local_data(iiB,jjB)=fm_mat_S%local_data(iiB,jjB)*laplace_transf

          END DO
        END DO
      END IF

      t_start=m_walltime()
      CALL cp_gemm(transa="T",transb="N",m=dimen_RI,n=dimen_RI,k=dimen_ia,alpha=1.0_dp,&
                      matrix_a=fm_mat_S,matrix_b=fm_mat_S,beta=0.0_dp,&
                      matrix_c=fm_mat_Q_gemm,error=error)
      t_end=m_walltime()
      actual_flop_rate=2.0_dp*REAL(dimen_ia,KIND=dp)*dimen_RI*REAL(dimen_RI,KIND=dp)/(MAX(0.01_dp,t_end-t_start))
      IF(para_env_RPA%mepos==0) my_flop_rate=my_flop_rate+actual_flop_rate
      my_num_dgemm_call=my_num_dgemm_call+1

      ! copy/redistribute fm_mat_Q_gemm to fm_mat_Q
      CALL cp_fm_set_all(matrix=fm_mat_Q,alpha=0.0_dp,error=error)
      CALL fm_redistribute(fm_mat_Q_gemm,fm_mat_Q,RPA_proc_map,buffer_rec,buffer_send,&
                           number_of_rec,number_of_send,&
                           map_send_size,map_rec_size,local_size_source,para_env_RPA,error)

      IF(my_open_shell) THEN
        !XXX  the same for the beta spin
        !XXX ! copy fm_mat_S into fm_mat_G
        !XXX CALL cp_fm_to_fm(source=fm_mat_S_beta,destination=fm_mat_G_beta,error=error)
        !XXX ! get info of fm_mat_G_beta
        !XXX CALL cp_fm_get_info(matrix=fm_mat_G_beta,&
        !XXX                     nrow_local=nrow_local,&
        !XXX                     ncol_local=ncol_local,&
        !XXX                     row_indices=row_indices,&
        !XXX                     col_indices=col_indices,&
        !XXX                     error=error)

        ! the same for the beta spin
        CALL cp_fm_get_info(matrix=fm_mat_S_beta,&
                            nrow_local=nrow_local,&
                            ncol_local=ncol_local,&
                            row_indices=row_indices,&
                            col_indices=col_indices,&
                            error=error)
        ! update G matrix with the new value of w and a
        IF(first_cycle) THEN
          DO jjB=1, ncol_local
            j_global=col_indices(jjB)
            DO iiB=1, nrow_local
              i_global=row_indices(iiB)

              iocc=MAX(1,i_global-1)/virtual_beta+1
              avirt=i_global-(iocc-1)*virtual_beta

              laplace_transf=SQRT(wj(jquad))*EXP(Eigenval_beta(iocc)*aj(jquad))*&
                                             EXP(-Eigenval_beta(avirt+homo_beta)*aj(jquad))
              laplace_transf=SQRT(laplace_transf)

              fm_mat_S_beta%local_data(iiB,jjB)=fm_mat_S_beta%local_data(iiB,jjB)*laplace_transf

            END DO
          END DO
        ELSE
          DO jjB=1, ncol_local
            j_global=col_indices(jjB)
            DO iiB=1, nrow_local
              i_global=row_indices(iiB)

              iocc=MAX(1,i_global-1)/virtual_beta+1
              avirt=i_global-(iocc-1)*virtual_beta

              laplace_transf=SQRT(wj(jquad)/w_old)*&
                             EXP(Eigenval_beta(iocc)*(aj(jquad)-alpha_old))*&
                             EXP(-Eigenval_beta(avirt+homo_beta)*(aj(jquad)-alpha_old))
              laplace_transf=SQRT(laplace_transf)

              fm_mat_S_beta%local_data(iiB,jjB)=fm_mat_S_beta%local_data(iiB,jjB)*laplace_transf

            END DO
          END DO
        END IF

        t_start=m_walltime()
        CALL cp_gemm(transa="T",transb="N",m=dimen_RI,n=dimen_RI,k=dimen_ia_beta,alpha=1.0_dp,&
                        matrix_a=fm_mat_S_beta,matrix_b=fm_mat_S_beta,beta=0.0_dp,&
                        matrix_c=fm_mat_Q_gemm_beta,error=error)
        t_end=m_walltime()
        actual_flop_rate=2.0_dp*REAL(dimen_ia_beta,KIND=dp)*dimen_RI*REAL(dimen_RI,KIND=dp)/(MAX(0.01_dp,t_end-t_start))
        IF(para_env_RPA%mepos==0) my_flop_rate=my_flop_rate+actual_flop_rate
        my_num_dgemm_call=my_num_dgemm_call+1

        ! copy/redistribute fm_mat_Q_gemm to fm_mat_Q
        CALL cp_fm_set_all(matrix=fm_mat_Q_beta,alpha=0.0_dp,error=error)
        CALL fm_redistribute(fm_mat_Q_gemm_beta,fm_mat_Q_beta,RPA_proc_map,buffer_rec_beta,buffer_send_beta,&
                             number_of_rec_beta,number_of_send_beta,&
                             map_send_size_beta,map_rec_size_beta,local_size_source_beta,para_env_RPA,error)

      END IF

      ! get info of fm_mat_Q
      CALL cp_fm_get_info(matrix=fm_mat_Q,&
                          nrow_local=nrow_local,&
                          ncol_local=ncol_local,&
                          row_indices=row_indices,&
                          col_indices=col_indices,&
                          error=error)

      ! calcualte the trace of the product Q*Q
      trace_XX=0.0_dp
      DO jjB=1, ncol_local
        j_global=col_indices(jjB)
        DO iiB=1, nrow_local
          i_global=row_indices(iiB)
          IF(my_open_shell) THEN
            trace_XX=trace_XX+fm_mat_Q%local_data(iiB,jjB)*fm_mat_Q_beta%local_data(iiB,jjB)
          ELSE
            trace_XX=trace_XX+fm_mat_Q%local_data(iiB,jjB)*fm_mat_Q%local_data(iiB,jjB)
          END IF
        END DO
      END DO

      Emp2=Emp2-trace_XX

      ! save alpha, w and reset the first_cycle flag
      first_cycle=.FALSE.
      alpha_old=aj(jquad)
      w_old=wj(jquad)
    END DO

    CALL mp_sum(Emp2,para_env%group)

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

    DEALLOCATE(aj)
    DEALLOCATE(wj)

    CALL timestop(handle)

    END SUBROUTINE laplace_minimax_approx

END MODULE mp2_laplace
