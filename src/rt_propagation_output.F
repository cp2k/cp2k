!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routine for the real time propagation output.
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

MODULE rt_propagation_output

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_create, cp_dbcsr_desymmetrize, &
       cp_dbcsr_get_occupation, cp_dbcsr_init, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_iterator,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale_and_add
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_double,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_constants,                 ONLY: ehrenfest,&
                                             real_time_propagation
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_max
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_types,                     ONLY: write_rt_mos_to_restart
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_post_dftb,                ONLY: scf_post_calculation_dftb
  USE qs_scf_post_gpw,                 ONLY: qs_scf_post_moments,&
                                             write_available_results
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_output'

  PUBLIC :: rt_prop_output,&
            rt_convergence,&
            rt_convergence_density

CONTAINS


  SUBROUTINE rt_prop_output(qs_env,run_type,delta_iter,used_time,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER                                  :: run_type
    REAL(dp), OPTIONAL                       :: delta_iter, used_time
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: im, ispin, n_electrons, &
                                                output_unit, re
    REAL(dp)                                 :: orthonormality, sparsity, &
                                                tot_rho_r
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, rho_new
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(section_vals_type), POINTER         :: dft_section, input, &
                                                rtp_section

    NULLIFY(logger)

    logger => cp_error_get_logger(error)
    CALL get_qs_env(qs_env=qs_env,rtp=rtp,matrix_s=matrix_s,input=input,rho=rho,particle_set=particle_set,&
         atomic_kind_set=atomic_kind_set,error=error)
    rtp_section => section_vals_get_subs_vals(input,"DFT%REAL_TIME_PROPAGATION",error=error)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
         nelectron=n_electrons)
    n_electrons = n_electrons - qs_env%dft_control%charge

    tot_rho_r = accurate_sum(rho%tot_rho_r)

    output_unit=cp_print_key_unit_nr(logger,rtp_section,"PRINT%PROGRAM_RUN_INFO",&
            extension=".scfLog",error=error)




    IF(output_unit>0) THEN
       WRITE (output_unit,FMT="(/,(T3,A,T40,I5))")&
            "Information at iteration step:",rtp%iter
       WRITE (UNIT=output_unit,FMT="((T3,A,T41,2F20.10))")&
            "Total electronic density (r-space): ",&
            tot_rho_r,&
            tot_rho_r + &
            REAL(n_electrons,dp)
       WRITE (UNIT=output_unit,FMT="((T3,A,T61,F20.10))")&
               "Total energy:",rtp%energy_new
       IF(run_type==ehrenfest)&
            WRITE (UNIT=output_unit,FMT="((T3,A,T61,F20.10))")&
            "Energy difference to previous iteration step:",rtp%energy_new-rtp%energy_old
       IF(run_type==real_time_propagation)&
            WRITE (UNIT=output_unit,FMT="((T3,A,T61,F20.10))")&
            "Energy difference to initial state:",rtp%energy_new-rtp%energy_old
       IF(PRESENT(delta_iter))&
       WRITE (UNIT=output_unit,FMT="((T3,A,T61,E20.6))")&
            "Convergence:",delta_iter
       IF(rtp%converged)THEN
           IF(run_type==real_time_propagation)&
               WRITE (UNIT=output_unit,FMT="((T3,A,T61,F12.2))")&
               "Time needed for propagation:",used_time
           WRITE (UNIT=output_unit,FMT="(/,(T3,A,3X,F16.10))")&
           "CONVERGENCE REACHED",rtp%energy_new-rtp%energy_old
       END IF
    END IF

    IF(rtp%converged) THEN
       IF(rtp%linear_scaling) THEN
          IF(output_unit>0) THEN
             CALL get_rtp(rtp=rtp,rho_new=rho_new,error=error)
             !TODO: Output the idempotency
             DO ispin=1,SIZE(rho_new)/2
                re=2*ispin-1
                im=2*ispin
                sparsity=(cp_dbcsr_get_occupation(rho_new(re)%matrix)+&
                         cp_dbcsr_get_occupation(rho_new(im)%matrix))/2
                WRITE(output_unit,FMT="((T3,A,I1,A,T61,F20.10))") "Sparsity of P for spin ", ispin, ":",sparsity
             ENDDO
          ENDIF
       ELSE
          CALL get_rtp(rtp=rtp,mos_new=mos_new,error=error)
          CALL rt_calculate_orthonormality(output_unit,orthonormality,&
               mos_new,matrix_s(1)%matrix,error)
          IF(output_unit>0)&
               WRITE(output_unit,FMT="(/,(T3,A,T60,F20.10))")&
               "Max deviation from orthonormalization:",orthonormality
       ENDIF
    END IF

    IF(output_unit>0)&
       CALL m_flush(output_unit)
    CALL cp_print_key_finished_output(output_unit,logger,rtp_section,&
         "PRINT%PROGRAM_RUN_INFO", error=error)

    IF(rtp%converged)THEN
       CALL make_moment(qs_env,error)
       CALL write_available_results(qs_env=qs_env,error=error)
       dft_section =>  section_vals_get_subs_vals(input,"DFT",error=error)
       !TODO: make a restart for density propagation
       IF (.NOT. rtp%linear_scaling) THEN
          CALL get_rtp(rtp=rtp,mos_new=mos_new,error=error)
          CALL write_rt_mos_to_restart(qs_env%mos,mos_new,particle_set,dft_section,atomic_kind_set,error)
       ENDIF
    END IF

  END SUBROUTINE rt_prop_output



! *****************************************************************************
!> \brief computes the effective orthonormality of a set of mos given an s-matrix
!>        orthonormality is the max deviation from unity of the C^T S C
!> \author Florian Schiffmann (02.09)
! *****************************************************************************
  SUBROUTINE rt_calculate_orthonormality(unit_nr,orthonormality,mos_new,matrix_s,error)
    INTEGER                                  :: unit_nr
    REAL(KIND=dp)                            :: orthonormality
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rt_calculate_orthonormality', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, im, ispin, j, k, &
                                                n, ncol_local, nrow_local, &
                                                nspin, re
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: alpha, max_alpha, max_beta
    TYPE(cp_fm_struct_type), POINTER         :: tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: overlap_re, svec_im, svec_re

   NULLIFY(tmp_fm_struct, svec_im, svec_re, overlap_re)

   CALL timeset(routineN,handle)

   nspin=SIZE(mos_new)/2
   max_alpha=0.0_dp
   max_beta=0.0_dp
   DO ispin=1,nspin
      re=ispin*2-1
      im=ispin*2
      ! get S*C
      CALL cp_fm_create(svec_re,mos_new(im)%matrix%matrix_struct,error=error)
      CALL cp_fm_create(svec_im,mos_new(im)%matrix%matrix_struct,error=error)
      CALL cp_fm_get_info(mos_new(im)%matrix,&
           nrow_global=n,ncol_global=k,error=error)
      CALL cp_dbcsr_sm_fm_multiply(matrix_s,mos_new(re)%matrix,&
           svec_re,k,error=error)
      CALL cp_dbcsr_sm_fm_multiply(matrix_s,mos_new(im)%matrix,&
           svec_im,k,error=error)

      ! get C^T (S*C)
      CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=k,ncol_global=k, &
           para_env=mos_new(re)%matrix%matrix_struct%para_env, &
           context=mos_new(re)%matrix%matrix_struct%context,error=error)
      CALL cp_fm_create(overlap_re,tmp_fm_struct,error=error)

      CALL cp_fm_struct_release(tmp_fm_struct,error=error)

      CALL cp_fm_gemm('T','N',k,k,n,1.0_dp, mos_new(re)%matrix,&
           svec_re,0.0_dp,overlap_re,error=error)
      CALL cp_fm_gemm('T','N',k,k,n,1.0_dp, mos_new(im)%matrix,&
           svec_im,1.0_dp,overlap_re,error=error)

      CALL cp_fm_release(svec_re,error=error)
      CALL cp_fm_release(svec_im,error=error)

      CALL cp_fm_get_info(overlap_re,nrow_local=nrow_local,ncol_local=ncol_local, &
           row_indices=row_indices,col_indices=col_indices,error=error)
      DO i=1,nrow_local
         DO j=1,ncol_local
            alpha=overlap_re%local_data(i,j)
            IF (row_indices(i).EQ.col_indices(j)) alpha=alpha-1.0_dp
            max_alpha=MAX(max_alpha,ABS(alpha))
         ENDDO
      ENDDO
      CALL cp_fm_release(overlap_re,error=error)
   ENDDO
   CALL mp_max(max_alpha,mos_new(1)%matrix%matrix_struct%para_env%group)
   CALL mp_max(max_beta,mos_new(1)%matrix%matrix_struct%para_env%group)
   orthonormality=max_alpha
   IF(unit_nr>0)THEN

   END IF


   CALL timestop(handle)

 END SUBROUTINE rt_calculate_orthonormality

! *****************************************************************************
!> \brief computs the convergence criterium for RTP and EMD
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

 SUBROUTINE rt_convergence(qs_env,delta_mos,delta_eps,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: delta_mos
    REAL(dp)                                 :: delta_eps
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rt_convergence', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, icol, im, ispin, &
                                                j, lcol, lrow, nao, newdim, &
                                                nmo, nspin, re
    LOGICAL                                  :: double_col, double_row
    REAL(KIND=dp)                            :: alpha, max_alpha
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrixv_s
    TYPE(cp_dbcsr_type), POINTER             :: matrix_s
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new
    TYPE(cp_fm_struct_type), POINTER         :: newstruct, newstruct1, &
                                                tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: work, work1, work2
    TYPE(rt_prop_type), POINTER              :: rtp

   NULLIFY(tmp_fm_struct)

   CALL timeset(routineN,handle)

   CALL get_qs_env(qs_env=qs_env,matrix_s=matrixv_s,rtp=rtp,error=error)
   CALL get_rtp(rtp=rtp,mos_new=mos_new,error=error)
   matrix_s=>matrixv_s(1)%matrix

   nspin=SIZE(delta_mos)/2
   max_alpha=0.0_dp

   DO i=1,SIZE(mos_new)
      CALL cp_fm_scale_and_add(-1.0_dp,delta_mos(i)%matrix,1.0_dp,mos_new(i)%matrix,error)
   END DO

   DO ispin=1,nspin
      re=ispin*2-1
      im=ispin*2

      double_col=.TRUE.
      double_row=.FALSE.
      CALL cp_fm_struct_double(newstruct,&
           delta_mos(re)%matrix%matrix_struct,&
           delta_mos(re)%matrix%matrix_struct%context,&
           double_col,&
           double_row,&
           error)

      CALL cp_fm_create(work,matrix_struct=newstruct,error=error)
      CALL cp_fm_create(work1,matrix_struct=newstruct,error=error)

      CALL cp_fm_get_info(delta_mos(re)%matrix,ncol_local=lcol,ncol_global=nmo,&
           nrow_global=nao,error=error)
      CALL cp_fm_get_info(work,ncol_global=newdim,error=error)

      CALL cp_fm_set_all(work,0.0_dp,0.0_dp,error)

      DO icol=1,lcol
         work%local_data(:,icol)=delta_mos(re)%matrix%local_data(:,icol)
         work%local_data(:,icol+lcol)=delta_mos(im)%matrix%local_data(:,icol)
      END DO


      CALL cp_dbcsr_sm_fm_multiply(matrix_s, work, work1, ncol=newdim, error=error)

      CALL cp_fm_release(work,error)

      CALL cp_fm_struct_create(tmp_fm_struct,nrow_global=nmo,ncol_global=nmo, &
           para_env=delta_mos(re)%matrix%matrix_struct%para_env, &
           context=delta_mos(re)%matrix%matrix_struct%context,error=error)
      CALL cp_fm_struct_double(newstruct1,&
           tmp_fm_struct,&
           delta_mos(re)%matrix%matrix_struct%context,&
           double_col,&
           double_row,&
           error)

      CALL cp_fm_create(work,matrix_struct=newstruct1,error=error)
      CALL cp_fm_create(work2,matrix_struct=newstruct1,error=error)


      CALL cp_fm_gemm("T","N",nmo,newdim,nao,1.0_dp,delta_mos(re)%matrix,&
           work1,0.0_dp,work,error=error)

      CALL cp_fm_gemm("T","N",nmo,newdim,nao,1.0_dp,delta_mos(im)%matrix,&
           work1,0.0_dp,work2,error=error)

      CALL cp_fm_get_info(work,nrow_local=lrow,error=error)
      DO i=1,lrow
         DO j=1,lcol
            alpha=SQRT((work%local_data(i,j)+work2%local_data(i,j+lcol))**2+&
                 (work%local_data(i,j+lcol)-work2%local_data(i,j))**2)
            max_alpha=MAX(max_alpha,ABS(alpha))
         ENDDO
      ENDDO

      CALL cp_fm_release(work,error)
      CALL cp_fm_release(work1,error)
      CALL cp_fm_release(work2,error)
      CALL cp_fm_struct_release(tmp_fm_struct,error=error)
      CALL cp_fm_struct_release(newstruct,error=error)
      CALL cp_fm_struct_release(newstruct1,error=error)

   ENDDO

   CALL mp_max(max_alpha,delta_mos(1)%matrix%matrix_struct%para_env%group)
   delta_eps=SQRT(max_alpha)

   CALL timestop(handle)

 END SUBROUTINE rt_convergence

! *****************************************************************************
!> \brief computs the convergence criterium for RTP and EMD based on the density matrix
!> \author Samuel Andermatt (02.14)
! *****************************************************************************

 SUBROUTINE rt_convergence_density(qs_env,delta_P,delta_eps,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: delta_P
    REAL(dp)                                 :: delta_eps
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rt_convergence_density', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: col_atom, handle, i, ispin, &
                                                row_atom
    LOGICAL                                  :: failure
    REAL(dp)                                 :: alpha, max_alpha
    REAL(dp), DIMENSION(:), POINTER          :: block_values
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_new
    TYPE(cp_dbcsr_type), POINTER             :: tmp
    TYPE(rt_prop_type), POINTER              :: rtp

   CALL timeset(routineN,handle)

   failure=.FALSE.

   CALL get_qs_env(qs_env=qs_env,rtp=rtp,error=error)
   CALL get_rtp(rtp=rtp,rho_new=rho_new,error=error)

   DO i=1,SIZE(rho_new)
      CALL cp_dbcsr_add(delta_P(i)%matrix,rho_new(i)%matrix,one,-one,error=error)
   ENDDO
   !get the maximum value of delta_P
   DO i=1,SIZE(delta_P)
      !square all entries of both matrices
       CALL cp_dbcsr_iterator_start(iter, delta_P(i)%matrix)
       DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
          CALL cp_dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
          block_values=block_values*block_values
       END DO
       CALL cp_dbcsr_iterator_stop (iter)
   END DO
   NULLIFY(tmp)
   ALLOCATE(tmp)
   CALL cp_dbcsr_init(tmp,error=error)
   CALL cp_dbcsr_create(tmp,template=delta_P(1)%matrix,matrix_type="N",error=error)
   DO ispin=1,SIZE(delta_P)/2
      CALL cp_dbcsr_desymmetrize(delta_P(2*ispin-1)%matrix,tmp,error=error)
      CALL cp_dbcsr_add(delta_P(2*ispin)%matrix,tmp,one,one,error=error)
   END DO
   !the absolute values are now in the even entries of delta_P
   max_alpha=zero
   DO ispin=1,SIZE(delta_P)/2
       CALL cp_dbcsr_iterator_start(iter, delta_P(2*ispin)%matrix)
       DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
          CALL cp_dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
          alpha=MAXVAL(block_values)
          IF(alpha>max_alpha) max_alpha=alpha
       END DO
       CALL cp_dbcsr_iterator_stop (iter)
   END DO
   CALL mp_max(max_alpha,delta_P(1)%matrix%matrix%m%dist%d%mp_env%mp%mp_group)
   delta_eps=SQRT(max_alpha)
   CALL cp_dbcsr_deallocate_matrix(tmp,error=error)
   CALL timestop(handle)

 END SUBROUTINE rt_convergence_density

! *****************************************************************************
!> \brief interface to qs_moments. Does only work for nonperiodic dipole
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE make_moment(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_moment', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, output_unit
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: dft_section, input

    CALL timeset(routineN,handle)

    NULLIFY(input)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    CALL get_qs_env(qs_env,input=input,error=error)
    dft_section => section_vals_get_subs_vals(qs_env%input,"DFT",error=error)
    IF(qs_env%dft_control%qs_control%dftb)THEN
       CALL scf_post_calculation_dftb(dft_section, qs_env=qs_env, error=error) 
    ELSE
       CALL qs_scf_post_moments(input, logger, qs_env, output_unit, error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE make_moment






END MODULE rt_propagation_output
