!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routine for the real time propagation output.
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

MODULE rt_propagation_output
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_allocate_matrix_set, cp_dbcsr_binary_write, &
       cp_dbcsr_checksum, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_deallocate_matrix, cp_dbcsr_deallocate_matrix_set, &
       cp_dbcsr_desymmetrize, cp_dbcsr_distribution, cp_dbcsr_filter, &
       cp_dbcsr_get_occupation, cp_dbcsr_init, cp_dbcsr_init_p, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_p_type, cp_dbcsr_set, cp_dbcsr_type, &
       dbcsr_distribution_mp, dbcsr_mp_group
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
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
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
  USE input_constants,                 ONLY: ehrenfest,&
                                             real_time_propagation
  USE input_section_types,             ONLY: section_get_ivals,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_max
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_copy
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_release
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                             qs_kind_type
  USE qs_linres_current,               ONLY: calculate_jrho_resp
  USE qs_linres_types,                 ONLY: current_env_type
  USE qs_mo_io,                        ONLY: write_rt_mos_to_restart
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_scf_post_dftb,                ONLY: scf_post_calculation_dftb
  USE qs_scf_post_gpw,                 ONLY: qs_scf_post_moments,&
                                             write_available_results,&
                                             write_mo_free_results
  USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                             qs_subsys_type
  USE realspace_grid_types,            ONLY: realspace_grid_desc_type
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_output'

  PUBLIC :: rt_prop_output,&
            rt_convergence,&
            rt_convergence_density,&
            report_density_occupation

CONTAINS


! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param run_type ...
!> \param delta_iter ...
!> \param used_time ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rt_prop_output(qs_env,run_type,delta_iter,used_time,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(in)                      :: run_type
    REAL(dp), INTENT(in), OPTIONAL           :: delta_iter, used_time
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rt_prop_output', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n_electrons, output_unit
    REAL(dp)                                 :: orthonormality, tot_rho_r
    REAL(KIND=dp), DIMENSION(:), POINTER     :: qs_tot_rho_r
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, rho_new
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(section_vals_type), POINTER         :: dft_section, input, &
                                                rtp_section

    NULLIFY(logger, dft_control)
    
    logger => cp_error_get_logger(error)
    CALL get_qs_env(qs_env,&
                    rtp=rtp,&
                    matrix_s=matrix_s,&
                    input=input,&
                    rho=rho,&
                    particle_set=particle_set,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    dft_control=dft_control,&
                    error=error)

    rtp_section => section_vals_get_subs_vals(input,"DFT%REAL_TIME_PROPAGATION",error=error)

    CALL get_qs_kind_set(qs_kind_set, nelectron=n_electrons)
    n_electrons = n_electrons - dft_control%charge

    CALL qs_rho_get(rho_struct=rho,tot_rho_r=qs_tot_rho_r,error=error)

    tot_rho_r = accurate_sum(qs_tot_rho_r)

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
       IF(.NOT.rtp%linear_scaling) THEN
          CALL get_rtp(rtp=rtp,mos_new=mos_new,error=error)
          CALL rt_calculate_orthonormality(orthonormality,&
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
       dft_section =>  section_vals_get_subs_vals(input,"DFT",error=error)
       IF (rtp%linear_scaling) THEN
          CALL get_rtp(rtp=rtp,rho_new=rho_new,error=error)
          IF (rtp%write_restart) CALL write_rt_p_to_restart(rho_new,error) !At some point a print section should control this output
          IF(.NOT.dft_control%qs_control%dftb) THEN
             !Not sure if these things could also work with dftb or not
             CALL write_mo_free_results(qs_env,error)
             IF(BTEST(cp_print_key_should_output(logger%iter_info,&
                dft_section,"REAL_TIME_PROPAGATION%PRINT%CURRENT",error=error),cp_p_file)) THEN
                   CALL rt_current(qs_env,rho_new(2)%matrix,dft_section,error) !not yet implemented for spin polarized systems and full matrix rtp
             ENDIF
          ENDIF
       ELSE
          CALL get_rtp(rtp=rtp,mos_new=mos_new,error=error)
          IF(.NOT.dft_control%qs_control%dftb) THEN
             CALL write_available_results(qs_env=qs_env,error=error)
          ENDIF
          CALL write_rt_mos_to_restart(qs_env%mos,mos_new,particle_set,&
                            dft_section,atomic_kind_set,qs_kind_set,error)
       ENDIF
    ENDIF
    
    rtp%energy_old=rtp%energy_new
    
    CALL cp_assert(rtp%converged.OR.rtp%iter<dft_control%rtp_control%max_iter,cp_failure_level,cp_assertion_failed,&
         routineP,"EMD did not converge, either increase MAX_ITER "//&
         "or use a smaller TIMESTEP",&
         error)

  END SUBROUTINE rt_prop_output



! *****************************************************************************
!> \brief computes the effective orthonormality of a set of mos given an s-matrix
!>        orthonormality is the max deviation from unity of the C^T S C
!> \param orthonormality ...
!> \param mos_new ...
!> \param matrix_s ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************
  SUBROUTINE rt_calculate_orthonormality(orthonormality,mos_new,matrix_s,error)
    REAL(KIND=dp), INTENT(out)               :: orthonormality
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

      CALL cp_gemm('T','N',k,k,n,1.0_dp, mos_new(re)%matrix,&
           svec_re,0.0_dp,overlap_re,error=error)
      CALL cp_gemm('T','N',k,k,n,1.0_dp, mos_new(im)%matrix,&
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

   CALL timestop(handle)

 END SUBROUTINE rt_calculate_orthonormality

! *****************************************************************************
!> \brief computs the convergence criterium for RTP and EMD
!> \param rtp ...
!> \param matrix_s Overlap matrix without the derivatives
!> \param delta_mos ...
!> \param delta_eps ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

 SUBROUTINE rt_convergence(rtp,matrix_s,delta_mos,delta_eps,error)

    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(cp_dbcsr_type), POINTER             :: matrix_s
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: delta_mos
    REAL(dp), INTENT(out)                    :: delta_eps
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rt_convergence', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, i, icol, im, ispin, &
                                                j, lcol, lrow, nao, newdim, &
                                                nmo, nspin, re
    LOGICAL                                  :: double_col, double_row
    REAL(KIND=dp)                            :: alpha, max_alpha
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new
    TYPE(cp_fm_struct_type), POINTER         :: newstruct, newstruct1, &
                                                tmp_fm_struct
    TYPE(cp_fm_type), POINTER                :: work, work1, work2

   NULLIFY(tmp_fm_struct)

   CALL timeset(routineN,handle)

   CALL get_rtp(rtp=rtp,mos_new=mos_new,error=error)

   nspin=SIZE(delta_mos)/2
   max_alpha=0.0_dp

   DO i=1,SIZE(mos_new)
      CALL cp_fm_scale_and_add(-one,delta_mos(i)%matrix,one,mos_new(i)%matrix,error)
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

      CALL cp_fm_set_all(work,zero,zero,error)

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


      CALL cp_gemm("T","N",nmo,newdim,nao,one,delta_mos(re)%matrix,&
           work1,zero,work,error=error)

      CALL cp_gemm("T","N",nmo,newdim,nao,one,delta_mos(im)%matrix,&
           work1,zero,work2,error=error)

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
!> \param rtp ...
!> \param delta_P ...
!> \param delta_eps ...
!> \param error ...
!> \author Samuel Andermatt (02.14)
! *****************************************************************************

 SUBROUTINE rt_convergence_density(rtp,delta_P,delta_eps,error)

    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: delta_P
    REAL(dp), INTENT(out)                    :: delta_eps
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

   CALL timeset(routineN,handle)

   failure=.FALSE.

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
   CALL mp_max(max_alpha,dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(delta_P(1)%matrix))))
   delta_eps=SQRT(max_alpha)
   CALL cp_dbcsr_deallocate_matrix(tmp,error=error)
   CALL timestop(handle)

 END SUBROUTINE rt_convergence_density

! *****************************************************************************
!> \brief interface to qs_moments. Does only work for nonperiodic dipole
!> \param qs_env ...
!> \param error ...
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE make_moment(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_moment', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, output_unit
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(section_vals_type), POINTER         :: dft_section, input

    CALL timeset(routineN,handle)

    NULLIFY(input, dft_control)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    CALL get_qs_env(qs_env, input=input, dft_control=dft_control, error=error)
    dft_section => section_vals_get_subs_vals(qs_env%input,"DFT",error=error)
    IF(dft_control%qs_control%dftb)THEN
       CALL scf_post_calculation_dftb(dft_section, qs_env=qs_env, error=error) 
    ELSE
       CALL qs_scf_post_moments(input, logger, qs_env, output_unit, error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE make_moment

! *****************************************************************************
!> \brief Reports the sparsity pattern of the complex density matrix
!> \param filter_eps ...
!> \param rho ...
!> \param error ...
!> \author Samuel Andermatt (09.14)
! *****************************************************************************

  SUBROUTINE report_density_occupation(filter_eps,rho,error)

    REAL(KIND=dp)                            :: filter_eps
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'report_density_occupation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, im, ispin, re, &
                                                unit_nr
    REAL(KIND=dp)                            :: eps, occ
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: tmp
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    unit_nr = cp_logger_get_default_io_unit(logger)
    NULLIFY(tmp)
    CALL cp_dbcsr_allocate_matrix_set(tmp,SIZE(rho),error=error)
    DO i=1,SIZE(rho)
      CALL cp_dbcsr_init_p(tmp(i)%matrix,error=error)
      CALL cp_dbcsr_create(tmp(i)%matrix,template=rho(i)%matrix,error=error)
      CALL cp_dbcsr_copy(tmp(i)%matrix,rho(i)%matrix,error=error)
    ENDDO
    DO ispin=1,SIZE(rho)/2
      re=2*ispin-1
      im=2*ispin
      eps=MAX(filter_eps,10E-12_dp)
      DO WHILE(eps<1.1_dp)
         CALL cp_dbcsr_filter(tmp(re)%matrix,eps,error=error)
         occ=cp_dbcsr_get_occupation(tmp(re)%matrix)
         IF(unit_nr>0) WRITE(unit_nr,FMT="((T3,A,I1,A,F15.12,A,T61,F20.10))") "Occupation of rho spin ",&
                       ispin," eps ", eps, " real: ", occ
         eps=eps*10 
      ENDDO
      eps=MAX(filter_eps,10E-12_dp)
      DO WHILE(eps<1.1_dp)
         CALL cp_dbcsr_filter(tmp(im)%matrix,eps,error=error)
         occ=cp_dbcsr_get_occupation(tmp(im)%matrix)
         IF(unit_nr>0) WRITE(unit_nr,FMT="((T3,A,I1,A,F15.12,A,T61,F20.10))") "Occupation of rho spin ",&
                       ispin," eps ", eps, " imag: " ,occ
         eps=eps*10 
      ENDDO
    ENDDO
    CALL cp_dbcsr_deallocate_matrix_set(tmp,error=error)
    CALL timestop(handle)

  END SUBROUTINE report_density_occupation

! *****************************************************************************
!> \brief Writes the density matrix and the atomic positions to a restart file
!> \param rho_new ...
!> \param error ...
!> \author Samuel Andermatt (09.14)
! *****************************************************************************

  SUBROUTINE write_rt_p_to_restart(rho_new,error)

    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_new
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_rt_p_to_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name, project_name
    INTEGER                                  :: handle, im, ispin, re, unit_nr
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: cs_pos
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure = .FALSE.
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF

    project_name = logger%iter_info%project_name
    DO ispin=1,SIZE(rho_new)/2
       re=2*ispin-1
       im=2*ispin
       WRITE(file_name,'(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_RE",ispin,"_RESTART.dm"
       cs_pos = cp_dbcsr_checksum (rho_new(re)%matrix, pos=.TRUE., error=error)
       IF (unit_nr>0) THEN
          WRITE(unit_nr,'(T2,A,E20.8)') "Writing restart DM "//TRIM(file_name)//" with checksum: ",cs_pos
       ENDIF
       CALL cp_dbcsr_binary_write(rho_new(re)%matrix,file_name,error)
       WRITE(file_name,'(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_IM",ispin,"_RESTART.dm"
       cs_pos = cp_dbcsr_checksum (rho_new(im)%matrix, pos=.TRUE., error=error)
       IF (unit_nr>0) THEN
          WRITE(unit_nr,'(T2,A,E20.8)') "Writing restart DM "//TRIM(file_name)//" with checksum: ",cs_pos
       ENDIF
       CALL cp_dbcsr_binary_write(rho_new(im)%matrix,file_name,error)
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE write_rt_p_to_restart
  
! *****************************************************************************
!> \brief Collocation of the current and printing of it in a cube file
!> \param qs_env ...
!> \param P_im ...
!> \param dft_section ...
!> \param error ...
!> \author Samuel Andermatt (06.15)
! *****************************************************************************
  SUBROUTINE rt_current(qs_env,P_im,dft_section,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_type), POINTER             :: P_im
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rt_current', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=2)                         :: sdir
    CHARACTER(len=7)                         :: ext
    INTEGER                                  :: dir, handle, print_unit
    INTEGER, DIMENSION(:), POINTER           :: stride(:)
    TYPE(cp_dbcsr_type), POINTER             :: zero
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(current_env_type)                   :: current_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, rho_r
    TYPE(pw_p_type), POINTER                 :: gs, rs
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    CALL get_qs_env(qs_env=qs_env,subsys=subsys,pw_env=pw_env,rho=rho,error=error)
    CALL pw_env_get(pw_env, auxbas_rs_desc=auxbas_rs_desc,auxbas_pw_pool=auxbas_pw_pool,error=error)

    NULLIFY(zero,rs,gs)
    ALLOCATE(zero)
    CALL cp_dbcsr_init(zero,error=error)
    CALL cp_dbcsr_create(zero,template=P_im,error=error)
    CALL cp_dbcsr_copy(zero,P_im,error=error)
    CALL cp_dbcsr_set(zero,0.0_dp,error=error)
    current_env%gauge=-1
    
    NULLIFY(stride)
    ALLOCATE(stride(3))
    
    DO dir=1,3
       NULLIFY(rs,gs)
       ALLOCATE(rs)
       ALLOCATE(gs)
       CALL qs_rho_get(rho,rho_r=rho_r,rho_g=rho_g,error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool,rs%pw,use_data=REALDATA3D,in_space=REALSPACE,error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool,gs%pw,use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
       CALL pw_copy(rho_r(1)%pw,rs%pw,error)
       CALL pw_copy(rho_g(1)%pw,gs%pw,error)
    
       CALL calculate_jrho_resp(zero,P_im,zero,zero,dir,dir,rs,gs,qs_env,current_env,error=error)
    
       CALL qs_subsys_get(subsys,particles=particles,error=error) 
    
       stride=section_get_ivals(dft_section,"REAL_TIME_PROPAGATION%PRINT%CURRENT%STRIDE",error=error)
       
       IF(dir==1) THEN
          sdir="-x"
       ELSEIF(dir==2) THEN
          sdir="-y"
       ELSE
          sdir="-z"
       ENDIF
       
       ext=sdir//".cube"
       
       print_unit=cp_print_key_unit_nr(logger,dft_section,"REAL_TIME_PROPAGATION%PRINT%CURRENT",&
                                       extension=ext, file_status="REPLACE", file_action="WRITE",&
                                       log_filename=.FALSE.,error=error)
       
       CALL cp_pw_to_cube(rs%pw,print_unit,"EMD current",particles=particles,stride=stride,error=error)
    
       CALL pw_release(rs%pw,error=error)
       CALL pw_release(gs%pw,error=error)
       DEALLOCATE(rs)
       DEALLOCATE(gs)
    END DO
    
    CALL cp_dbcsr_deallocate_matrix(zero,error=error)

    DEALLOCATE(stride)
    
    CALL timestop(handle)

  END SUBROUTINE rt_current
  

END MODULE rt_propagation_output
