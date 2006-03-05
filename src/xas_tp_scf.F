!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****s* cp2k/xas_tp_scf   [1.0] *
!!
!!   NAME
!!      xas_methods
!!
!!   FUNCTION
!!      xas_scf for the tp method
!!      It is repeaated for every atom that have to be excited
!!
!!   AUTHOR
!!     MI (05.2005)
!!
!!   MODIFICATION HISTORY
!!     created 05.2005
!!
!!   SOURCE
!******************************************************************************


MODULE xas_tp_scf
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_init_random,&
                                             cp_fm_p_type,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: core_guess,&
                                             do_loc_none,&
                                             ot_precond_full_kinetic,&
                                             xas_dscf,&
                                             xas_scf_general,&
                                             xas_scf_ot,&
                                             xas_tp_xfh,&
                                             xas_tp_xhh
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE particle_types,                  ONLY: particle_type
  USE preconditioner,                  ONLY: destroy_preconditioner,&
                                             init_preconditioner,&
                                             make_preconditioner,&
                                             preconditioner_type
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_diis,                         ONLY: qs_diis_b_step
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_loc_control,                  ONLY: localized_wfn_control_type
  USE qs_loc_methods,                  ONLY: optimize_loc_berry
  USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                             qs_loc_env_new_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type,&
                                             set_mo_set,&
                                             write_mo_set
  USE qs_ot_eigensolver,               ONLY: ot_eigensolver
  USE qs_ot_scf,                       ONLY: ot_scf_mini
  USE qs_rho_atom_methods,             ONLY: calculate_rho_atom_coeff
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_methods,                  ONLY: eigensolver,&
                                             scf_env_density_mixing
  USE qs_scf_types,                    ONLY: ao_mixing_nr,&
                                             qs_scf_env_type
  USE scf_control_types,               ONLY: scf_control_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type,&
                                             real_matrix_type,&
                                             transfer_matrix
  USE termination,                     ONLY: external_control
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_control,                     ONLY: xas_control_type
  USE xas_env_types,                   ONLY: get_xas_env,&
                                             xas_environment_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tp_scf'  

! *** Public subroutines ***

  PUBLIC :: xas_do_tp_scf 
 
!****************************************************************************

CONTAINS


!!****f* xas_methods/xas_do_tp_scf [1.0] *
!!
!!   NAME
!!     xas_do_tp_scf
!!
!!   FUNCTION
!!     perform an scf loop to calculate the xas spectrum 
!!     given by the excitation of a inner state of a selected atom
!!     by using the transition potential method
!!
!!   NOTES
!!
!!   ARGUMENTS
!!     - xas_env: the environment for XAS  calculations 
!!     - scf_env: the scf_env where to perform the scf procedure
!!     - qs_env: the qs_env, the scf_env and xas_env live in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     05.2005 created [MI]
!!
!!*** **********************************************************************

  SUBROUTINE xas_do_tp_scf (xas_control,xas_env,iatom,scf_env,qs_env,globenv,error)

    TYPE(xas_control_type)                   :: xas_control
    TYPE(xas_environment_type), POINTER      :: xas_env
    INTEGER, INTENT(IN)                      :: iatom
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xas_do_tp_scf', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, handle2, homo, i, ic_max, ir_max, ispin, istat, &
      istate, iter_count, m, n, nao, nmo, nvirtual1, nvirtual2, output_unit, &
      uno_iter
    INTEGER, DIMENSION(:), POINTER           :: col_indices
    LOGICAL :: diis_step, do_level_shift, energy_only, failure, gapw, &
      id_equal, ionode, should_stop, uniform_occupation, use_cholesky, &
      use_jacobi
    REAL(dp), DIMENSION(6)                   :: weights
    REAL(dp), DIMENSION(:, :), POINTER       :: centers_wfn
    REAL(KIND=dp)                            :: a_max, b_max, diis_error, &
                                                IP_energy, ra(3), rac(3), &
                                                rc(3), t1, t2, tot1_h, &
                                                tot1_s, uno_eps
    REAL(KIND=dp), DIMENSION(:), POINTER     :: all_evals, eigenvalues, &
                                                occupation_numbers, uno_evals
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer, vecbuffer2
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: moloc_coeff, scf_work1
    TYPE(cp_fm_p_type), DIMENSION(:, :), &
      POINTER                                :: op_fm_set
    TYPE(cp_fm_type), POINTER                :: all_vectors, excvec_coeff, &
                                                excvec_overlap, fm_work, &
                                                mo_coeff, ortho, scf_work2, &
                                                uno_orbs
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(preconditioner_type), POINTER       :: preconditioner
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_ks, &
                                                matrix_p, matrix_s
    TYPE(real_matrix_p_type), &
      DIMENSION(:, :), POINTER               :: op_sm_set
    TYPE(real_matrix_type), POINTER          :: orthogonality_metric
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: dft_section, input, &
                                                scf_section

    CALL timeset(routineN,"I"," ",handle)
    NULLIFY(dft_control,matrix_h,matrix_s,matrix_ks,matrix_p)
    NULLIFY(rho,energy,scf_control,logger, ks_env,mos,atomic_kind_set)
    NULLIFY(particle_set,vecbuffer,vecbuffer2, occupation_numbers)
    NULLIFY( scf_work1,ortho,scf_work2,excvec_coeff,excvec_overlap,mo_coeff)
    NULLIFY(fm_work, col_indices, eigenvalues)
    NULLIFY(centers_wfn,qs_loc_env,localized_wfn_control)
    NULLIFY(op_sm_set,op_fm_set,moloc_coeff,cell)
    NULLIFY( dft_section, scf_section,input,qs_charges)

    logger => cp_error_get_logger(error)
    t1 = m_walltime()
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(xas_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(xas_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    
    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         dft_control=dft_control,&
         scf_control=scf_control,&
         input=input)


    ionode = logger%para_env%source==logger%para_env%mepos
    energy_only = .FALSE.
    dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
    scf_section => section_vals_get_subs_vals(dft_section,"SCF",error=error)
    output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".scfLog",error=error)
    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
            "XAS_TP_SCF WAVEFUNCTION OPTIMIZATION"
    END IF

!   True the GAPW method is used
    gapw = dft_control%qs_control%gapw

    CALL get_qs_env(qs_env=qs_env,&
         matrix_h=matrix_h,&
         matrix_s=matrix_s,energy=energy,&
         particle_set=particle_set,&
         qs_charges=qs_charges,&
         ks_env=ks_env)

    scf_work1 => scf_env%scf_work1
    scf_work2 => scf_env%scf_work2
    ortho => scf_env%ortho

    CALL get_qs_env(qs_env,matrix_ks=matrix_ks,rho=rho,mos=mos)
    matrix_p => rho%rho_ao
    nao = mos(1)%mo_set%nao
    ALLOCATE(vecbuffer(1,nao),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    vecbuffer = 0.0_dp
    ALLOCATE(vecbuffer2(1,xas_control%nexc_search),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    vecbuffer2 = 0.0_dp

    ! do some assertions here on these matrices having the same structure, 
    ! as is currently required
    DO ispin=1,SIZE(matrix_s)
      id_equal=(matrix_s(ispin)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id)
      CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    ENDDO
    DO ispin=1,SIZE(matrix_p)
      id_equal=(matrix_p(ispin)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id)
      CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    ENDDO
    DO ispin=1,SIZE(matrix_h)
      id_equal=(matrix_h(ispin)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id)
      CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    ENDDO
    DO ispin=1,SIZE(matrix_ks)
      id_equal=(matrix_ks(ispin)%matrix%sparsity_id==matrix_s(1)%matrix%sparsity_id)
      CPPrecondition(id_equal,cp_failure_level,routineP,error,failure)
    ENDDO
    ! end sparsity check

    ! Some elements from the xas_env
     CALL get_xas_env(xas_env=xas_env,excvec_coeff=excvec_coeff,nvirtual2=nvirtual2,&
          excvec_overlap=excvec_overlap,fm_work=fm_work,centers_wfn=centers_wfn)
     CPPrecondition(ASSOCIATED(excvec_overlap),cp_failure_level,routineP,error,failure)

    iter_count = 0
    diis_step = .FALSE.
    use_jacobi = .FALSE.

    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,&
            FMT="(/,T3,A,T9,A,T36,A,T49,A,T68,A,/,T3,A)")&
            "Step","Update method","Time","Convergence","Total energy",&
            REPEAT("-",77)
    END IF
    

!   *** SCF loop ***

    scf_loop: DO
       CALL timeset("scf_iter","I"," ",handle2)

       IF (ionode) CALL m_flush(output_unit)

       iter_count = iter_count + 1

! ** here qs_env%rho%rho_r and qs_env%rho%rho_g should be up to date

       CALL qs_ks_update_qs_env(ks_env,qs_env=qs_env,&
            error=error,&
            calculate_forces=.FALSE.,just_energy=energy_only)

       scf_env%p_mix = xas_env%p_mix
       scf_env%iter_param = 0.0_dp
       IF (scf_env%mixing_method==ao_mixing_nr) scf_env%iter_param = scf_env%p_mix

       SELECT CASE (xas_env%scf_method)
         CASE DEFAULT
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                 routineP,"unknown scf method method:"//&
                 cp_to_string(xas_env%scf_method),error,failure)

         CASE(xas_scf_general) ! diagonalisation (default)

            DO ispin=1,dft_control%nspins 
             CALL copy_sm_to_fm(matrix_ks(ispin)%matrix,scf_work1(ispin)%matrix)
            ENDDO

            IF (iter_count > 1) THEN
               CALL qs_diis_b_step(scf_env%scf_diis_buffer,mos,scf_work1,&
                    scf_work2,scf_env%iter_delta,diis_error,diis_step,&
                    xas_env%eps_diis,s_matrix=matrix_s,scf_section=scf_section,&
                    error=error)
            END IF

            do_level_shift = ((scf_control%level_shift /= 0.0_dp).AND.&
               ((scf_control%density_guess == core_guess).OR.(iter_count > 1)))

            IF (diis_step) THEN
               scf_env%p_mix = 1.0_dp
               scf_env%iter_param = diis_error
               IF (use_jacobi) THEN
                  scf_env%iter_method = "DIIS/Jacobi"
               ELSE
                  scf_env%iter_method = "DIIS/Diag."
               END IF
            ELSE
               IF (use_jacobi) THEN
                  scf_env%iter_method = "Mixing/Jacobi"
               ELSE
                  scf_env%iter_method = "Mixing/Diag."
               END IF
            END IF

            IF ((iter_count > 1).AND.(scf_env%iter_delta < scf_control%eps_jacobi)) THEN
               use_jacobi = .TRUE.
            ELSE
               use_cholesky = scf_control%use_cholesky
               use_jacobi = .FALSE.
            END IF

            scf_env%iter_delta=0.0_dp

            DO ispin=1,dft_control%nspins

               CALL eigensolver(scf_work1(ispin)%matrix,&
                    mos(ispin)%mo_set,ortho,scf_work2,&
                    do_level_shift,scf_control%level_shift,&
                    use_cholesky=use_cholesky,&
                    work_syevx=scf_control%work_syevx,&
                    use_jacobi=use_jacobi,&
                    jacobi_threshold=scf_control%jacobi_threshold,&
                    smear=0.0_dp, &
                    error=error)

               IF(ispin == 1) THEN
! ** use the maximum overlap criterion to find the index of the excited orbital
                 CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff,nao=nao,nmo=nmo,&
                      eigenvalues=eigenvalues,homo=homo)
!TEST : Calculate the centers of the wavefunctions to check the localization
                 IF(.FALSE.) THEN
                   CALL get_xas_env(xas_env=xas_env, qs_loc_env=qs_loc_env)
                   CALL get_qs_loc_env(qs_loc_env=qs_loc_env,&
                        localized_wfn_control=localized_wfn_control,&
                        moloc_coeff=moloc_coeff,&
                        op_sm_set=op_sm_set,op_fm_set=op_fm_set,&
                        cell=cell, weights=weights)
 
                   CALL cp_fm_to_fm(mo_coeff,moloc_coeff(1)%matrix,xas_control%nexc_search,1,1 )
                   CALL optimize_loc_berry(do_loc_none, localized_wfn_control, &
                        moloc_coeff(1)%matrix, op_sm_set, &
                        op_fm_set, cell, weights, 1, ionode, error=error)

                   ra(1:3) = particle_set(iatom)%r(1:3)
                   DO istate = 1,xas_control%nexc_search
                      centers_wfn(1,istate) = localized_wfn_control%centers_set(1)%array(1,istate)
                      centers_wfn(2,istate) = localized_wfn_control%centers_set(1)%array(2,istate)
                      centers_wfn(3,istate) = localized_wfn_control%centers_set(1)%array(3,istate)
                      rc(1:3) = centers_wfn(1:3,istate)
                      rac = pbc(ra,rc,cell)
                   END DO 
                 END IF
!TEST

                 CALL cp_sm_fm_multiply(matrix_s(1)%matrix,mo_coeff,fm_work,ncol=nmo)
                 CALL cp_fm_gemm("T","N",1,xas_control%nexc_search,nao,1.0_dp,excvec_coeff,&
                      fm_work,0.0_dp,excvec_overlap,b_first_col=1,error=error)
                 CALL cp_fm_get_info ( matrix=excvec_overlap, col_indices = col_indices,&
                        nrow_global = m, ncol_global = n )
!        CALL cp_fm_get_info ( matrix=excvec_coeff, nrow_global = m, ncol_global = n )
!        CALL cp_fm_get_info ( matrix=fm_work, nrow_global = m, ncol_global = n )
                 CALL cp_fm_get_submatrix(excvec_overlap,vecbuffer2,1,1,&
                      1,xas_control%nexc_search,transpose=.FALSE.,error=error)

                 uniform_occupation = .FALSE.
                 CALL get_mo_set(mos(1)%mo_set,&
                        occupation_numbers=occupation_numbers)

                 b_max = 0.0_dp 

                 DO i = 1,xas_control%nexc_search
                   a_max = ABS(vecbuffer2(1,i))
                   IF(a_max > b_max)THEN
                      ic_max = i 
                      b_max = a_max
                   ENDIF
                 END DO

                 IF(ic_max /= dft_control% xas_estate) THEN
                   ir_max = dft_control% xas_estate
                   dft_control% xas_estate = ic_max

                   occupation_numbers(dft_control%xas_estate) = xas_control%occ_estate
                   occupation_numbers(ir_max) = 1.0_dp
                   occupation_numbers(homo) =  xas_control%occ_homo
                 END IF

                 ! IOnization Potential
                 IP_energy = eigenvalues(dft_control% xas_estate)

                 CALL set_mo_set(mos(1)%mo_set,&
                      uniform_occupation=uniform_occupation)
                 CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,dft_control%xas_estate,&
                      nao,1,transpose=.TRUE.,error=error)
                 CALL cp_fm_set_submatrix(excvec_coeff,vecbuffer,1,1,&
                      nao,1,transpose=.TRUE.,error=error)
               END IF 

               IF(xas_control%smear/=0.0_dp) THEN
                  CALL xas_smearing(mos(ispin)%mo_set,xas_control,dft_control%xas_estate,ispin,ionode,error=error)
               END IF

               CALL write_mo_set(mos(ispin)%mo_set,atomic_kind_set,particle_set,&
                    4,6,dft_section,globenv)

               CALL calculate_density_matrix(mos(ispin)%mo_set,&
                    scf_env%p_mix_new(ispin)%matrix,error=error)

            ENDDO

         CASE(xas_scf_ot) ! orbital transforms
! in case of LSD the first spin qs_ot_env will drive the minimization
! in the case of a restricted calculation, it will make sure the spin orbitals are equal
            orthogonality_metric=>matrix_s(1)%matrix

            CALL ot_scf_mini(mos,qs_env%mo_derivs,orthogonality_metric, &
                             energy%total, energy_only,scf_env%iter_delta, &
                             scf_env%qs_ot_env)

            DO ispin=1,dft_control%nspins
              IF(ispin == 1) THEN
! ** use the maximum overlap criterion to find the index of the excited orbital
                 CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff,nao=nao,nmo=nmo,&
                      eigenvalues=eigenvalues,homo=homo)

                 CALL cp_sm_fm_multiply(matrix_s(1)%matrix,mo_coeff,fm_work,ncol=nmo)

                 CALL cp_fm_gemm("T","N",1,xas_control%nexc_search,nao,1.0_dp,excvec_coeff,&
                      fm_work,0.0_dp,excvec_overlap,b_first_col=1,error=error)
                 CALL cp_fm_get_info ( matrix=excvec_overlap, col_indices = col_indices)

                 CALL cp_fm_get_submatrix(excvec_overlap,vecbuffer2,1,1,&
                      1,xas_control%nexc_search,transpose=.TRUE.,error=error)

                 uniform_occupation = .FALSE.
                 CALL get_mo_set(mos(1)%mo_set,&
                        occupation_numbers=occupation_numbers)

                 b_max = 0.0_dp
 
                 DO i = 1,xas_control%nexc_search
                   a_max = ABS(vecbuffer2(1,i))
                   IF(a_max > b_max)THEN
                      ic_max = i 
                      b_max = a_max
                   ENDIF
                 END DO

                 IF(ic_max /= dft_control% xas_estate) THEN
                   ir_max = dft_control% xas_estate
                   dft_control% xas_estate = ic_max

                   occupation_numbers(dft_control%xas_estate) = xas_control%occ_estate
                   occupation_numbers(ir_max) = 1.0_dp
                   occupation_numbers(homo) =  xas_control%occ_homo
                 END IF

                 ! IOnization Potential
                 IP_energy = eigenvalues(dft_control% xas_estate)

                 CALL set_mo_set(mos(1)%mo_set,&
                      uniform_occupation=uniform_occupation)
                 CALL cp_fm_get_submatrix(mo_coeff,vecbuffer,1,dft_control%xas_estate,&
                      nao,1,transpose=.TRUE.,error=error)
                 CALL cp_fm_set_submatrix(excvec_coeff,vecbuffer,1,1,&
                      nao,1,transpose=.TRUE.,error=error)
               END IF

               CALL calculate_density_matrix(mos(ispin)%mo_set,&
                     rho%rho_ao(ispin)%matrix,error=error)
            ENDDO

            scf_env%iter_method=scf_env%qs_ot_env(1)%OT_METHOD_FULL
            scf_env%iter_param=scf_env%qs_ot_env(1)%ds_min

       END SELECT


       IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
          CALL scf_env_density_mixing(scf_env%p_mix_new,scf_env%mixing_method,&
                                    scf_env%p_mix,scf_env%iter_delta,&
                                    qs_env=qs_env,error=error)
       ENDIF


       t2 = m_walltime()

       IF (output_unit>0.and.scf_env%print_iter_line) THEN
          WRITE (UNIT=output_unit,&
               FMT="(T2,I5,2X,A,T22,E10.2,T32,F8.2,T40,2F20.10)")&
               iter_count,TRIM(scf_env%iter_method),&
               scf_env%iter_param,t2 - t1,scf_env%iter_delta,energy%total
       END IF

! ** convergence check
       CALL external_control(should_stop,"SCF",globenv,error)
       IF (scf_env%iter_delta < xas_control%eps_scf) THEN
          IF (output_unit>0) THEN
             WRITE(UNIT=output_unit,FMT="(/,T3,A,I5,A/)")&
                  "*** SCF run converged in ",iter_count," steps ***"
          END IF
          CALL timestop(0.0_dp,handle2)
          EXIT scf_loop
       ELSE IF (should_stop.OR.&
            iter_count == xas_control%max_scf) THEN
          IF (output_unit>0) THEN
             WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
                  "*** SCF run NOT converged ***"
          END IF
          CALL timestop(0.0_dp,handle2)
          EXIT scf_loop
       END IF

    IF (.NOT.BTEST(cp_print_key_should_output(logger%iter_info,&
         scf_section,"PRINT%ITERATION_INFO/TIME_CUMUL",error=error),cp_p_file)) t1 = m_walltime()

!   *** mixing methods have the new density matrix in p_mix_new
       IF (scf_env%mixing_method.EQ.ao_mixing_nr) THEN
           DO ispin=1,dft_control%nspins
              CALL transfer_matrix(scf_env%p_mix_new(ispin)%matrix,&
                                       rho%rho_ao(ispin)%matrix)
           END DO
       ENDIF


! ** update qs_env%rho
       CALL qs_rho_update_rho(rho, qs_env=qs_env, gapw=gapw, error=error)
       IF(gapw) THEN
         CALL calculate_rho_atom_coeff(qs_env,rho%rho_ao)
       ENDIF

       CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)      

       CALL timestop(0.0_dp,handle2)

    END DO scf_loop

    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
            "Total electronic density (r-space): ",&
            SUM(rho%tot_rho_r),&
            SUM(rho%tot_rho_r)+ REAL(scf_env%nelectron,dp),&
            "Total core charge density (r-space):",&
            qs_charges%total_rho_core_rspace,&
            qs_charges%total_rho_core_rspace - REAL(scf_env%nelectron+dft_control%charge,dp)
       IF(gapw) THEN
          tot1_h =  qs_charges%total_rho1_hard(1)
          tot1_s =  qs_charges%total_rho1_soft(1)
          DO ispin=2,dft_control%nspins
            tot1_h = tot1_h + qs_charges%total_rho1_hard(ispin)
            tot1_s = tot1_s + qs_charges%total_rho1_soft(ispin)
          END DO
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
                "Hard and soft densities (Lebedev):",&
                tot1_h, tot1_s
          WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
               "Total Rho_soft + Rho1_hard - Rho1_soft (r-space): ",&
               SUM(rho%tot_rho_r)+ tot1_h - tot1_s ,&
               "Total charge density (r-space):      ",&
               SUM(rho%tot_rho_r)+ tot1_h - tot1_s &
               + qs_charges%total_rho_core_rspace,&
               "Total Rho_soft + Rho0_soft (g-space):",&
               qs_charges%total_rho_gspace
       ELSE
         WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
            "Total charge density (r-space):     ",&
            SUM(rho%tot_rho_r)+&
            qs_charges%total_rho_core_rspace,&
            "Total charge density (g-space):     ",qs_charges%total_rho_gspace
       END IF

         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Overlap energy of the core charge distribution:",energy%core_overlap,&
            "Self energy of the core charge distribution:   ",energy%core_self,&
            "Core Hamiltonian energy:                       ",energy%core,&
            "Hartree energy:                                ",energy%hartree,&
            "Exchange-correlation energy:                   ",energy%exc
         IF(gapw) THEN
           WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
              "GAPW| Exc from hard and soft atomic rho1:      ",energy%exc1,&
              "GAPW| local Eh = 1 center integrals:           ",energy%hartree_1c,&
              "GAPW| local Eh = 2 center integrals:           ",energy%hartree_2c,&
              "GAPW| local Eh = 3 centers integrals:          ",energy%hartree_3c
         END IF
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Total energy:                                  ",energy%total
         WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
            "Ionization potential of the excited atom:      ",IP_energy
       CALL m_flush(output_unit)
    END IF

    CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,error=error)

    DEALLOCATE(vecbuffer,vecbuffer2,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF(.NOT.xas_env%scf_method==xas_scf_general) THEN
    ! Rotate the wfn to get the eugenstate of the KS hamiltonian
    ! Only ispin=1 should be needed
      DO ispin=1,dft_control%nspins
        CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff,nao=nao,nmo=nmo,&
             eigenvalues=eigenvalues,homo=homo)
        CALL calculate_subspace_eigenvalues(mo_coeff,&
             matrix_ks(ispin)%matrix,eigenvalues, &
             para_env=qs_env%para_env, &
             do_rotation=.TRUE.)
      END DO  ! ispin
    END IF

    !Calculate the virtual states from the KS matrix matrix_ks(1)
    IF(nvirtual2 .GT. 0) THEN
      NULLIFY(mo_coeff,uno_orbs,uno_evals,preconditioner)
      ALLOCATE(preconditioner,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      CALL init_preconditioner(preconditioner,para_env=qs_env%para_env,&
           blacs_env=qs_env%blacs_env)
      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff,nao=nao,nmo=nmo)
      CALL make_preconditioner(preconditioner, &
                  precon_type=ot_precond_full_kinetic, &
                  matrix_h=matrix_ks(1)%matrix,matrix_s=matrix_s(1)%matrix,&
                  matrix_t=qs_env%kinetic(1)%matrix, &
                  mo_coeff=mo_coeff,energy_gap=0.2_dp)

      CALL get_xas_env(xas_env=xas_env,unoccupied_orbs=uno_orbs,&
           unoccupied_evals=uno_evals,unoccupied_eps=uno_eps,unoccupied_max_iter=uno_iter)
      CALL cp_fm_init_random(uno_orbs,nvirtual2)
      CALL ot_eigensolver(matrix_h=matrix_ks(1)%matrix,matrix_s=matrix_s(1)%matrix, &
           matrix_c=uno_orbs,matrix_orthogonal_space=mo_coeff,&
           preconditioner=preconditioner,eps_gradient=uno_eps,&
           iter_max=uno_iter,globenv=globenv,size_ortho_space=nmo)
      CALL calculate_subspace_eigenvalues(uno_orbs,matrix_ks(1)%matrix,&
           uno_evals,para_env=qs_env%para_env,do_rotation=.TRUE.)
      CALL destroy_preconditioner(preconditioner)
      DEALLOCATE(preconditioner,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    ! Prapare arrays for the calculation of the spectrum
    IF(.NOT. xas_control%xas_method == xas_dscf) THEN
      ! Copy the final vectors in the array
      NULLIFY(all_vectors,all_evals)
      CALL get_xas_env(xas_env=xas_env,all_vectors=all_vectors,&
           all_evals=all_evals,nvirtual1=nvirtual1)
      CALL get_mo_set(mos(1)%mo_set, eigenvalues=eigenvalues,mo_coeff=mo_coeff,&
           nao=nao,nmo=nmo,homo=homo)

      IF(xas_control%emission) THEN
        CALL cp_fm_to_fm(mo_coeff,all_vectors,ncol=nvirtual1,&
             source_start=1,target_start=1)
      ELSE
        CALL cp_fm_to_fm(mo_coeff,all_vectors,ncol=nvirtual1,&
             source_start=homo+1,target_start=1)
      END IF
      IF(nvirtual2 .GT. 0) THEN
        CALL cp_fm_to_fm(uno_orbs,all_vectors,ncol=nvirtual2,&
             source_start=1,target_start=1+nvirtual1)
      END IF

      IF(xas_control%emission) THEN
        DO istate = 1,nvirtual1
          all_evals(istate) = eigenvalues(istate)
        ENDDO
      ELSE
        DO istate = 1,nvirtual1
          all_evals(istate) = eigenvalues(homo+istate)
        ENDDO
      END IF
      DO istate = 1,nvirtual2
        all_evals(istate+nvirtual1) = uno_evals(istate)
      END DO 
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
         "PRINT%PROGRAM_RUN_INFO", error=error)    
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE xas_do_tp_scf

!!*** **********************************************************************

  SUBROUTINE xas_smearing(mo_set,xas_control,estate,ispin,ionode,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(xas_control_type)                   :: xas_control
    INTEGER, INTENT(IN)                      :: estate, ispin
    LOGICAL, INTENT(IN)                      :: ionode
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: homo, imin, imo, lfomo, &
                                                nelectron, nmo, nomo
    REAL(dp)                                 :: e1, e2, edelta, edist, nelec, &
                                                nelec0, occ_estate, smear
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues, occupation

    CALL get_mo_set(mo_set,eigenvalues=eigenvalues,&
         occupation_numbers=occupation,nelectron=nelectron,nmo=nmo)

    nomo = nelectron
    IF(ispin == 1) THEN
      IF(xas_control%xas_method==xas_dscf .OR. xas_control%xas_method==xas_tp_xhh &
         .OR.  xas_control%xas_method==xas_tp_xfh) THEN
         nomo  = nomo + 1
      END IF    
      occ_estate = occupation(estate)
      occupation(1:nomo) = 1.0_dp
      occupation(estate) = occ_estate
      occupation(nomo+1:nmo) = 0.0_dp
      IF(xas_control%xas_method==xas_dscf .OR. xas_control%xas_method==xas_tp_xhh &
         .OR.  xas_control%xas_method==xas_tp_xfh) THEN
        occupation(nomo) = 1.0_dp - occ_estate
        nelec0 = 0.0_dp
        DO imo = 1,nomo
          nelec0 = nelec0 + occupation(imo)
        END DO
      END IF    
    ELSE
      occupation(1:nomo) = 1.0_dp
      occupation(nomo+1:nmo) = 0.0_dp
    END IF
 
    imin = estate+1
    smear = xas_control%smear

    IF(xas_control%fermidist) THEN
       ! not implemented yet
    ELSE
      e1 = eigenvalues(nomo) - 0.5_dp*smear
      e2 = eigenvalues(nomo) + 0.5_dp*smear

      DO imo=imin,nomo
        IF (eigenvalues(imo) > e1) THEN
          lfomo = imo
          EXIT
        END IF
      END DO

      IF(e2 > eigenvalues(nmo)) THEN
        IF(ionode) WRITE(6,*) "WARNING might be too few states for the smearing !"
      END IF
      DO imo=nmo,nomo,-1
        IF (eigenvalues(imo) < e2) THEN
          homo = imo
          EXIT
        END IF
      END DO

!     *** Get the number of electrons to be smeared ***

      edist = 0.0_dp
      nelec = 0.0_dp

      DO imo=lfomo,homo
        nelec = nelec + occupation(imo)
        edist = edist + ABS(e2 - eigenvalues(imo))
      END DO

!     *** Smear electrons inside the window ***

      DO imo=lfomo,homo
        edelta = ABS(e2 - eigenvalues(imo))
        occupation(imo) = MIN(1.0_dp,nelec*edelta/edist)
        nelec = nelec - occupation(imo)
        edist = edist - edelta
      END DO

!     *** Check, if the smearing involves more than one MO ***

      IF (lfomo == homo) THEN
        homo = nomo
        lfomo = nomo + 1
      END IF

    END IF
  
    IF(ispin==1) THEN
      occupation(estate) = occ_estate
      IF(xas_control%xas_method==xas_dscf .OR. xas_control%xas_method==xas_tp_xhh &
           .OR.  xas_control%xas_method==xas_tp_xfh) THEN
        nelec = 0.0_dp
        DO imo = 1,homo
          nelec = nelec + occupation(imo)
        END DO
        IF(nelec/=nelec0 .AND. ionode)&
            WRITE(6,*) "WARNING different number of electrons !", nelec, nelec0
      END IF
    END IF

    CALL set_mo_set(mo_set,homo=homo,lfomo=lfomo,&
                      uniform_occupation=.FALSE.)
    
  END SUBROUTINE xas_smearing
!***************************************************************************

END MODULE xas_tp_scf
